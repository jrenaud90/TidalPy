#pragma once
/*
 * solidliquid_.hpp: c_SolidLiquidLayer, a thermo-mechanical layer with phase changes, built on c_PhysicsLayer.
 *
 * Adds the optional cooling and radiogenics sub-models, and the conductive and adiabatic calculations that need
 * the layer's geometry or solved profile. The thermal constants those read (conductivity, expansivity, heat
 * capacity) belong to the material, the layer's EOS model, like every other material property. All MKS.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::SolidLiquidLayer (102)
 *   payload:
 *     [all c_BaseLayer fields: same byte layout as the BaseLayer binary payload]
 *     [all c_PhysicsLayer additions: love_numbers k/h/l re+im (6×8), the three classification flags,
 *      temperature, use_thermal_eos, use_heating]
 *     eos_model       presence flag (uint8_t, 1) + (if present) its binary record
 *     shear_rheology  presence flag (uint8_t, 1) + (if present) its binary record
 *     bulk_rheology   presence flag (uint8_t, 1) + (if present) its binary record
 *     cooling         presence flag (uint8_t, 1) + (if present) its binary record
 *     radiogenics     presence flag (uint8_t, 1) + (if present) its binary record
 *   The attached material EOS (which carries its own viscosity and partial-melt models), rheology, cooling, and
 *   radiogenics models are serialized recursively: the five presence flags belong to this payload and each
 *   nested model follows as its own record. The EOS profile data is not serialized; re-run the world EOS solve
 *   after loading.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

#include "physics_.hpp"
#include "cooling_.hpp"
#include "radiogenics_.hpp"

namespace tidalpy {

// Construction parameters for c_SolidLiquidLayer. It adds no scalars to c_PhysicsConfig: what distinguishes the
// class is the cooling and radiogenics models it can hold.
struct c_SolidLiquidConfig : public c_PhysicsConfig {};

class c_SolidLiquidLayer : public c_PhysicsLayer {
public:
    c_SolidLiquidLayer() = default;

    explicit c_SolidLiquidLayer(const c_SolidLiquidConfig& cfg)
        : c_PhysicsLayer(cfg)
    {}

    ~c_SolidLiquidLayer() override = default;

    // The unique_ptr members here and in c_PhysicsLayer delete the implicit copy assignment that Cython's stack
    // allocation needs. Cython temporaries always have null sub-model pointers, so resetting on copy is safe.
    c_SolidLiquidLayer& operator=(const c_SolidLiquidLayer& other) noexcept {
        if (this != &other) {
            c_PhysicsLayer::operator=(other);
            this->p_cooling.reset();
            this->p_radiogenics.reset();
        }
        return *this;
    }
    c_SolidLiquidLayer& operator=(c_SolidLiquidLayer&&) noexcept = default;

    // Thermal constants of the material, read from the layer's EOS model (NaN when none is attached).
    double get_thermal_conductivity() const noexcept {
        return this->p_eos ? this->p_eos->get_thermal_conductivity() : TidalPyConstants::d_NAN;
    }
    double get_thermal_expansion() const noexcept {
        return this->p_eos ? this->p_eos->get_thermal_expansion() : TidalPyConstants::d_NAN;
    }
    double get_heat_capacity() const noexcept {
        return this->p_eos ? this->p_eos->get_heat_capacity() : TidalPyConstants::d_NAN;
    }

    uint32_t get_layer_class_id() const noexcept override {
        return static_cast<uint32_t>(BinaryClassID::SolidLiquidLayer);
    }

    // Thermal transport (const, MKS). Each reads the material's constants, which do not vary with temperature or
    // pressure; what the layer adds is its geometry and its solved profile. The temperature and pressure arguments
    // are accepted but currently unused.

    // The material's constant thermal conductivity k [W/(m K)]; the temperature is unused.
    double calc_thermal_conductivity(double /*temperature*/) const noexcept {
        return this->get_thermal_conductivity();
    }

    // kappa = k / (rho c_p)  [m^2/s] from the material's constants at the layer's bulk density (its mass over its
    // volume); the temperature is unused.
    double calc_thermal_diffusivity(double /*temperature*/) const noexcept {
        return this->p_eos ? this->p_eos->calc_thermal_diffusivity(this->get_density_bulk())
                           : TidalPyConstants::d_NAN;
    }

    // Adiabatic temperature gradient [K/m]: dT/dr = alpha T g / c_p, with the material's constant alpha and c_p at
    // the given temperature; the pressure is unused. Gravity comes from the EOS profile at the layer's outer
    // boundary, read under the owning world's call lock; 0.0 when that profile is unpopulated.
    double calc_adiabatic_temperature_gradient(double temperature,
                                               double /*pressure*/) const noexcept {
        double g = 0.0;
        {
            const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
            if (this->p_eos_data.is_populated()) {
                double state[C_EOS_DY_VALUES];
                this->p_eos_data.evaluate(this->p_radius, state);
                g = state[C_EOS_GRAVITY_INDEX];
            }
        }
        if (g <= 0.0 || temperature <= 0.0) { return 0.0; }
        return this->get_thermal_expansion() * temperature * g / this->get_heat_capacity();
    }

    // Conductive heat flux [W/m^2]: q = k (T_base - T_top) / thickness; 0.0 for a zero-thickness layer.
    double calc_heat_flux_conductive(double temperature_base,
                                     double temperature_top) const noexcept {
        if (this->p_thickness <= 0.0) { return 0.0; }
        return this->get_thermal_conductivity()
               * (temperature_base - temperature_top)
               / this->p_thickness;
    }

    // Radiogenic heating [W] from the attached model; 0.0 when none is attached.
    double calc_radiogenic_heating(double time, double mass) const noexcept {
        if (!this->p_radiogenics) { return 0.0; }
        return this->p_radiogenics->calc_heating(time, mass);
    }

    // Sub-model setters (transfer ownership; each registers this layer as the observer).
    void set_cooling(std::unique_ptr<c_CoolingBase> cooling) {
        this->p_cooling = std::move(cooling);
        if (this->p_cooling) { this->p_cooling->set_layer_ptr(this); }
    }

    void set_radiogenics(std::unique_ptr<c_RadiogenicsBase> radiogenics) {
        this->p_radiogenics = std::move(radiogenics);
        if (this->p_radiogenics) { this->p_radiogenics->set_layer_ptr(this); }
    }

    bool get_cooling_set()     const noexcept { return this->p_cooling     != nullptr; }
    bool get_radiogenics_set() const noexcept { return this->p_radiogenics != nullptr; }

    // Non-owning observer pointers (nullptr if unset).
    c_CoolingBase*     get_cooling_model()     const noexcept { return this->p_cooling.get(); }
    c_RadiogenicsBase* get_radiogenics_model() const noexcept { return this->p_radiogenics.get(); }

    void write_binary(std::ostream& out) const override {
        write_binary_header(
            out, static_cast<uint32_t>(BinaryClassID::SolidLiquidLayer),
            this->p_base_fields_bytes() + this->p_physics_fields_bytes()
                + optional_binary_flag_bytes()             // material EOS model presence flag
                + this->physics_models_presence_bytes()    // shear and bulk rheology presence flags
                + 2 * optional_binary_flag_bytes());       // cooling and radiogenics presence flags
        this->p_write_base_fields(out);
        this->p_write_physics_fields(out);
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write SolidLiquidLayer binary data");
        }
        this->write_eos_model_binary(out);
        this->write_physics_models_binary(out);
        this->write_submodels_binary(out);
    }

    void read_binary(std::istream& in, bool force = false) override {
        this->p_begin_read_binary(in, force);
        this->p_read_base_fields(in);
        this->p_read_physics_fields(in);
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read SolidLiquidLayer binary data");
        }
        this->read_eos_model_binary(in, force);
        this->read_physics_models_binary(in, force);
        this->read_submodels_binary(in, force);
        this->update_physicals();
    }

protected:
    // Recursive (de)serialization of the optional cooling and radiogenics models, mirroring
    // c_PhysicsLayer::write_physics_models_binary: a presence flag each, followed when set by the model's own
    // record. On read the concrete model is rebuilt through the cooling and radiogenics binary-dispatch
    // factories and re-registered as this layer's observer.
    void write_submodels_binary(std::ostream& out) const {
        write_optional_binary(out, this->p_cooling);
        write_optional_binary(out, this->p_radiogenics);
    }

    void read_submodels_binary(std::istream& in, bool force) {
        this->p_cooling =
            read_optional_binary<c_CoolingBase>(in, force, c_cooling_from_binary);
        if (this->p_cooling) { this->p_cooling->set_layer_ptr(this); }
        this->p_radiogenics =
            read_optional_binary<c_RadiogenicsBase>(in, force, c_radiogenics_from_binary);
        if (this->p_radiogenics) { this->p_radiogenics->set_layer_ptr(this); }
    }

    std::unique_ptr<c_CoolingBase>      p_cooling;
    std::unique_ptr<c_RadiogenicsBase>  p_radiogenics;
};

} // namespace tidalpy
