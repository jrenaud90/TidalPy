#pragma once
/*
 * solidliquid_.hpp: c_SolidLiquidLayer, a thermo-mechanical layer with phase changes, built on c_BaseLayer.
 *
 * Adds the optional cooling and radiogenics sub-models, and the conductive and adiabatic calculations that need
 * the layer's geometry or solved profile. The thermal constants those read (conductivity, expansivity, heat
 * capacity) belong to the material, the layer's EOS model, like every other material property. All MKS.
 *
 * Binary payload: the c_BaseLayer payload, then the cooling and radiogenics models, each behind a presence flag.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

#include "base_.hpp"
#include "cooling_.hpp"
#include "radiogenics_.hpp"

namespace tidalpy {

// Construction parameters for c_SolidLiquidLayer. It adds no scalars to c_BaseLayerConfig: what distinguishes the
// class is the cooling and radiogenics models it can hold.
struct c_SolidLiquidConfig : public c_BaseLayerConfig {};

class c_SolidLiquidLayer : public c_BaseLayer {
public:
    c_SolidLiquidLayer() = default;

    explicit c_SolidLiquidLayer(const c_SolidLiquidConfig& cfg)
        : c_BaseLayer(cfg)
    {}

    ~c_SolidLiquidLayer() override = default;

    // The unique_ptr members here and in c_BaseLayer delete the implicit copy assignment that Cython's stack
    // allocation needs. Cython temporaries always have null sub-model pointers, so resetting on copy is safe.
    c_SolidLiquidLayer& operator=(const c_SolidLiquidLayer& other) noexcept {
        if (this != &other) {
            c_BaseLayer::operator=(other);
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

    uint32_t get_binary_class_id() const override {
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

    // Sub-model setters (transfer ownership; each registers this layer as the observer). A thermal EOS solve reads
    // both, so each makes the owning world forget its solved structure (c_LayerOwner).
    void set_cooling(std::unique_ptr<c_CoolingBase> cooling) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_cooling = std::move(cooling);
        if (this->p_cooling) { this->p_cooling->set_layer_ptr(this); }
        this->p_update_owner_after_change();
    }

    void set_radiogenics(std::unique_ptr<c_RadiogenicsBase> radiogenics) {
        const c_WorldCallLock call_lock(this->p_owner_call_mutex.get());
        this->p_radiogenics = std::move(radiogenics);
        if (this->p_radiogenics) { this->p_radiogenics->set_layer_ptr(this); }
        this->p_update_owner_after_change();
    }

    bool get_cooling_set()     const noexcept { return this->p_cooling     != nullptr; }
    bool get_radiogenics_set() const noexcept { return this->p_radiogenics != nullptr; }

    // Non-owning observer pointers (nullptr if unset).
    c_CoolingBase*     get_cooling_model()     const noexcept { return this->p_cooling.get(); }
    c_RadiogenicsBase* get_radiogenics_model() const noexcept { return this->p_radiogenics.get(); }

    // What load_binary reads a file into first (c_TidalPyBaseClass::make_binary_scratch).
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_SolidLiquidLayer>();
    }

protected:
    void p_write_payload(std::ostream& out) const override {
        c_BaseLayer::p_write_payload(out);
        write_optional_binary(out, this->p_cooling);
        write_optional_binary(out, this->p_radiogenics);
    }

    void p_read_payload(std::istream& in, bool force) override {
        c_BaseLayer::p_read_payload(in, force);
        this->p_cooling = read_optional_binary<c_CoolingBase>(in, force, c_cooling_from_binary);
        if (this->p_cooling) { this->p_cooling->set_layer_ptr(this); }
        this->p_radiogenics = read_optional_binary<c_RadiogenicsBase>(in, force, c_radiogenics_from_binary);
        if (this->p_radiogenics) { this->p_radiogenics->set_layer_ptr(this); }
    }

    std::unique_ptr<c_CoolingBase>      p_cooling;
    std::unique_ptr<c_RadiogenicsBase>  p_radiogenics;
};

} // namespace tidalpy
