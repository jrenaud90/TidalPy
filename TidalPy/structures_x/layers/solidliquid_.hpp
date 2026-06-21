#pragma once
/*
 * solidliquid_.hpp — c_SolidLiquidLayer: thermo-mechanical layer with phase changes.
 *
 * Inherits c_PhysicsLayer (structures_x/layers/physics_.hpp).
 * Adds thermal properties, Arrhenius viscosity, melt-fraction tracking,
 * and optional cooling / radiogenics sub-models.
 *
 * When no cooling or radiogenics object is attached the corresponding
 * calc_* methods return 0.0.
 *
 * All spatial fields are in MKS units.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::SolidLiquidLayer (102)
 *   payload:
 *     [all c_BaseLayer fields — same byte layout as BaseLayer binary payload]
 *     [all c_PhysicsLayer additions — shear modulus, bulk modulus,
 *      shear viscosity, bulk viscosity, love_numbers k/h/l re+im (10×8)]
 *     thermal_conductivity_ref_w_mk  (double, 8)
 *     thermal_expansion_ref_1_k      (double, 8)
 *     heat_capacity_ref_j_kgk        (double, 8)
 *     activation_energy_j_mol        (double, 8)
 *     activation_volume_m3_mol       (double, 8)
 *     solidus_temperature_k          (double, 8)
 *     liquidus_temperature_k         (double, 8)
 *     melt_fraction_exponent         (double, 8)
 *     reference_density_kg_m3        (double, 8)
 *     reference_temperature_k        (double, 8)
 *     melt_viscosity_reduction       (double, 8)
 *     shear_rheology  presence flag (uint8_t, 1) + (if present) its binary record
 *     bulk_rheology   presence flag (uint8_t, 1) + (if present) its binary record
 *     cooling         presence flag (uint8_t, 1) + (if present) its binary record
 *     radiogenics     presence flag (uint8_t, 1) + (if present) its binary record
 *   Attached rheology, cooling, and radiogenics models ARE serialized recursively
 *   (presence flag + the model's own binary record); the four presence flags are
 *   part of this payload, each nested model record follows as a separate record.
 *   EOS profile data is NOT serialized.
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

// -------------------------------------------------------------------------------
// c_SolidLiquidConfig — construction parameters for c_SolidLiquidLayer.
// Extends c_PhysicsConfig with thermal and melt-fraction parameters.
// -------------------------------------------------------------------------------
struct c_SolidLiquidConfig : public c_PhysicsConfig {
    double thermal_conductivity_ref_w_mk = 4.0;       // [W/m/K]
    double thermal_expansion_ref_1_k     = 3.0e-5;    // [1/K]
    double heat_capacity_ref_j_kgk       = 1200.0;    // [J/(kg·K)]
    double activation_energy_j_mol       = 300.0e3;   // [J/mol]
    double activation_volume_m3_mol      = 5.0e-6;    // [m³/mol]
    double solidus_temperature_k         = 1600.0;    // [K]
    double liquidus_temperature_k        = 2000.0;    // [K]
    double melt_fraction_exponent        = 1.0;       // [dimensionless]
    double reference_density_kg_m3       = 3500.0;    // [kg/m³]
    double reference_temperature_k       = 1600.0;    // [K] Arrhenius reference
    double melt_viscosity_reduction      = 25.0;      // [dimensionless] exp coefficient
};

// -------------------------------------------------------------------------------
// c_SolidLiquidLayer
// -------------------------------------------------------------------------------
class c_SolidLiquidLayer : public c_PhysicsLayer {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_SolidLiquidLayer() = default;

    explicit c_SolidLiquidLayer(const c_SolidLiquidConfig& cfg)
        : c_PhysicsLayer(cfg),
          p_thermal_conductivity_ref(cfg.thermal_conductivity_ref_w_mk),
          p_thermal_expansion_ref(cfg.thermal_expansion_ref_1_k),
          p_heat_capacity_ref(cfg.heat_capacity_ref_j_kgk),
          p_activation_energy(cfg.activation_energy_j_mol),
          p_activation_volume(cfg.activation_volume_m3_mol),
          p_solidus_temperature(cfg.solidus_temperature_k),
          p_liquidus_temperature(cfg.liquidus_temperature_k),
          p_melt_fraction_exponent(cfg.melt_fraction_exponent),
          p_reference_density(cfg.reference_density_kg_m3),
          p_reference_temperature(cfg.reference_temperature_k),
          p_melt_viscosity_reduction(cfg.melt_viscosity_reduction)
    {}

    ~c_SolidLiquidLayer() override = default;

    // unique_ptr members in this class and in c_PhysicsLayer (via parent) delete
    // implicit copy-assignment. Cython stack allocation requires it; Cython
    // temporaries always have null sub-model pointers, so resetting is safe.
    c_SolidLiquidLayer& operator=(const c_SolidLiquidLayer& other) noexcept {
        if (this != &other) {
            c_PhysicsLayer::operator=(other);
            this->p_thermal_conductivity_ref  = other.p_thermal_conductivity_ref;
            this->p_thermal_expansion_ref     = other.p_thermal_expansion_ref;
            this->p_heat_capacity_ref         = other.p_heat_capacity_ref;
            this->p_activation_energy         = other.p_activation_energy;
            this->p_activation_volume         = other.p_activation_volume;
            this->p_solidus_temperature       = other.p_solidus_temperature;
            this->p_liquidus_temperature      = other.p_liquidus_temperature;
            this->p_melt_fraction_exponent    = other.p_melt_fraction_exponent;
            this->p_reference_density         = other.p_reference_density;
            this->p_reference_temperature     = other.p_reference_temperature;
            this->p_melt_viscosity_reduction  = other.p_melt_viscosity_reduction;
            this->p_cooling.reset();
            this->p_radiogenics.reset();
        }
        return *this;
    }
    c_SolidLiquidLayer& operator=(c_SolidLiquidLayer&&) noexcept = default;

    // -----------------------------------------------------------------------
    // Thermal property getters (const, MKS)
    // -----------------------------------------------------------------------
    double get_thermal_conductivity_ref()  const noexcept { return this->p_thermal_conductivity_ref; }
    double get_thermal_expansion_ref()     const noexcept { return this->p_thermal_expansion_ref; }
    double get_heat_capacity_ref()         const noexcept { return this->p_heat_capacity_ref; }
    double get_activation_energy()         const noexcept { return this->p_activation_energy; }
    double get_activation_volume()         const noexcept { return this->p_activation_volume; }
    double get_solidus_temperature()       const noexcept { return this->p_solidus_temperature; }
    double get_liquidus_temperature()      const noexcept { return this->p_liquidus_temperature; }
    double get_melt_fraction_exponent()    const noexcept { return this->p_melt_fraction_exponent; }
    double get_reference_density()         const noexcept { return this->p_reference_density; }
    double get_reference_temperature()     const noexcept { return this->p_reference_temperature; }
    double get_melt_viscosity_reduction()  const noexcept { return this->p_melt_viscosity_reduction; }

    uint32_t get_layer_class_id() const noexcept override {
        return static_cast<uint32_t>(BinaryClassID::SolidLiquidLayer);
    }

    // -----------------------------------------------------------------------
    // calc_melt_fraction [dimensionless, 0..1]
    //
    // Linear interpolation between solidus and liquidus, raised to
    // melt_fraction_exponent.  Pressure dependence on the melt curve is
    // deferred to Phase 9 (EOSHandler integration).
    // -----------------------------------------------------------------------
    double calc_melt_fraction(double temperature_k, double /*pressure_pa*/) const noexcept {
        const double dT = this->p_liquidus_temperature - this->p_solidus_temperature;
        if (dT <= 0.0) { return temperature_k >= this->p_solidus_temperature ? 1.0 : 0.0; }
        const double tau = (temperature_k - this->p_solidus_temperature) / dT;
        const double tau_clamped = std::max(0.0, std::min(1.0, tau));
        return std::pow(tau_clamped, this->p_melt_fraction_exponent);
    }

    // -----------------------------------------------------------------------
    // calc_viscosity [Pa·s]
    //
    // Arrhenius temperature and pressure dependence relative to the reference
    // viscosity (p_viscosity_static_pas) at (p_reference_temperature, P=0):
    //
    //   η(T,P) = η_ref * exp((E_a + P·V_a)/(R·T) − E_a/(R·T_ref))
    //
    // Partial-melt exponential suppression (Roscoe-type):
    //
    //   η_eff = η(T,P) * exp(−C·φ)
    //
    // The exponent is clamped to [−100, 100] to prevent overflow/underflow.
    // -----------------------------------------------------------------------
    double calc_viscosity(double temperature_k, double pressure_pa) const noexcept {
        if (temperature_k <= 0.0 || tidalpy_config_ptr == nullptr) { return this->p_shear_viscosity_static_pas; }
        const double R = tidalpy_config_ptr->d_R;
        const double exponent =
            (this->p_activation_energy + pressure_pa * this->p_activation_volume)
                / (R * temperature_k)
            - this->p_activation_energy
                / (R * this->p_reference_temperature);
        const double eta = this->p_shear_viscosity_static_pas
                           * std::exp(std::clamp(exponent, -100.0, 100.0));
        const double phi = calc_melt_fraction(temperature_k, pressure_pa);
        return eta * std::exp(std::clamp(-this->p_melt_viscosity_reduction * phi, -100.0, 0.0));
    }

    // -----------------------------------------------------------------------
    // calc_shear_modulus [Pa]
    //
    // Melt-fraction-reduced shear modulus: G = G_static * (1 − φ).
    // Returns 0 when fully molten.
    // -----------------------------------------------------------------------
    double calc_shear_modulus(double temperature_k, double pressure_pa) const noexcept {
        const double phi = calc_melt_fraction(temperature_k, pressure_pa);
        return this->p_shear_modulus_static_pa * (1.0 - phi);
    }

    // -----------------------------------------------------------------------
    // Thermal transport (const, MKS)
    // -----------------------------------------------------------------------
    double calc_thermal_conductivity(double /*temperature_k*/) const noexcept {
        return this->p_thermal_conductivity_ref;
    }

    // κ = k / (ρ_ref · c_p)  [m²/s]
    double calc_thermal_diffusivity(double temperature_k) const noexcept {
        const double k   = calc_thermal_conductivity(temperature_k);
        const double rho = (this->p_reference_density > 0.0) ? this->p_reference_density : 1.0;
        const double cp  = (this->p_heat_capacity_ref > 0.0)  ? this->p_heat_capacity_ref  : 1.0;
        return k / (rho * cp);
    }

    // -----------------------------------------------------------------------
    // calc_adiabatic_temperature_gradient [K/m]
    //
    // dT/dr_adiabatic = α · T · g / c_p
    //
    // Gravity is taken from the EOS profile at the layer's outer boundary
    // if available; otherwise 0.0 is returned.
    // -----------------------------------------------------------------------
    double calc_adiabatic_temperature_gradient(double temperature_k,
                                               double /*pressure_pa*/) const noexcept {
        double g = 0.0;
        if (this->p_eos_data.is_populated()) {
            g = this->p_eos_data.get_gravity(this->p_radius);
        }
        if (g <= 0.0 || temperature_k <= 0.0) { return 0.0; }
        return this->p_thermal_expansion_ref * temperature_k * g / this->p_heat_capacity_ref;
    }

    // -----------------------------------------------------------------------
    // calc_heat_flux_conductive [W/m²]
    //
    // q = k · (T_base − T_top) / thickness
    // Returns 0.0 when layer has zero thickness.
    // -----------------------------------------------------------------------
    double calc_heat_flux_conductive(double temperature_base_k,
                                     double temperature_top_k) const noexcept {
        if (this->p_thickness <= 0.0) { return 0.0; }
        return this->p_thermal_conductivity_ref
               * (temperature_base_k - temperature_top_k)
               / this->p_thickness;
    }

    // -----------------------------------------------------------------------
    // calc_radiogenic_heating [W]
    //
    // Delegates to p_radiogenics if set; otherwise returns 0.0.
    // -----------------------------------------------------------------------
    double calc_radiogenic_heating(double time_s, double mass_kg) const noexcept {
        if (!this->p_radiogenics) { return 0.0; }
        return this->p_radiogenics->calc_heating(time_s, mass_kg);
    }

    // -----------------------------------------------------------------------
    // Sub-model setters (non-const; transfer ownership via unique_ptr)
    // -----------------------------------------------------------------------
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

    // -----------------------------------------------------------------------
    // Binary I/O
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        const auto     name_len = static_cast<uint32_t>(this->p_name.size());
        const auto     mat_len  = static_cast<uint32_t>(this->p_material_name.size());
        const uint64_t payload  =
            sizeof(double)   * 2 +           // p_radius, p_mass
            sizeof(uint32_t) + name_len +    // name length + bytes
            sizeof(int32_t)  +               // layer_index
            sizeof(double)   +               // radius_inner_m
            sizeof(uint32_t) + mat_len +     // material_name length + bytes
            sizeof(uint8_t)  +               // is_tidal
            sizeof(double)   +               // tidal_scale
            sizeof(double)   * 10 +          // shear/bulk modulus, shear/bulk viscosity, love_numbers k/h/l re+im
            sizeof(double)   * 11 +          // SolidLiquidLayer thermal fields
            this->rheology_presence_bytes() +    // shear + bulk rheology presence flags
            2 * optional_binary_flag_bytes();    // cooling + radiogenics presence flags

        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::SolidLiquidLayer), payload);

        // c_BaseLayer fields
        out.write(reinterpret_cast<const char*>(&this->p_radius), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_mass),   sizeof(double));
        out.write(reinterpret_cast<const char*>(&name_len),       sizeof(uint32_t));
        if (name_len > 0) { out.write(this->p_name.data(), name_len); }
        const int32_t idx = static_cast<int32_t>(this->p_layer_index);
        out.write(reinterpret_cast<const char*>(&idx),                  sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&this->p_radius_inner), sizeof(double));
        out.write(reinterpret_cast<const char*>(&mat_len),              sizeof(uint32_t));
        if (mat_len > 0) { out.write(this->p_material_name.data(), mat_len); }
        const uint8_t is_tidal_byte = static_cast<uint8_t>(this->p_is_tidal);
        out.write(reinterpret_cast<const char*>(&is_tidal_byte),        sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&this->p_tidal_scale),  sizeof(double));

        // c_PhysicsLayer fields
        out.write(reinterpret_cast<const char*>(&this->p_shear_modulus_static_pa),    sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_bulk_modulus_static_pa),     sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_shear_viscosity_static_pas), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_bulk_viscosity_static_pas),  sizeof(double));
        auto write_complex = [&](const std::complex<double>& c) {
            const double re = c.real(), im = c.imag();
            out.write(reinterpret_cast<const char*>(&re), sizeof(double));
            out.write(reinterpret_cast<const char*>(&im), sizeof(double));
        };
        write_complex(this->p_love_numbers.k);
        write_complex(this->p_love_numbers.h);
        write_complex(this->p_love_numbers.l);

        // c_SolidLiquidLayer fields
        out.write(reinterpret_cast<const char*>(&this->p_thermal_conductivity_ref), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_thermal_expansion_ref),    sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_heat_capacity_ref),        sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_activation_energy),        sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_activation_volume),        sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_solidus_temperature),      sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_liquidus_temperature),     sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_melt_fraction_exponent),   sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_reference_density),        sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_reference_temperature),    sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_melt_viscosity_reduction), sizeof(double));

        if (!out) {
            throw std::runtime_error("TidalPy: failed to write SolidLiquidLayer binary data");
        }

        // Attached sub-models (presence flag + recursive record each).
        this->write_rheology_binary(out);    // inherited from c_PhysicsLayer
        this->write_submodels_binary(out);   // cooling + radiogenics
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);

        // c_BaseLayer fields
        in.read(reinterpret_cast<char*>(&this->p_radius), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_mass),   sizeof(double));

        uint32_t name_len = 0;
        in.read(reinterpret_cast<char*>(&name_len), sizeof(uint32_t));
        this->p_name.resize(name_len);
        if (name_len > 0) { in.read(this->p_name.data(), name_len); }

        int32_t idx = 0;
        in.read(reinterpret_cast<char*>(&idx), sizeof(int32_t));
        this->p_layer_index = static_cast<int>(idx);

        in.read(reinterpret_cast<char*>(&this->p_radius_inner), sizeof(double));

        uint32_t mat_len = 0;
        in.read(reinterpret_cast<char*>(&mat_len), sizeof(uint32_t));
        this->p_material_name.resize(mat_len);
        if (mat_len > 0) { in.read(this->p_material_name.data(), mat_len); }

        uint8_t is_tidal_byte = 0;
        in.read(reinterpret_cast<char*>(&is_tidal_byte), sizeof(uint8_t));
        this->p_is_tidal = static_cast<bool>(is_tidal_byte);

        in.read(reinterpret_cast<char*>(&this->p_tidal_scale), sizeof(double));

        // c_PhysicsLayer fields
        in.read(reinterpret_cast<char*>(&this->p_shear_modulus_static_pa),    sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_bulk_modulus_static_pa),     sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_shear_viscosity_static_pas), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_bulk_viscosity_static_pas),  sizeof(double));
        auto read_complex = [&](std::complex<double>& c) {
            double re = 0.0, im = 0.0;
            in.read(reinterpret_cast<char*>(&re), sizeof(double));
            in.read(reinterpret_cast<char*>(&im), sizeof(double));
            c = std::complex<double>(re, im);
        };
        read_complex(this->p_love_numbers.k);
        read_complex(this->p_love_numbers.h);
        read_complex(this->p_love_numbers.l);

        // c_SolidLiquidLayer fields
        in.read(reinterpret_cast<char*>(&this->p_thermal_conductivity_ref), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_thermal_expansion_ref),    sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_heat_capacity_ref),        sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_activation_energy),        sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_activation_volume),        sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_solidus_temperature),      sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_liquidus_temperature),     sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_melt_fraction_exponent),   sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_reference_density),        sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_reference_temperature),    sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_melt_viscosity_reduction), sizeof(double));

        if (!in) {
            throw std::runtime_error("TidalPy: failed to read SolidLiquidLayer binary data");
        }

        // Attached sub-models (presence flag + recursive record each).
        this->read_rheology_binary(in, force);    // inherited from c_PhysicsLayer
        this->read_submodels_binary(in, force);   // cooling + radiogenics

        this->update_physicals();
    }

protected:
    // -----------------------------------------------------------------------
    // Recursive (de)serialization of the optional cooling/radiogenics models.
    // Mirrors c_PhysicsLayer::write_rheology_binary/read_rheology_binary.  Each
    // model is written as a presence flag followed, when set, by the model's own
    // binary record; on read the correct concrete model is rebuilt via the
    // cooling/radiogenics binary-dispatch factories and re-registered as this
    // layer's observer.
    // -----------------------------------------------------------------------
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

    double p_thermal_conductivity_ref = 4.0;       // [W/m/K]
    double p_thermal_expansion_ref    = 3.0e-5;    // [1/K]
    double p_heat_capacity_ref        = 1200.0;    // [J/(kg·K)]
    double p_activation_energy        = 300.0e3;   // [J/mol]
    double p_activation_volume        = 5.0e-6;    // [m³/mol]
    double p_solidus_temperature      = 1600.0;    // [K]
    double p_liquidus_temperature     = 2000.0;    // [K]
    double p_melt_fraction_exponent   = 1.0;       // [dimensionless]
    double p_reference_density        = 3500.0;    // [kg/m³]
    double p_reference_temperature    = 1600.0;    // [K]
    double p_melt_viscosity_reduction = 25.0;      // [dimensionless]

    std::unique_ptr<c_CoolingBase>      p_cooling;
    std::unique_ptr<c_RadiogenicsBase>  p_radiogenics;
};

} // namespace tidalpy
