#pragma once
/*
 * gas_.hpp — c_GasLayer: ideal-gas/fluid layer with reduced feature set.
 *
 * Inherits c_PhysicsLayer (structures_x/layers/physics_.hpp).
 * Adds ideal-gas thermodynamic properties: adiabatic lapse rate, scale height,
 * pressure, and sound speed.  No phase changes, no solidus/liquidus, no
 * cooling or radiogenics sub-models.
 *
 * All spatial fields are in MKS units.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::GasLayer (103)
 *   payload:
 *     [all c_BaseLayer fields — same byte layout as BaseLayer binary payload]
 *     [all c_PhysicsLayer additions — shear modulus, bulk modulus,
 *      shear viscosity, bulk viscosity, love_numbers k/h/l re+im (10×8)]
 *     mean_molecular_weight_kg_mol  (double, 8)
 *     adiabatic_index               (double, 8)
 *     reference_temperature_k       (double, 8)
 *     reference_density_kg_m3       (double, 8)
 *     shear_rheology presence flag (uint8_t, 1) + (if present) its binary record
 *     bulk_rheology  presence flag (uint8_t, 1) + (if present) its binary record
 *   Attached rheology objects (inherited from c_PhysicsLayer) ARE serialized
 *   recursively; the two presence flags are part of this payload, each nested
 *   model record follows as a separate record.
 *   EOS profile data is NOT serialized.
 */

#include <cmath>
#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

#include "physics_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_GasConfig — construction parameters for c_GasLayer.
// Extends c_PhysicsConfig with ideal-gas thermodynamic fields.
// -------------------------------------------------------------------------------
struct c_GasConfig : public c_PhysicsConfig {
    double mean_molecular_weight_kg_mol = 2.0e-3;    // [kg/mol] hydrogen default
    double adiabatic_index              = 1.4;       // γ = c_p/c_v [dimensionless]
    double reference_temperature_k      = 300.0;     // [K]
    double reference_density_kg_m3      = 1.0;       // [kg/m³]
};

// -------------------------------------------------------------------------------
// c_GasLayer
// -------------------------------------------------------------------------------
class c_GasLayer : public c_PhysicsLayer {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_GasLayer() = default;

    explicit c_GasLayer(const c_GasConfig& cfg)
        : c_PhysicsLayer(cfg),
          p_mean_molecular_weight(cfg.mean_molecular_weight_kg_mol),
          p_adiabatic_index(cfg.adiabatic_index),
          p_reference_temperature(cfg.reference_temperature_k),
          p_reference_density(cfg.reference_density_kg_m3)
    {}

    ~c_GasLayer() override = default;

    // unique_ptr members inherited from c_PhysicsLayer delete implicit
    // copy-assignment; define it explicitly for Cython stack allocation.
    // Cython temporaries always have null rheology pointers, so reset() is safe.
    c_GasLayer& operator=(const c_GasLayer& other) noexcept {
        if (this != &other) {
            c_PhysicsLayer::operator=(other);
            this->p_mean_molecular_weight = other.p_mean_molecular_weight;
            this->p_adiabatic_index       = other.p_adiabatic_index;
            this->p_reference_temperature = other.p_reference_temperature;
            this->p_reference_density     = other.p_reference_density;
        }
        return *this;
    }
    c_GasLayer& operator=(c_GasLayer&&) noexcept = default;

    // -----------------------------------------------------------------------
    // Property getters (const, MKS)
    // -----------------------------------------------------------------------
    uint32_t get_layer_class_id() const noexcept override {
        return static_cast<uint32_t>(BinaryClassID::GasLayer);
    }
    double get_mean_molecular_weight() const noexcept { return this->p_mean_molecular_weight; }
    double get_adiabatic_index()       const noexcept { return this->p_adiabatic_index; }
    double get_reference_temperature() const noexcept { return this->p_reference_temperature; }
    double get_reference_density()     const noexcept { return this->p_reference_density; }

    // -----------------------------------------------------------------------
    // calc_adiabatic_lapse_rate [K/m]
    //
    // Dry adiabatic lapse rate for an ideal gas:
    //   Γ = g * (γ - 1) * M / (γ * R)
    //
    // where γ is the adiabatic index, M is mean molecular weight [kg/mol],
    // and R is the universal gas constant [J/(mol·K)].
    // Returns 0.0 when the config pointer is unavailable or inputs are invalid.
    // -----------------------------------------------------------------------
    double calc_adiabatic_lapse_rate(double gravity_m_s2) const noexcept {
        if (gravity_m_s2 <= 0.0 || tidalpy_config_ptr == nullptr) { return 0.0; }
        const double R = tidalpy_config_ptr->d_R;
        if (R <= 0.0 || this->p_adiabatic_index <= 1.0) { return 0.0; }
        return gravity_m_s2 * (this->p_adiabatic_index - 1.0) * this->p_mean_molecular_weight
               / (this->p_adiabatic_index * R);
    }

    // -----------------------------------------------------------------------
    // calc_scale_height [m]
    //
    // Barometric (pressure) scale height:
    //   H = R * T / (g * M)
    //
    // Returns 0.0 when inputs are non-positive or config is unavailable.
    // -----------------------------------------------------------------------
    double calc_scale_height(double temperature_k, double gravity_m_s2) const noexcept {
        if (temperature_k <= 0.0 || gravity_m_s2 <= 0.0
                || this->p_mean_molecular_weight <= 0.0
                || tidalpy_config_ptr == nullptr) {
            return 0.0;
        }
        const double R = tidalpy_config_ptr->d_R;
        return R * temperature_k / (gravity_m_s2 * this->p_mean_molecular_weight);
    }

    // -----------------------------------------------------------------------
    // calc_pressure_ideal_gas [Pa]
    //
    // Ideal gas law: P = ρ * R * T / M
    //
    // Returns 0.0 when inputs are non-positive or config is unavailable.
    // -----------------------------------------------------------------------
    double calc_pressure_ideal_gas(double temperature_k,
                                   double density_kg_m3) const noexcept {
        if (temperature_k <= 0.0 || density_kg_m3 <= 0.0
                || this->p_mean_molecular_weight <= 0.0
                || tidalpy_config_ptr == nullptr) {
            return 0.0;
        }
        const double R = tidalpy_config_ptr->d_R;
        return density_kg_m3 * R * temperature_k / this->p_mean_molecular_weight;
    }

    // -----------------------------------------------------------------------
    // calc_sound_speed [m/s]
    //
    // Adiabatic sound speed for an ideal gas:
    //   c_s = sqrt(γ * R * T / M)
    //
    // Returns 0.0 when inputs are invalid or config is unavailable.
    // -----------------------------------------------------------------------
    double calc_sound_speed(double temperature_k) const noexcept {
        if (temperature_k <= 0.0 || this->p_mean_molecular_weight <= 0.0
                || tidalpy_config_ptr == nullptr) {
            return 0.0;
        }
        const double R = tidalpy_config_ptr->d_R;
        return std::sqrt(this->p_adiabatic_index * R * temperature_k
                         / this->p_mean_molecular_weight);
    }

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
            sizeof(double)   * 4 +           // GasLayer fields
            this->rheology_presence_bytes(); // shear + bulk rheology presence flags

        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::GasLayer), payload);

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

        // c_GasLayer fields
        out.write(reinterpret_cast<const char*>(&this->p_mean_molecular_weight), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_adiabatic_index),       sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_reference_temperature), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_reference_density),     sizeof(double));

        if (!out) {
            throw std::runtime_error("TidalPy: failed to write GasLayer binary data");
        }

        // Attached rheology models (presence flag + recursive record each).
        this->write_rheology_binary(out);
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

        // c_GasLayer fields
        in.read(reinterpret_cast<char*>(&this->p_mean_molecular_weight), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_adiabatic_index),       sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_reference_temperature), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_reference_density),     sizeof(double));

        if (!in) {
            throw std::runtime_error("TidalPy: failed to read GasLayer binary data");
        }

        // Attached rheology models (presence flag + recursive record each).
        this->read_rheology_binary(in, force);

        this->update_physicals();
    }

protected:
    double p_mean_molecular_weight = 2.0e-3;   // [kg/mol]
    double p_adiabatic_index       = 1.4;       // γ = c_p/c_v [dimensionless]
    double p_reference_temperature = 300.0;     // [K]
    double p_reference_density     = 1.0;       // [kg/m³]
};

} // namespace tidalpy
