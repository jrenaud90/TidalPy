#pragma once
/*
 * gas_.hpp: c_GasLayer, an ideal-gas layer built on c_PhysicsLayer.
 *
 * Adds the ideal-gas parameters an attached EOS model reads (mean molecular weight, adiabatic index, and the
 * reference state). No phase changes, no solidus or liquidus, and no cooling or radiogenics sub-models. All MKS.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::GasLayer (103)
 *   payload:
 *     [all c_BaseLayer fields: same byte layout as the BaseLayer binary payload]
 *     [all c_PhysicsLayer additions: love_numbers k/h/l re+im (6×8), the three classification flags,
 *      temperature, use_thermal_eos, use_heating]
 *     mean_molecular_weight  (double, 8)
 *     adiabatic_index               (double, 8)
 *     reference_temperature       (double, 8)
 *     reference_density       (double, 8)
 *     eos_model       presence flag (uint8_t, 1) + (if present) its binary record
 *     shear_rheology  presence flag (uint8_t, 1) + (if present) its binary record
 *     bulk_rheology   presence flag (uint8_t, 1) + (if present) its binary record
 *   The attached material EOS model and the two rheology models are serialized recursively: the three
 *   presence flags belong to this payload and each nested model follows as its own record. The EOS profile data
 *   is not serialized; re-run the world EOS solve after loading.
 */

#include <cmath>
#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

#include "physics_.hpp"

namespace tidalpy {

// Construction parameters for c_GasLayer: c_PhysicsConfig plus the ideal-gas thermodynamic fields.
struct c_GasConfig : public c_PhysicsConfig {
    double mean_molecular_weight = 2.0e-3;    // [kg/mol] hydrogen default
    double adiabatic_index       = 1.4;       // γ = c_p/c_v [dimensionless]
    double reference_temperature = 300.0;     // [K]
    double reference_density     = 1.0;       // [kg/m³]

    // A gas carries no shear stress, so the radial solver treats a gas layer as a (static) liquid layer.
    c_GasConfig() { this->is_solid = false; }
};

class c_GasLayer : public c_PhysicsLayer {
public:
    // Construction
    c_GasLayer() = default;

    explicit c_GasLayer(const c_GasConfig& cfg)
        : c_PhysicsLayer(cfg),
          p_mean_molecular_weight(cfg.mean_molecular_weight),
          p_adiabatic_index(cfg.adiabatic_index),
          p_reference_temperature(cfg.reference_temperature),
          p_reference_density(cfg.reference_density)
    {}

    ~c_GasLayer() override = default;

    // The unique_ptr members inherited from c_PhysicsLayer delete the implicit copy assignment, so an explicit
    // one is needed for Cython stack allocation. Cython temporaries always have null model pointers.
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

    // Property getters (const, MKS)
    uint32_t get_layer_class_id() const noexcept override {
        return static_cast<uint32_t>(BinaryClassID::GasLayer);
    }
    double get_mean_molecular_weight() const noexcept { return this->p_mean_molecular_weight; }
    double get_adiabatic_index()       const noexcept { return this->p_adiabatic_index; }
    double get_reference_temperature() const noexcept { return this->p_reference_temperature; }
    double get_reference_density()     const noexcept { return this->p_reference_density; }

    // Binary I/O
    void write_binary(std::ostream& out) const override {
        const auto     name_len = static_cast<uint32_t>(this->p_name.size());
        const auto     mat_len  = static_cast<uint32_t>(this->p_material_name.size());
        const uint64_t payload  =
            sizeof(double)   * 2 +           // p_radius, p_mass
            sizeof(uint32_t) + name_len +    // name length + bytes
            sizeof(int32_t)  +               // layer_index
            sizeof(double)   +               // radius_inner
            sizeof(uint32_t) + mat_len +     // material_name length + bytes
            sizeof(uint8_t)  * 2 +           // is_tidal, is_volume_fixed
            sizeof(double)   +               // tidal_scale
            sizeof(uint8_t)  +               // tidal_scale_method
            sizeof(double)   * 6 +           // love_numbers k/h/l re+im
            sizeof(uint8_t)  * 3 +           // is_solid, is_static, is_incompressible
            material_law_bytes() +           // temperature, use_thermal_eos, use_heating
            sizeof(double)   * 4 +           // GasLayer fields
            optional_binary_flag_bytes() +         // material EOS model presence flag
            this->physics_models_presence_bytes(); // shear and bulk rheology presence flags

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
        const uint8_t is_volume_fixed_byte = static_cast<uint8_t>(this->p_is_volume_fixed);
        out.write(reinterpret_cast<const char*>(&is_tidal_byte),        sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&is_volume_fixed_byte), sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&this->p_tidal_scale),  sizeof(double));
        const uint8_t scale_method_byte = static_cast<uint8_t>(this->p_tidal_scale_method);
        out.write(reinterpret_cast<const char*>(&scale_method_byte), sizeof(uint8_t));

        // c_PhysicsLayer fields
        auto write_complex = [&](const std::complex<double>& c) {
            const double re = c.real(), im = c.imag();
            out.write(reinterpret_cast<const char*>(&re), sizeof(double));
            out.write(reinterpret_cast<const char*>(&im), sizeof(double));
        };
        write_complex(this->p_love_numbers.k);
        write_complex(this->p_love_numbers.h);
        write_complex(this->p_love_numbers.l);

        // Radial-solver layer classification flags (mirrors c_PhysicsLayer's layout).
        const uint8_t is_solid_byte          = static_cast<uint8_t>(this->p_is_solid);
        const uint8_t is_static_byte         = static_cast<uint8_t>(this->p_is_static);
        const uint8_t is_incompressible_byte = static_cast<uint8_t>(this->p_is_incompressible);
        out.write(reinterpret_cast<const char*>(&is_solid_byte),          sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&is_static_byte),         sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&is_incompressible_byte), sizeof(uint8_t));
        this->write_material_law_binary(out);

        // c_GasLayer fields
        out.write(reinterpret_cast<const char*>(&this->p_mean_molecular_weight), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_adiabatic_index),       sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_reference_temperature), sizeof(double));
        out.write(reinterpret_cast<const char*>(&this->p_reference_density),     sizeof(double));

        if (!out) {
            throw std::runtime_error("TidalPy: failed to write GasLayer binary data");
        }

        this->write_eos_model_binary(out);
        this->write_physics_models_binary(out);
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
        uint8_t is_volume_fixed_byte = 0;
        in.read(reinterpret_cast<char*>(&is_volume_fixed_byte), sizeof(uint8_t));
        this->p_is_volume_fixed = static_cast<bool>(is_volume_fixed_byte);

        in.read(reinterpret_cast<char*>(&this->p_tidal_scale), sizeof(double));

        uint8_t scale_method_byte = 0;
        in.read(reinterpret_cast<char*>(&scale_method_byte), sizeof(uint8_t));
        this->p_tidal_scale_method = static_cast<c_TidalScaleMethod>(scale_method_byte);

        // c_PhysicsLayer fields
        auto read_complex = [&](std::complex<double>& c) {
            double re = 0.0, im = 0.0;
            in.read(reinterpret_cast<char*>(&re), sizeof(double));
            in.read(reinterpret_cast<char*>(&im), sizeof(double));
            c = std::complex<double>(re, im);
        };
        read_complex(this->p_love_numbers.k);
        read_complex(this->p_love_numbers.h);
        read_complex(this->p_love_numbers.l);

        // Radial-solver layer classification flags (mirrors c_PhysicsLayer's layout).
        uint8_t is_solid_byte = 0;
        uint8_t is_static_byte = 0;
        uint8_t is_incompressible_byte = 0;
        in.read(reinterpret_cast<char*>(&is_solid_byte),          sizeof(uint8_t));
        in.read(reinterpret_cast<char*>(&is_static_byte),         sizeof(uint8_t));
        in.read(reinterpret_cast<char*>(&is_incompressible_byte), sizeof(uint8_t));
        this->p_is_solid          = static_cast<bool>(is_solid_byte);
        this->p_is_static         = static_cast<bool>(is_static_byte);
        this->p_is_incompressible = static_cast<bool>(is_incompressible_byte);
        this->read_material_law_binary(in);

        // c_GasLayer fields
        in.read(reinterpret_cast<char*>(&this->p_mean_molecular_weight), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_adiabatic_index),       sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_reference_temperature), sizeof(double));
        in.read(reinterpret_cast<char*>(&this->p_reference_density),     sizeof(double));

        if (!in) {
            throw std::runtime_error("TidalPy: failed to read GasLayer binary data");
        }

        this->read_eos_model_binary(in, force);
        this->read_physics_models_binary(in, force);

        this->update_physicals();
    }

protected:
    double p_mean_molecular_weight = 2.0e-3;   // [kg/mol]
    double p_adiabatic_index       = 1.4;       // γ = c_p/c_v [dimensionless]
    double p_reference_temperature = 300.0;     // [K]
    double p_reference_density     = 1.0;       // [kg/m³]
};

} // namespace tidalpy
