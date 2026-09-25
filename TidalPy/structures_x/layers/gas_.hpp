#pragma once
/*
 * gas_.hpp: c_GasLayer, an ideal-gas layer built on c_PhysicsLayer.
 *
 * Adds ideal-gas parameters (mean molecular weight, adiabatic index, and a reference state), stored and serialized
 * for a future gas description; nothing reads them yet, and the layer's density comes from its material's law. No
 * phase changes, no solidus or liquidus, and no cooling or radiogenics sub-models. All MKS.
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

    void write_binary(std::ostream& out) const override {
        write_binary_header(
            out, static_cast<uint32_t>(BinaryClassID::GasLayer),
            this->p_base_fields_bytes() + this->p_physics_fields_bytes()
                + sizeof(double) * 4                       // the ideal-gas fields
                + optional_binary_flag_bytes()             // material EOS model presence flag
                + this->physics_models_presence_bytes());  // shear and bulk rheology presence flags
        this->p_write_base_fields(out);
        this->p_write_physics_fields(out);
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
        this->p_begin_read_binary(in, force);
        this->p_read_base_fields(in);
        this->p_read_physics_fields(in);
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
