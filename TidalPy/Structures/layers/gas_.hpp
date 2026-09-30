#pragma once
/*
 * gas_.hpp: c_GasLayer, an ideal-gas layer built on c_PhysicsLayer.
 *
 * Adds ideal-gas parameters (mean molecular weight, adiabatic index, and a reference state), stored and serialized
 * for a future gas description; nothing reads them yet, and the layer's density comes from its material's law. No
 * phase changes, no solidus or liquidus, and no cooling or radiogenics sub-models. All MKS.
 *
 * Binary payload: the c_PhysicsLayer payload, then the mean molecular weight, adiabatic index, reference temperature,
 * and reference density.
 */

#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
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
    uint32_t get_binary_class_id() const override {
        return static_cast<uint32_t>(BinaryClassID::GasLayer);
    }
    double get_mean_molecular_weight() const noexcept { return this->p_mean_molecular_weight; }
    double get_adiabatic_index()       const noexcept { return this->p_adiabatic_index; }
    double get_reference_temperature() const noexcept { return this->p_reference_temperature; }
    double get_reference_density()     const noexcept { return this->p_reference_density; }

    // What load_binary reads a file into first (c_TidalPyBaseClass::make_binary_scratch).
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_GasLayer>();
    }

protected:
    void p_write_payload(std::ostream& out) const override {
        c_PhysicsLayer::p_write_payload(out);
        const double gas_fields[4] = {
            this->p_mean_molecular_weight, this->p_adiabatic_index, this->p_reference_temperature,
            this->p_reference_density};
        out.write(reinterpret_cast<const char*>(gas_fields), sizeof(gas_fields));
    }

    void p_read_payload(std::istream& in, bool force) override {
        c_PhysicsLayer::p_read_payload(in, force);
        double gas_fields[4] = {0.0, 0.0, 0.0, 0.0};
        in.read(reinterpret_cast<char*>(gas_fields), sizeof(gas_fields));
        this->p_mean_molecular_weight = gas_fields[0];
        this->p_adiabatic_index       = gas_fields[1];
        this->p_reference_temperature = gas_fields[2];
        this->p_reference_density     = gas_fields[3];
    }

    double p_mean_molecular_weight = 2.0e-3;   // [kg/mol]
    double p_adiabatic_index       = 1.4;       // γ = c_p/c_v [dimensionless]
    double p_reference_temperature = 300.0;     // [K]
    double p_reference_density     = 1.0;       // [kg/m³]
};

} // namespace tidalpy
