#pragma once
/*
 * viscosity_.hpp — TidalPy (solid/liquid) viscosity models.
 *
 * Inherits c_ViscosityBase (viscosity_base_.hpp), which itself inherits
 * c_PhysicsBase. Each model implements calc_viscosity(temperature, pressure),
 * returning the dynamic viscosity [Pa·s].
 *
 * Models (with config aliases handled by the factory):
 *   c_ConstantViscosity   (alias "const")      — temperature/pressure independent.
 *   c_ReferenceViscosity  (alias "ref")        — relative-activation law.
 *   c_ArrheniusViscosity  (alias "arr")        — Arrhenius flow law.
 *
 * All quantities MKS. The math mirrors the validated legacy implementation in
 * TidalPy/rheology/viscosity/viscosity_models.py. The molar gas constant R comes
 * from the shared TidalPy config (tidalpy_config_ptr->d_R).
 *
 * References
 * ----------
 * - Moore (2006) — Arrhenius flow law (activation energy/volume).
 * - Henning (2009) — reference-viscosity (relative activation) law.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::<Model> (801-803)
 *   payload: model_name length (uint32_t) | model_name bytes | model params (doubles)
 *   Constant writes 1, Reference writes 4, Arrhenius writes 8 scalars.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants_.hpp"      // tidalpy_config_ptr->d_R, TidalPyConstants::d_EPS/d_INF
#include "viscosity_base_.hpp"
#include "../Utilities_x/math_x/numerics_.hpp"  // c_safe_pow

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_ViscosityConfig — combined construction parameters for all viscosity models.
// Each model reads only the fields it needs.
// -------------------------------------------------------------------------------
struct c_ViscosityConfig {
    // Constant / Reference.
    double reference_viscosity   = 1.0e22;   // [Pa·s]
    double reference_temperature = 1000.0;   // [K]  (reference model)

    // Reference / Arrhenius.
    double molar_activation_energy = 3.0e5;  // E_a [J/mol]
    double molar_activation_volume = 0.0;    // V_a [m^3/mol]

    // Arrhenius-only.
    double arrhenius_coeff = 1.0;            // pre-exponential A
    double stress          = 1.0;            // applied stress σ [Pa]
    double stress_expo     = 1.0;            // stress exponent n (power 1−n)
    double grain_size      = 1.0e-3;         // grain size d [m]
    double grain_size_expo = 0.0;            // grain-size exponent m
    bool   additional_temp_dependence = false;  // multiply by T if true
};

// -------------------------------------------------------------------------------
// Lower-case a model name for case-insensitive factory lookup.
// -------------------------------------------------------------------------------
inline std::string visc_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// =====================================================================================================================
// Viscosity models
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_ConstantViscosity — viscosity independent of temperature and pressure
// (alias "const"/"constant").
// -------------------------------------------------------------------------------
class c_ConstantViscosity : public c_ViscosityBase {
public:
    c_ConstantViscosity() : c_ViscosityBase("constant") {}
    explicit c_ConstantViscosity(const c_ViscosityConfig& cfg)
        : c_ViscosityBase("constant"),
          p_reference_viscosity(cfg.reference_viscosity) {}
    ~c_ConstantViscosity() override = default;

    double get_reference_viscosity() const noexcept { return this->p_reference_viscosity; }

    double calc_viscosity(double /*temperature_k*/, double /*pressure_pa*/) const override {
        return this->p_reference_viscosity;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::ConstantViscosity),
                                   {this->p_reference_viscosity});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 1);
        this->p_reference_viscosity = params[0];
    }

protected:
    double p_reference_viscosity = 1.0e22;
};

// -------------------------------------------------------------------------------
// c_ReferenceViscosity — relative-activation law (alias "ref"/"reference"):
//
//   η = η_ref · exp( ((E_a + P·V_a) / R) · (1/T − 1/T_ref) )
// -------------------------------------------------------------------------------
class c_ReferenceViscosity : public c_ViscosityBase {
public:
    c_ReferenceViscosity() : c_ViscosityBase("reference") {}
    explicit c_ReferenceViscosity(const c_ViscosityConfig& cfg)
        : c_ViscosityBase("reference"),
          p_reference_viscosity(cfg.reference_viscosity),
          p_reference_temperature(cfg.reference_temperature),
          p_molar_activation_energy(cfg.molar_activation_energy),
          p_molar_activation_volume(cfg.molar_activation_volume) {}
    ~c_ReferenceViscosity() override = default;

    double get_reference_viscosity()      const noexcept { return this->p_reference_viscosity; }
    double get_reference_temperature()    const noexcept { return this->p_reference_temperature; }
    double get_molar_activation_energy()  const noexcept { return this->p_molar_activation_energy; }
    double get_molar_activation_volume()  const noexcept { return this->p_molar_activation_volume; }

    double calc_viscosity(double temperature_k, double pressure_pa) const override {
        const double R = tidalpy_config_ptr->d_R;
        // A non-positive temperature is the cold limit: effectively rigid (infinite viscosity),
        // which the rheology models treat as a purely elastic response.
        if (temperature_k <= TidalPyConstants::d_EPS
            || this->p_reference_temperature <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_INF;
        }
        const double delta_inv_temp = (1.0 / temperature_k) - (1.0 / this->p_reference_temperature);
        const double exponent =
            ((this->p_molar_activation_energy + pressure_pa * this->p_molar_activation_volume) / R)
            * delta_inv_temp;
        // Plain exp: an overflowing (very cold) exponent saturates to the same rigid limit.
        return this->p_reference_viscosity * std::exp(exponent);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::ReferenceViscosity),
            {this->p_reference_viscosity, this->p_reference_temperature,
             this->p_molar_activation_energy, this->p_molar_activation_volume});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 4);
        this->p_reference_viscosity     = params[0];
        this->p_reference_temperature   = params[1];
        this->p_molar_activation_energy = params[2];
        this->p_molar_activation_volume = params[3];
    }

protected:
    double p_reference_viscosity     = 1.0e22;
    double p_reference_temperature   = 1000.0;
    double p_molar_activation_energy = 3.0e5;
    double p_molar_activation_volume = 0.0;
};

// -------------------------------------------------------------------------------
// c_ArrheniusViscosity — Arrhenius flow law (alias "arr"/"arrhenius"):
//
//   η = A · σ^(1−n) · d^m · exp( (E_a + P·V_a) / (R·T) )    [· T if additional_temp_dependence]
// -------------------------------------------------------------------------------
class c_ArrheniusViscosity : public c_ViscosityBase {
public:
    c_ArrheniusViscosity() : c_ViscosityBase("arrhenius") {}
    explicit c_ArrheniusViscosity(const c_ViscosityConfig& cfg)
        : c_ViscosityBase("arrhenius"),
          p_arrhenius_coeff(cfg.arrhenius_coeff),
          p_stress(cfg.stress),
          p_stress_expo(cfg.stress_expo),
          p_grain_size(cfg.grain_size),
          p_grain_size_expo(cfg.grain_size_expo),
          p_molar_activation_energy(cfg.molar_activation_energy),
          p_molar_activation_volume(cfg.molar_activation_volume),
          p_additional_temp_dependence(cfg.additional_temp_dependence) {}
    ~c_ArrheniusViscosity() override = default;

    double get_arrhenius_coeff()            const noexcept { return this->p_arrhenius_coeff; }
    double get_stress()                     const noexcept { return this->p_stress; }
    double get_stress_expo()                const noexcept { return this->p_stress_expo; }
    double get_grain_size()                 const noexcept { return this->p_grain_size; }
    double get_grain_size_expo()            const noexcept { return this->p_grain_size_expo; }
    double get_molar_activation_energy()    const noexcept { return this->p_molar_activation_energy; }
    double get_molar_activation_volume()    const noexcept { return this->p_molar_activation_volume; }
    bool   get_additional_temp_dependence() const noexcept { return this->p_additional_temp_dependence; }

    double calc_viscosity(double temperature_k, double pressure_pa) const override {
        const double R = tidalpy_config_ptr->d_R;
        // A non-positive temperature is the cold limit: effectively rigid (infinite viscosity),
        // which the rheology models treat as a purely elastic response.
        if (temperature_k <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_INF;
        }
        const double exponent =
            (this->p_molar_activation_energy + pressure_pa * this->p_molar_activation_volume)
            / (R * temperature_k);
        // Plain exp: an overflowing (very cold) exponent saturates to the same rigid limit.
        double viscosity = this->p_arrhenius_coeff
                         * c_safe_pow(this->p_stress, 1.0 - this->p_stress_expo)
                         * c_safe_pow(this->p_grain_size, this->p_grain_size_expo)
                         * std::exp(exponent);
        if (this->p_additional_temp_dependence) {
            viscosity *= temperature_k;
        }
        return viscosity;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::ArrheniusViscosity),
            {this->p_arrhenius_coeff, this->p_stress, this->p_stress_expo,
             this->p_grain_size, this->p_grain_size_expo,
             this->p_molar_activation_energy, this->p_molar_activation_volume,
             this->p_additional_temp_dependence ? 1.0 : 0.0});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 8);
        this->p_arrhenius_coeff            = params[0];
        this->p_stress                     = params[1];
        this->p_stress_expo                = params[2];
        this->p_grain_size                 = params[3];
        this->p_grain_size_expo            = params[4];
        this->p_molar_activation_energy    = params[5];
        this->p_molar_activation_volume    = params[6];
        this->p_additional_temp_dependence = (params[7] != 0.0);
    }

protected:
    double p_arrhenius_coeff            = 1.0;
    double p_stress                     = 1.0;
    double p_stress_expo                = 1.0;
    double p_grain_size                 = 1.0e-3;
    double p_grain_size_expo            = 0.0;
    double p_molar_activation_energy    = 3.0e5;
    double p_molar_activation_volume    = 0.0;
    bool   p_additional_temp_dependence = false;
};

// =====================================================================================================================
// Factory
// =====================================================================================================================

enum class c_ViscosityModel : uint8_t {
    Arrhenius = 0,
    Reference = 1,
    Constant  = 2,
};

// Map a (case-insensitive) model name or alias to a c_ViscosityModel enum value.
// Throws std::invalid_argument on an unknown name.
inline c_ViscosityModel c_viscosity_model_from_name(const std::string& model_name) {
    const std::string name = visc_to_lower(model_name);
    if (name == "arrhenius" || name == "arr")   { return c_ViscosityModel::Arrhenius; }
    if (name == "reference" || name == "ref")   { return c_ViscosityModel::Reference; }
    if (name == "constant"  || name == "const") { return c_ViscosityModel::Constant; }
    throw std::invalid_argument("TidalPy: unknown viscosity model name '" + model_name + "'");
}

// Build the viscosity model named by the enum; returns an owning unique_ptr.
inline std::unique_ptr<c_ViscosityBase> c_find_viscosity(
        c_ViscosityModel model, const c_ViscosityConfig& cfg) {
    switch (model) {
        case c_ViscosityModel::Arrhenius: return std::make_unique<c_ArrheniusViscosity>(cfg);
        case c_ViscosityModel::Reference: return std::make_unique<c_ReferenceViscosity>(cfg);
        case c_ViscosityModel::Constant:  return std::make_unique<c_ConstantViscosity>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_ViscosityModel enum value");
}

// Name overload.
inline std::unique_ptr<c_ViscosityBase> c_find_viscosity(
        const std::string& model_name, const c_ViscosityConfig& cfg) {
    return c_find_viscosity(c_viscosity_model_from_name(model_name), cfg);
}

// Reconstruct a viscosity model from a binary stream (peek class id -> build -> read).
inline std::unique_ptr<c_ViscosityBase> c_viscosity_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_ViscosityBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::ArrheniusViscosity: model = std::make_unique<c_ArrheniusViscosity>(); break;
        case BinaryClassID::ReferenceViscosity: model = std::make_unique<c_ReferenceViscosity>(); break;
        case BinaryClassID::ConstantViscosity:  model = std::make_unique<c_ConstantViscosity>();  break;
        default:
            throw std::runtime_error("TidalPy: unknown viscosity class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
