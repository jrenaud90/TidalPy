#pragma once
/* TidalPy (solid and liquid) viscosity models. All quantities MKS.
 *
 * Each model declares its parameters in one table (c_SpecModel, spec_model_.hpp), which gives its constructor,
 * validation, config entries, and binary record.
 *
 * References
 * ----------
 * - Moore (2006): Arrhenius flow law (activation energy and volume).
 * - Henning (2009): reference-viscosity (relative activation) law.
 */

#include <cmath>
#include <istream>
#include <memory>
#include <string>
#include <vector>

#include "constants_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"
#include "viscosity_base_.hpp"
#include "../Utilities/math/numerics_.hpp"  // c_safe_pow

namespace tidalpy {

// The gas constant R [J/mol/K] from the shared runtime config; NaN when the pointer was never wired, so a missing
// initialization shows up in the viscosity.
inline double c_viscosity_gas_constant() noexcept {
    return (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_R : TidalPyConstants::d_NAN;
}

// Viscosity independent of temperature and pressure (alias "const").
class c_ConstantViscosity final : public c_SpecModel<c_ConstantViscosity, c_ViscosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ConstantViscosity;

    static const std::vector<c_ParamSpec<c_ConstantViscosity>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_ConstantViscosity>> specs = {
            {"reference_viscosity", "reference_viscosity_pas", &c_ConstantViscosity::p_reference_viscosity, 1.0e22,
             c_ParamBounds::Positive, "Viscosity [Pa s]."},
        };
        return specs;
    }

    c_ConstantViscosity() : c_ConstantViscosity(c_ParamMap{}) {}
    explicit c_ConstantViscosity(const c_ParamMap& params) : c_SpecModel("constant") { this->p_initialize(params); }

    double calc_viscosity(double /*temperature*/, double /*pressure*/) const noexcept override {
        return this->p_reference_viscosity;
    }

protected:
    double p_reference_viscosity = 0.0;
};

// Relative-activation law (alias "ref"):
//   eta = eta_ref * exp( (E_a / R) (1/T - 1/T_ref) + P V_a / (R T) )
// The activation energy is anchored at the reference temperature; the reference viscosity is at zero pressure, so a
// positive activation volume always raises the viscosity, by exp(P V_a / (R T)).
class c_ReferenceViscosity final : public c_SpecModel<c_ReferenceViscosity, c_ViscosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ReferenceViscosity;

    static const std::vector<c_ParamSpec<c_ReferenceViscosity>>& parameter_specs() {
        using Self = c_ReferenceViscosity;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"reference_viscosity", "reference_viscosity_pas", &Self::p_reference_viscosity, 1.0e22,
             c_ParamBounds::Positive, "Viscosity at the reference temperature and zero pressure [Pa s]."},
            {"reference_temperature", "reference_temperature_k", &Self::p_reference_temperature, 1000.0,
             c_ParamBounds::Positive, "Reference temperature [K]."},
            {"molar_activation_energy", "molar_activation_energy_j_mol", &Self::p_molar_activation_energy, 3.0e5,
             c_ParamBounds::NonNegative, "Molar activation energy E_a [J mol-1]."},
            {"molar_activation_volume", "molar_activation_volume_m3_mol", &Self::p_molar_activation_volume, 0.0,
             c_ParamBounds::Finite, "Molar activation volume V_a [m3 mol-1]."},
        };
        return specs;
    }

    c_ReferenceViscosity() : c_ReferenceViscosity(c_ParamMap{}) {}
    explicit c_ReferenceViscosity(const c_ParamMap& params) : c_SpecModel("reference") {
        this->p_initialize(params);
    }

    double calc_viscosity(double temperature, double pressure) const noexcept override {
        // Cold limit: rigid, which the rheology models read as a purely elastic response.
        if (temperature <= TidalPyConstants::d_EPS) { return TidalPyConstants::d_INF; }
        const double R = c_viscosity_gas_constant();
        const double exponent =
            (this->p_molar_activation_energy / R) * ((1.0 / temperature) - (1.0 / this->p_reference_temperature))
            + (pressure * this->p_molar_activation_volume) / (R * temperature);
        // Plain exp: an overflowing (very cold) exponent saturates to that same rigid limit.
        return this->p_reference_viscosity * std::exp(exponent);
    }

protected:
    double p_reference_viscosity     = 0.0;
    double p_reference_temperature   = 0.0;
    double p_molar_activation_energy = 0.0;
    double p_molar_activation_volume = 0.0;
};

// Arrhenius flow law (alias "arr"):
//   eta = A * sigma^(1-n) * d^m * exp( (E_a + P * V_a) / (R * T) ), times T when additional_temp_dependence is set.
class c_ArrheniusViscosity final : public c_SpecModel<c_ArrheniusViscosity, c_ViscosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ArrheniusViscosity;

    static const std::vector<c_ParamSpec<c_ArrheniusViscosity>>& parameter_specs() {
        using Self = c_ArrheniusViscosity;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"arrhenius_coeff", "arrhenius_coeff", &Self::p_arrhenius_coeff, 1.0,
             c_ParamBounds::Positive, "Pre-exponential coefficient A; its units follow from the exponents."},
            {"stress", "stress_pa", &Self::p_stress, 1.0,
             c_ParamBounds::Positive, "Applied stress sigma [Pa]."},
            {"stress_expo", "stress_expo", &Self::p_stress_expo, 1.0,
             c_ParamBounds::Finite, "Stress exponent n: 1 for diffusion creep, above 1 for dislocation creep."},
            {"grain_size", "grain_size_m", &Self::p_grain_size, 1.0e-3,
             c_ParamBounds::Positive, "Grain size d [m]."},
            {"grain_size_expo", "grain_size_expo", &Self::p_grain_size_expo, 0.0,
             c_ParamBounds::Finite, "Grain-size exponent m; 0 removes the grain-size dependence."},
            {"molar_activation_energy", "molar_activation_energy_j_mol", &Self::p_molar_activation_energy, 3.0e5,
             c_ParamBounds::NonNegative, "Molar activation energy E_a [J mol-1]."},
            {"molar_activation_volume", "molar_activation_volume_m3_mol", &Self::p_molar_activation_volume, 0.0,
             c_ParamBounds::Finite, "Molar activation volume V_a [m3 mol-1]."},
            {"additional_temp_dependence", "additional_temp_dependence", &Self::p_additional_temp_dependence, 0.0,
             c_ParamBounds::Any, "Multiply the law by T, as the Goldsby and Kohlstedt form does."},
        };
        return specs;
    }

    c_ArrheniusViscosity() : c_ArrheniusViscosity(c_ParamMap{}) {}
    explicit c_ArrheniusViscosity(const c_ParamMap& params) : c_SpecModel("arrhenius") {
        this->p_initialize(params);
    }

    double calc_viscosity(double temperature, double pressure) const noexcept override {
        // Cold limit: rigid, which the rheology models read as a purely elastic response.
        if (temperature <= TidalPyConstants::d_EPS) { return TidalPyConstants::d_INF; }
        const double R = c_viscosity_gas_constant();
        const double exponent =
            (this->p_molar_activation_energy + pressure * this->p_molar_activation_volume) / (R * temperature);
        double viscosity = this->p_prefactor * std::exp(exponent);
        if (this->p_additional_temp_dependence) { viscosity *= temperature; }
        return viscosity;
    }

protected:
    // A sigma^(1 - n) d^m, the temperature- and pressure-independent part of the law. It associates as
    // (A sigma^(1 - n)) d^m, the left-to-right order of the whole product, so forming it once gives the viscosity to
    // the last bit.
    void p_update_derived() noexcept override {
        this->p_prefactor = this->p_arrhenius_coeff
                          * c_safe_pow(this->p_stress, 1.0 - this->p_stress_expo)
                          * c_safe_pow(this->p_grain_size, this->p_grain_size_expo);
    }

    double p_arrhenius_coeff            = 0.0;
    double p_stress                     = 0.0;
    double p_stress_expo                = 0.0;
    double p_grain_size                 = 0.0;
    double p_grain_size_expo            = 0.0;
    double p_molar_activation_energy    = 0.0;
    double p_molar_activation_volume    = 0.0;
    bool   p_additional_temp_dependence = false;
    // Rebuilt from the parameters after every change; never serialized.
    double p_prefactor = 0.0;
};

inline const c_ModelRegistry<c_ViscosityBase>& c_viscosity_registry() {
    static const c_ModelRegistry<c_ViscosityBase> registry = {
        {{"arrhenius", "arr"},   BinaryClassID::ArrheniusViscosity, &c_make_entry<c_ViscosityBase, c_ArrheniusViscosity>},
        {{"reference", "ref"},   BinaryClassID::ReferenceViscosity, &c_make_entry<c_ViscosityBase, c_ReferenceViscosity>},
        {{"constant", "const"},  BinaryClassID::ConstantViscosity,  &c_make_entry<c_ViscosityBase, c_ConstantViscosity>},
    };
    return registry;
}

// The family's entry points, each one line over the generic registry functions.
inline std::unique_ptr<c_ViscosityBase> c_find_viscosity(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_viscosity_registry(), model_name, params);
}

inline std::unique_ptr<c_ViscosityBase> c_viscosity_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_viscosity_registry(), in, force);
}

inline std::string c_viscosity_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_viscosity_registry(), model_name);
}

inline std::vector<std::string> c_viscosity_model_names() {
    return c_model_names(c_viscosity_registry());
}

}  // namespace tidalpy
