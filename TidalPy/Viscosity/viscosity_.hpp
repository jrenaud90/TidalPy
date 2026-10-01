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

#include "binary_.hpp"
#include "constants_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"
#include "viscosity_base_.hpp"
#include "../Utilities/arrays/table_lookup_.hpp"
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
             c_ParamBounds::PositiveOrInfinite, "Viscosity [Pa s]; infinite for a rigid (purely elastic) material."},
        };
        return specs;
    }

    c_ConstantViscosity() : c_ConstantViscosity(c_ParamMap{}) {}
    explicit c_ConstantViscosity(const c_ParamMap& params) : c_SpecModel("constant") { this->p_initialize(params); }

    double calc_viscosity(const c_ThermoPoint& /*point*/) const noexcept override {
        return this->p_reference_viscosity;
    }

protected:
    double p_reference_viscosity = 0.0;
};

// Relative-activation law (alias "ref"), anchored at a reference temperature and pressure:
//   eta = eta_ref * exp( (E_a + P V_a) / (R T) - (E_a + P_ref V_a) / (R T_ref) )
// so eta = eta_ref at (T_ref, P_ref). With the default P_ref = 0 a positive activation volume always raises the
// viscosity; a deep layer's law is better anchored at a pressure inside the layer.
class c_ReferenceViscosity final : public c_SpecModel<c_ReferenceViscosity, c_ViscosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ReferenceViscosity;

    static const std::vector<c_ParamSpec<c_ReferenceViscosity>>& parameter_specs() {
        using Self = c_ReferenceViscosity;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"reference_viscosity", "reference_viscosity_pas", &Self::p_reference_viscosity, 1.0e22,
             c_ParamBounds::Positive, "Viscosity at the reference temperature and pressure [Pa s]."},
            {"reference_temperature", "reference_temperature_k", &Self::p_reference_temperature, 1000.0,
             c_ParamBounds::Positive, "Reference temperature [K]."},
            {"molar_activation_energy", "molar_activation_energy_j_mol", &Self::p_molar_activation_energy, 3.0e5,
             c_ParamBounds::NonNegative, "Molar activation energy E_a [J mol-1]."},
            {"molar_activation_volume", "molar_activation_volume_m3_mol", &Self::p_molar_activation_volume, 0.0,
             c_ParamBounds::Finite, "Molar activation volume V_a [m3 mol-1]."},
            {"reference_pressure", "reference_pressure_pa", &Self::p_reference_pressure, 0.0,
             c_ParamBounds::Finite, "Reference pressure [Pa]."},
        };
        return specs;
    }

    c_ReferenceViscosity() : c_ReferenceViscosity(c_ParamMap{}) {}
    explicit c_ReferenceViscosity(const c_ParamMap& params) : c_SpecModel("reference") {
        this->p_initialize(params);
    }

    double calc_viscosity(const c_ThermoPoint& point) const noexcept override {
        const double temperature = point.temperature;
        // Cold limit: rigid, which the rheology models read as a purely elastic response.
        if (temperature <= TidalPyConstants::d_EPS) { return TidalPyConstants::d_INF; }
        const double R = c_viscosity_gas_constant();
        const double exponent =
            (this->p_molar_activation_energy + point.pressure * this->p_molar_activation_volume) / (R * temperature)
            - (this->p_molar_activation_energy + this->p_reference_pressure * this->p_molar_activation_volume)
              / (R * this->p_reference_temperature);
        // Plain exp: an overflowing (very cold) exponent saturates to that same rigid limit.
        return this->p_reference_viscosity * std::exp(exponent);
    }

protected:
    double p_reference_viscosity     = 0.0;
    double p_reference_temperature   = 0.0;
    double p_molar_activation_energy = 0.0;
    double p_molar_activation_volume = 0.0;
    double p_reference_pressure      = 0.0;
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

    double calc_viscosity(const c_ThermoPoint& point) const noexcept override {
        const double temperature = point.temperature;
        // Cold limit: rigid, which the rheology models read as a purely elastic response.
        if (temperature <= TidalPyConstants::d_EPS) { return TidalPyConstants::d_INF; }
        const double R = c_viscosity_gas_constant();
        const double exponent =
            (this->p_molar_activation_energy + point.pressure * this->p_molar_activation_volume) / (R * temperature);
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

// A viscosity profile tabulated in radius (alias "interp"), linear between the table points and held at the end
// values beyond them. A seismic profile's quality factor travels in this table for the seismic_q rheology, which reads
// its viscosity input as Q.
class c_InterpolatedViscosity final : public c_SpecModel<c_InterpolatedViscosity, c_ViscosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::InterpolatedViscosity;

    static const std::vector<c_ParamSpec<c_InterpolatedViscosity>>& parameter_specs() {
        using Self = c_InterpolatedViscosity;
        static const std::vector<c_ParamSpec<Self>> specs = {
            // The default tables hold one point, a constant profile.
            {"radius", "radius_m", &Self::p_radius, 0.0, c_ParamBounds::Finite,
             "Table radii, ascending [m].", {0.0}},
            {"viscosity", "viscosity_pas", &Self::p_viscosity, 0.0, c_ParamBounds::Positive,
             "Viscosity at each table radius [Pa s].", {1.0e22}},
        };
        return specs;
    }

    c_InterpolatedViscosity() : c_InterpolatedViscosity(c_ParamMap{}) {}
    explicit c_InterpolatedViscosity(const c_ParamMap& params) : c_SpecModel("interpolate") {
        this->p_initialize(params);
    }

    double calc_viscosity(const c_ThermoPoint& point) const noexcept override {
        return this->p_lookup.interpolate(point.radius, this->p_radius, this->p_viscosity);
    }

protected:
    void p_validate() const override {
        c_check_table(this->p_describe(), this->p_radius, {&this->p_viscosity});
        if (this->p_viscosity.size() != this->p_radius.size()) {
            throw std::invalid_argument(this->p_describe() + " needs a 'viscosity_pas' value at every table radius.");
        }
    }
    void p_update_derived() noexcept override { this->p_lookup.build(this->p_radius); }

    std::vector<double> p_radius;
    std::vector<double> p_viscosity;
    c_TableLookup p_lookup;
};

// Defined below the registry; the composite reads its mechanisms' records through it.
inline std::unique_ptr<c_ViscosityBase> c_viscosity_from_binary(std::istream& in, bool force);

// Several deformation mechanisms acting in parallel (alias "parallel"): their strain rates add at a common stress,
// so 1 / eta = sum of 1 / eta_i, and the weakest mechanism dominates. Published ice and olivine flow laws combine
// diffusion creep, dislocation creep, and grain-boundary sliding this way; a mechanism whose law switches activation
// energy at a temperature (ice near 255 K) is two mechanisms. A rigid (infinite) mechanism adds nothing.
class c_CompositeViscosity final : public c_SpecModel<c_CompositeViscosity, c_ViscosityBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::CompositeViscosity;
    static const std::vector<c_ParamSpec<c_CompositeViscosity>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_CompositeViscosity>> specs;
        return specs;
    }

    using Mechanisms = std::vector<std::shared_ptr<const c_ViscosityBase>>;

    // Without mechanisms (a default instance, which a binary read fills) it holds one constant mechanism.
    c_CompositeViscosity() : c_CompositeViscosity(c_ParamMap{}) {}
    explicit c_CompositeViscosity(const c_ParamMap& params)
        : c_CompositeViscosity(params, Mechanisms{std::make_shared<const c_ConstantViscosity>()}) {}
    c_CompositeViscosity(const c_ParamMap& params, Mechanisms mechanisms)
        : c_SpecModel("composite"), p_mechanisms(std::move(mechanisms)) {
        this->p_initialize(params);
    }

    const Mechanisms& get_mechanisms() const noexcept { return this->p_mechanisms; }

    // The mechanisms as the shared base pointers the Python wrappers hold.
    std::vector<std::shared_ptr<c_PhysicsBase>> get_mechanism_models() const {
        std::vector<std::shared_ptr<c_PhysicsBase>> models;
        for (const auto& mechanism : this->p_mechanisms) { models.push_back(c_share_physics_of(mechanism)); }
        return models;
    }

    double calc_viscosity(const c_ThermoPoint& point) const noexcept override {
        double inverse = 0.0;
        for (const auto& mechanism : this->p_mechanisms) {
            const double viscosity = mechanism->calc_viscosity(point);
            if (std::isnan(viscosity)) { return TidalPyConstants::d_NAN; }
            if (viscosity <= 0.0) { return 0.0; }
            inverse += 1.0 / viscosity;
        }
        return (inverse > 0.0) ? 1.0 / inverse : TidalPyConstants::d_INF;
    }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_SpecModel::append_config_entries(out);
        std::vector<std::vector<c_ConfigEntry>> tables;
        for (const auto& mechanism : this->p_mechanisms) { tables.push_back(mechanism->get_config_entries()); }
        out.push_back(c_config_table_list("mechanisms", tables));
    }

protected:
    void p_validate() const override {
        if (this->p_mechanisms.empty()) {
            throw std::invalid_argument(this->p_describe() + " needs at least one mechanism.");
        }
        for (const auto& mechanism : this->p_mechanisms) {
            if (!mechanism) { throw std::invalid_argument(this->p_describe() + " was given an empty mechanism."); }
        }
    }

    // The parameters (none), then the mechanism count and each mechanism's own record.
    void p_write_payload(std::ostream& out) const override {
        c_SpecModel::p_write_payload(out);
        const uint32_t num_mechanisms = static_cast<uint32_t>(this->p_mechanisms.size());
        out.write(reinterpret_cast<const char*>(&num_mechanisms), sizeof(uint32_t));
        for (const auto& mechanism : this->p_mechanisms) { mechanism->write_binary(out); }
    }

    void p_read_payload(std::istream& in, bool force) override {
        c_SpecModel::p_read_payload(in, force);
        uint32_t num_mechanisms = 0;
        in.read(reinterpret_cast<char*>(&num_mechanisms), sizeof(uint32_t));
        if (!in) { throw std::runtime_error("TidalPy: failed to read a composite viscosity's mechanism count"); }
        check_binary_count(in, num_mechanisms, TIDALPY_BINARY_HEADER_BYTES, "composite viscosity mechanisms");
        Mechanisms mechanisms;
        for (uint32_t mechanism_i = 0; mechanism_i < num_mechanisms; ++mechanism_i) {
            mechanisms.push_back(std::shared_ptr<const c_ViscosityBase>(c_viscosity_from_binary(in, force)));
        }
        this->p_mechanisms = std::move(mechanisms);
        try {
            this->p_validate();
        }
        catch (const std::invalid_argument& mechanism_error) {
            throw std::runtime_error(std::string("TidalPy: corrupt binary data: ") + mechanism_error.what());
        }
    }

    Mechanisms p_mechanisms;
};

// A composite from mechanisms held as shared base pointers (the Python wrappers' form); each must be a viscosity model.
inline std::unique_ptr<c_ViscosityBase> c_make_composite_viscosity(
        const std::vector<std::shared_ptr<c_PhysicsBase>>& mechanisms) {
    c_CompositeViscosity::Mechanisms family_mechanisms;
    for (const auto& mechanism : mechanisms) {
        family_mechanisms.push_back(c_share_as<c_ViscosityBase>(mechanism, "a composite viscosity mechanism"));
    }
    return std::make_unique<c_CompositeViscosity>(c_ParamMap{}, std::move(family_mechanisms));
}

inline const c_ModelRegistry<c_ViscosityBase>& c_viscosity_registry() {
    static const c_ModelRegistry<c_ViscosityBase> registry = {
        {{"arrhenius", "arr"},   BinaryClassID::ArrheniusViscosity, &c_make_entry<c_ViscosityBase, c_ArrheniusViscosity>},
        {{"reference", "ref"},   BinaryClassID::ReferenceViscosity, &c_make_entry<c_ViscosityBase, c_ReferenceViscosity>},
        {{"constant", "const"},  BinaryClassID::ConstantViscosity,  &c_make_entry<c_ViscosityBase, c_ConstantViscosity>},
        {{"interpolate", "interp", "interpolated"},
         BinaryClassID::InterpolatedViscosity, &c_make_entry<c_ViscosityBase, c_InterpolatedViscosity>},
        {{"composite", "parallel"}, BinaryClassID::CompositeViscosity, &c_make_entry<c_ViscosityBase, c_CompositeViscosity>},
    };
    return registry;
}

// The family's entry points, each one line over the generic registry functions.
inline std::unique_ptr<c_ViscosityBase> c_find_viscosity(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_viscosity_registry(), model_name, params);
}

inline std::unique_ptr<c_ViscosityBase> c_viscosity_from_binary(std::istream& in, bool force) {
    return c_model_from_binary(c_viscosity_registry(), in, force);
}

inline std::string c_viscosity_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_viscosity_registry(), model_name);
}

inline std::vector<std::string> c_viscosity_model_names() {
    return c_model_names(c_viscosity_registry());
}

}  // namespace tidalpy
