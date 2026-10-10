#pragma once
/* Shear-modulus laws: a phase's static (unrelaxed) shear modulus [Pa] at a point (pressure, temperature, radius).
 * Frequency dependence is the rheology's job, and melt weakening the material's. All MKS.
 *
 * A law may give a negative or vanishing value far from its fit (a steep temperature derivative, say); the material
 * floors the result at [numerical] minimum_modulus.
 */

#include <cmath>
#include <istream>
#include <memory>
#include <string>
#include <vector>

#include "broadcast_.hpp"
#include "physics_base_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"
#include "thermo_point_.hpp"
#include "pressure_laws_.hpp"   // d_EOS_REFERENCE_TEMPERATURE
#include "../../constants_.hpp"
#include "../../Utilities/arrays/table_lookup_.hpp"

namespace tidalpy {

class c_ShearModulusBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "shear modulus";

    explicit c_ShearModulusBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_ShearModulusBase() override = default;

    // Static shear modulus [Pa] at a point.
    virtual double calc_shear_modulus(const c_ThermoPoint& point) const noexcept = 0;

    // Element-wise over pressure, temperature, and radius; each holds one value per point or a single value used at
    // every point.
    void calc_shear_modulus_vectorize(
            const std::vector<double>& pressure,
            const std::vector<double>& temperature,
            const std::vector<double>& radius,
            std::vector<double>& out_shear_modulus) const {
        const std::size_t num_points = c_broadcast_length(
            {pressure.size(), temperature.size(), radius.size()}, "calc_shear_modulus_vectorize");
        const std::size_t pressure_stride    = c_broadcast_stride(pressure.size());
        const std::size_t temperature_stride = c_broadcast_stride(temperature.size());
        const std::size_t radius_stride      = c_broadcast_stride(radius.size());
        out_shear_modulus.resize(num_points);
        c_ThermoPoint point;
        for (std::size_t i = 0; i < num_points; ++i) {
            point.pressure       = pressure[i * pressure_stride];
            point.temperature    = temperature[i * temperature_stride];
            point.radius         = radius[i * radius_stride];
            out_shear_modulus[i] = this->calc_shear_modulus(point);
        }
    }
};

// A constant shear modulus (alias "const").
class c_ConstantShearModulus final : public c_SpecModel<c_ConstantShearModulus, c_ShearModulusBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ConstantShearModulus;

    static const std::vector<c_ParamSpec<c_ConstantShearModulus>>& parameter_specs() {
        using Self = c_ConstantShearModulus;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"shear_modulus", "shear_modulus_pa", &Self::p_shear_modulus, 5.0e10, c_ParamBounds::NonNegative,
             "Shear modulus [Pa]."},
        };
        return specs;
    }

    c_ConstantShearModulus() : c_ConstantShearModulus(c_ParamMap{}) {}
    explicit c_ConstantShearModulus(const c_ParamMap& params) : c_SpecModel("constant") {
        this->p_initialize(params);
    }

    double calc_shear_modulus(const c_ThermoPoint& /*point*/) const noexcept override {
        return this->p_shear_modulus;
    }

protected:
    double p_shear_modulus = 0.0;
};

// Linear in pressure and temperature: mu = mu0 + mu'_P (P - P_ref) + mu'_T (T - T_ref). Without a finite temperature
// the temperature term is left out.
class c_LinearShearModulus final : public c_SpecModel<c_LinearShearModulus, c_ShearModulusBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::LinearShearModulus;

    static const std::vector<c_ParamSpec<c_LinearShearModulus>>& parameter_specs() {
        using Self = c_LinearShearModulus;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"shear_modulus", "shear_modulus_pa", &Self::p_shear_modulus, 5.0e10, c_ParamBounds::NonNegative,
             "Shear modulus at the reference pressure and temperature, mu0 [Pa]."},
            {"pressure_derivative", "pressure_derivative", &Self::p_pressure_derivative, 0.0, c_ParamBounds::Finite,
             "d mu / d P [dimensionless]."},
            {"temperature_derivative", "temperature_derivative_pa_k", &Self::p_temperature_derivative, 0.0,
             c_ParamBounds::Finite, "d mu / d T [Pa K-1]."},
            {"reference_pressure", "reference_pressure_pa", &Self::p_reference_pressure, 0.0, c_ParamBounds::Finite,
             "Pressure where mu0 applies [Pa]."},
            {"reference_temperature", "reference_temperature_k", &Self::p_reference_temperature,
             d_EOS_REFERENCE_TEMPERATURE, c_ParamBounds::Positive, "Temperature where mu0 applies [K]."},
        };
        return specs;
    }

    c_LinearShearModulus() : c_LinearShearModulus(c_ParamMap{}) {}
    explicit c_LinearShearModulus(const c_ParamMap& params) : c_SpecModel("linear") { this->p_initialize(params); }

    double calc_shear_modulus(const c_ThermoPoint& point) const noexcept override {
        double shear = this->p_shear_modulus + this->p_pressure_derivative * (point.pressure - this->p_reference_pressure);
        if (std::isfinite(point.temperature)) {
            shear += this->p_temperature_derivative * (point.temperature - this->p_reference_temperature);
        }
        return shear;
    }

protected:
    double p_shear_modulus          = 0.0;
    double p_pressure_derivative    = 0.0;
    double p_temperature_derivative = 0.0;
    double p_reference_pressure     = 0.0;
    double p_reference_temperature  = 0.0;
};

// A shear-modulus profile tabulated in radius (aliases "interp", "interpolated"), linear between the table points and
// held at the end values beyond them.
class c_InterpolatedShearModulus final : public c_SpecModel<c_InterpolatedShearModulus, c_ShearModulusBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::InterpolatedShearModulus;

    static const std::vector<c_ParamSpec<c_InterpolatedShearModulus>>& parameter_specs() {
        using Self = c_InterpolatedShearModulus;
        static const std::vector<c_ParamSpec<Self>> specs = {
            // The default tables hold one point, a constant profile.
            {"radius", "radius_m", &Self::p_radius, 0.0, c_ParamBounds::Finite, "Table radii, ascending [m].",
             {0.0}},
            {"shear_modulus", "shear_modulus_pa", &Self::p_shear_modulus, 0.0, c_ParamBounds::NonNegative,
             "Shear modulus at each table radius [Pa].", {5.0e10}},
        };
        return specs;
    }

    c_InterpolatedShearModulus() : c_InterpolatedShearModulus(c_ParamMap{}) {}
    explicit c_InterpolatedShearModulus(const c_ParamMap& params) : c_SpecModel("interpolate") {
        this->p_initialize(params);
    }

    double calc_shear_modulus(const c_ThermoPoint& point) const noexcept override {
        return this->p_lookup.interpolate(point.radius, this->p_radius, this->p_shear_modulus);
    }

protected:
    void p_validate() const override {
        c_check_table(this->p_describe(), this->p_radius, {&this->p_shear_modulus});
        if (this->p_shear_modulus.size() != this->p_radius.size()) {
            throw std::invalid_argument(this->p_describe() + " needs a 'shear_modulus_pa' value at every table radius.");
        }
    }
    void p_update_derived() noexcept override { this->p_lookup.build(this->p_radius); }

    std::vector<double> p_radius;
    std::vector<double> p_shear_modulus;
    c_TableLookup p_lookup;
};

inline const c_ModelRegistry<c_ShearModulusBase>& c_shear_modulus_registry() {
    static const c_ModelRegistry<c_ShearModulusBase> registry = {
        {{"constant", "const"}, BinaryClassID::ConstantShearModulus,
         &c_make_entry<c_ShearModulusBase, c_ConstantShearModulus>},
        {{"linear"}, BinaryClassID::LinearShearModulus, &c_make_entry<c_ShearModulusBase, c_LinearShearModulus>},
        {{"interpolate", "interp", "interpolated"}, BinaryClassID::InterpolatedShearModulus,
         &c_make_entry<c_ShearModulusBase, c_InterpolatedShearModulus>},
    };
    return registry;
}

inline std::unique_ptr<c_ShearModulusBase> c_find_shear_modulus(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_shear_modulus_registry(), model_name, params);
}

inline std::unique_ptr<c_ShearModulusBase> c_shear_modulus_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_shear_modulus_registry(), in, force);
}

inline std::string c_shear_modulus_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_shear_modulus_registry(), model_name);
}

inline std::vector<std::string> c_shear_modulus_model_names() {
    return c_model_names(c_shear_modulus_registry());
}

}  // namespace tidalpy
