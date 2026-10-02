#pragma once
/* Melting curves: a solidus or liquidus temperature [K] as a function of pressure [Pa]. A material with a liquid
 * phase holds two, and its melt fraction runs linearly between them. All MKS.
 *
 * The Simon and Glatzel (1929) law, T(P) = T0 (1 + (P - P_ref) / a)^(1 / c), fits most planetary melting data with
 * a = c = reference values from the fit; a negative a gives a melting temperature that falls with pressure (ice Ih).
 * Below P_ref the curve holds T0, and where 1 + (P - P_ref) / a reaches zero (past the end of a falling curve) the
 * temperature is zero: a material is only meant to be used inside the pressure range its curve was fitted over.
 *
 * References
 * ----------
 * - Simon and Glatzel (1929), Z. Anorg. Allg. Chem. 178, 309.
 * - Fiquet et al. (2010), Science 329, 1516; Andrault et al. (2011), EPSL 304, 251; Monteux, Andrault, and Samuel
 *   (2016), EPSL 448, 140: peridotite and chondritic-mantle melting fitted as two Simon-Glatzel branches.
 */

#include <algorithm>
#include <cmath>
#include <istream>
#include <memory>
#include <string>
#include <vector>

#include "physics_base_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"
#include "../constants_.hpp"
#include "../Utilities/arrays/table_lookup_.hpp"

namespace tidalpy {

class c_MeltingCurveBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "melting curve";

    explicit c_MeltingCurveBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_MeltingCurveBase() override = default;

    // Melting temperature [K] at a pressure [Pa]; NaN for a non-finite pressure.
    virtual double calc_melting_temperature(double pressure) const noexcept = 0;

    // Slope of the curve, dT_m/dP [K Pa-1], at a pressure [Pa]; zero where the curve is held flat, NaN for a
    // non-finite pressure. The adiabat of a melting range takes the latent heat's pressure term from it.
    virtual double calc_melting_slope(double pressure) const noexcept = 0;

    void calc_melting_temperature_vectorize(
            const std::vector<double>& pressure,
            std::vector<double>& out_temperature) const {
        out_temperature.resize(pressure.size());
        for (std::size_t i = 0; i < pressure.size(); ++i) {
            out_temperature[i] = this->calc_melting_temperature(pressure[i]);
        }
    }

    void calc_melting_slope_vectorize(
            const std::vector<double>& pressure,
            std::vector<double>& out_slope) const {
        out_slope.resize(pressure.size());
        for (std::size_t i = 0; i < pressure.size(); ++i) {
            out_slope[i] = this->calc_melting_slope(pressure[i]);
        }
    }
};

// T0 (1 + (P - P_ref) / a)^(1 / c), held at T0 below P_ref and zero where the base reaches zero.
inline double c_simon_glatzel(
        double pressure,
        double temperature,
        double simon_a,
        double simon_c,
        double reference_pressure) noexcept {
    if (!std::isfinite(pressure)) { return TidalPyConstants::d_NAN; }
    const double base = 1.0 + (std::max(pressure, reference_pressure) - reference_pressure) / simon_a;
    if (!(base > 0.0)) { return 0.0; }
    return temperature * std::pow(base, 1.0 / simon_c);
}

// The slope of c_simon_glatzel, dT/dP [K Pa-1]: T0 / (a c) (1 + (P - P_ref) / a)^(1 / c - 1), zero below P_ref and
// where the base reaches zero, as the curve is held there.
inline double c_simon_glatzel_slope(
        double pressure,
        double temperature,
        double simon_a,
        double simon_c,
        double reference_pressure) noexcept {
    if (!std::isfinite(pressure)) { return TidalPyConstants::d_NAN; }
    if (!(pressure > reference_pressure)) { return 0.0; }
    const double base = 1.0 + (pressure - reference_pressure) / simon_a;
    if (!(base > 0.0)) { return 0.0; }
    return temperature / (simon_a * simon_c) * std::pow(base, 1.0 / simon_c - 1.0);
}

// A melting temperature independent of pressure (alias "const").
class c_ConstantMeltingCurve final : public c_SpecModel<c_ConstantMeltingCurve, c_MeltingCurveBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ConstantMeltingCurve;

    static const std::vector<c_ParamSpec<c_ConstantMeltingCurve>>& parameter_specs() {
        using Self = c_ConstantMeltingCurve;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"temperature", "temperature_k", &Self::p_temperature, 1600.0, c_ParamBounds::Positive,
             "Melting temperature [K]."},
        };
        return specs;
    }

    c_ConstantMeltingCurve() : c_ConstantMeltingCurve(c_ParamMap{}) {}
    explicit c_ConstantMeltingCurve(const c_ParamMap& params) : c_SpecModel("constant") {
        this->p_initialize(params);
    }

    double calc_melting_temperature(double pressure) const noexcept override {
        return std::isfinite(pressure) ? this->p_temperature : TidalPyConstants::d_NAN;
    }
    double calc_melting_slope(double pressure) const noexcept override {
        return std::isfinite(pressure) ? 0.0 : TidalPyConstants::d_NAN;
    }

protected:
    double p_temperature = 0.0;
};

// One Simon and Glatzel branch.
class c_SimonGlatzelCurve final : public c_SpecModel<c_SimonGlatzelCurve, c_MeltingCurveBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::SimonGlatzelCurve;

    static const std::vector<c_ParamSpec<c_SimonGlatzelCurve>>& parameter_specs() {
        using Self = c_SimonGlatzelCurve;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"temperature", "temperature_k", &Self::p_temperature, 1600.0, c_ParamBounds::Positive,
             "Melting temperature at the reference pressure, T0 [K]."},
            {"simon_a", "simon_a_pa", &Self::p_simon_a, 1.0e9, c_ParamBounds::Finite,
             "Simon-Glatzel a [Pa]; negative for a curve that falls with pressure."},
            {"simon_c", "simon_c", &Self::p_simon_c, 5.0, c_ParamBounds::Positive,
             "Simon-Glatzel c [dimensionless]."},
            {"reference_pressure", "reference_pressure_pa", &Self::p_reference_pressure, 0.0, c_ParamBounds::Finite,
             "Pressure where T0 applies [Pa]."},
        };
        return specs;
    }

    c_SimonGlatzelCurve() : c_SimonGlatzelCurve(c_ParamMap{}) {}
    explicit c_SimonGlatzelCurve(const c_ParamMap& params) : c_SpecModel("simon_glatzel") {
        this->p_initialize(params);
    }

    double calc_melting_temperature(double pressure) const noexcept override {
        return c_simon_glatzel(
            pressure, this->p_temperature, this->p_simon_a, this->p_simon_c, this->p_reference_pressure);
    }
    double calc_melting_slope(double pressure) const noexcept override {
        return c_simon_glatzel_slope(
            pressure, this->p_temperature, this->p_simon_a, this->p_simon_c, this->p_reference_pressure);
    }

protected:
    void p_validate() const override {
        if (this->p_simon_a == 0.0) {
            throw std::invalid_argument(this->p_describe() + " needs a nonzero 'simon_a_pa'.");
        }
    }

    double p_temperature        = 0.0;
    double p_simon_a            = 0.0;
    double p_simon_c            = 0.0;
    double p_reference_pressure = 0.0;
};

// Two Simon and Glatzel branches joined at a transition pressure: the low branch below it, the high branch, written
// in the absolute pressure as the published fits are, above it.
class c_SimonGlatzel2Curve final : public c_SpecModel<c_SimonGlatzel2Curve, c_MeltingCurveBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::SimonGlatzel2Curve;

    static const std::vector<c_ParamSpec<c_SimonGlatzel2Curve>>& parameter_specs() {
        using Self = c_SimonGlatzel2Curve;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"temperature", "temperature_k", &Self::p_temperature, 1661.2, c_ParamBounds::Positive,
             "Low branch: melting temperature at the reference pressure, T0 [K]."},
            {"simon_a", "simon_a_pa", &Self::p_simon_a, 1.336e9, c_ParamBounds::Finite, "Low branch: a [Pa]."},
            {"simon_c", "simon_c", &Self::p_simon_c, 7.437, c_ParamBounds::Positive, "Low branch: c [dimensionless]."},
            {"reference_pressure", "reference_pressure_pa", &Self::p_reference_pressure, 0.0, c_ParamBounds::Finite,
             "Low branch: pressure where T0 applies [Pa]."},
            {"transition_pressure", "transition_pressure_pa", &Self::p_transition_pressure, 20.0e9,
             c_ParamBounds::Finite, "Pressure above which the high branch applies [Pa]."},
            {"high_temperature", "high_temperature_k", &Self::p_high_temperature, 2081.8, c_ParamBounds::Positive,
             "High branch: T0 [K]."},
            {"high_simon_a", "high_simon_a_pa", &Self::p_high_simon_a, 1.0169e11, c_ParamBounds::Finite,
             "High branch: a [Pa]."},
            {"high_simon_c", "high_simon_c", &Self::p_high_simon_c, 1.226, c_ParamBounds::Positive,
             "High branch: c [dimensionless]."},
            {"high_reference_pressure", "high_reference_pressure_pa", &Self::p_high_reference_pressure, 0.0,
             c_ParamBounds::Finite, "High branch: pressure where its T0 applies [Pa]."},
        };
        return specs;
    }

    c_SimonGlatzel2Curve() : c_SimonGlatzel2Curve(c_ParamMap{}) {}
    explicit c_SimonGlatzel2Curve(const c_ParamMap& params) : c_SpecModel("simon_glatzel_2") {
        this->p_initialize(params);
    }

    double calc_melting_temperature(double pressure) const noexcept override {
        if (std::isfinite(pressure) && (pressure > this->p_transition_pressure)) {
            return c_simon_glatzel(
                pressure, this->p_high_temperature, this->p_high_simon_a, this->p_high_simon_c,
                this->p_high_reference_pressure);
        }
        return c_simon_glatzel(
            pressure, this->p_temperature, this->p_simon_a, this->p_simon_c, this->p_reference_pressure);
    }
    double calc_melting_slope(double pressure) const noexcept override {
        if (std::isfinite(pressure) && (pressure > this->p_transition_pressure)) {
            return c_simon_glatzel_slope(
                pressure, this->p_high_temperature, this->p_high_simon_a, this->p_high_simon_c,
                this->p_high_reference_pressure);
        }
        return c_simon_glatzel_slope(
            pressure, this->p_temperature, this->p_simon_a, this->p_simon_c, this->p_reference_pressure);
    }

protected:
    void p_validate() const override {
        if ((this->p_simon_a == 0.0) || (this->p_high_simon_a == 0.0)) {
            throw std::invalid_argument(this->p_describe() + " needs nonzero 'simon_a_pa' and 'high_simon_a_pa'.");
        }
    }

    double p_temperature             = 0.0;
    double p_simon_a                 = 0.0;
    double p_simon_c                 = 0.0;
    double p_reference_pressure      = 0.0;
    double p_transition_pressure     = 0.0;
    double p_high_temperature        = 0.0;
    double p_high_simon_a            = 0.0;
    double p_high_simon_c            = 0.0;
    double p_high_reference_pressure = 0.0;
};

// A melting curve tabulated in pressure (aliases "interp", "interpolated"), linear between the points and held at the
// end values beyond them.
class c_InterpolatedMeltingCurve final : public c_SpecModel<c_InterpolatedMeltingCurve, c_MeltingCurveBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::InterpolatedMeltingCurve;

    static const std::vector<c_ParamSpec<c_InterpolatedMeltingCurve>>& parameter_specs() {
        using Self = c_InterpolatedMeltingCurve;
        static const std::vector<c_ParamSpec<Self>> specs = {
            // The default tables hold one point, a constant curve.
            {"pressure", "pressure_pa", &Self::p_pressure, 0.0, c_ParamBounds::Finite,
             "Table pressures, ascending [Pa].", {0.0}},
            {"temperature", "temperature_k", &Self::p_temperature, 0.0, c_ParamBounds::Positive,
             "Melting temperature at each table pressure [K].", {1600.0}},
        };
        return specs;
    }

    c_InterpolatedMeltingCurve() : c_InterpolatedMeltingCurve(c_ParamMap{}) {}
    explicit c_InterpolatedMeltingCurve(const c_ParamMap& params) : c_SpecModel("interpolate") {
        this->p_initialize(params);
    }

    double calc_melting_temperature(double pressure) const noexcept override {
        if (!std::isfinite(pressure)) { return TidalPyConstants::d_NAN; }
        return this->p_lookup.interpolate(pressure, this->p_pressure, this->p_temperature);
    }
    double calc_melting_slope(double pressure) const noexcept override {
        if (!std::isfinite(pressure)) { return TidalPyConstants::d_NAN; }
        return this->p_lookup.slope(pressure, this->p_pressure, this->p_temperature);
    }

protected:
    void p_validate() const override {
        c_check_table(this->p_describe(), this->p_pressure, {&this->p_temperature});
        if (this->p_temperature.size() != this->p_pressure.size()) {
            throw std::invalid_argument(this->p_describe() + " needs a 'temperature_k' value at every table pressure.");
        }
    }
    void p_update_derived() noexcept override { this->p_lookup.build(this->p_pressure); }

    std::vector<double> p_pressure;
    std::vector<double> p_temperature;
    c_TableLookup p_lookup;
};

inline const c_ModelRegistry<c_MeltingCurveBase>& c_melting_curve_registry() {
    static const c_ModelRegistry<c_MeltingCurveBase> registry = {
        {{"constant", "const"}, BinaryClassID::ConstantMeltingCurve,
         &c_make_entry<c_MeltingCurveBase, c_ConstantMeltingCurve>},
        {{"simon_glatzel", "simon-glatzel"}, BinaryClassID::SimonGlatzelCurve,
         &c_make_entry<c_MeltingCurveBase, c_SimonGlatzelCurve>},
        {{"simon_glatzel_2", "simon-glatzel-2"}, BinaryClassID::SimonGlatzel2Curve,
         &c_make_entry<c_MeltingCurveBase, c_SimonGlatzel2Curve>},
        {{"interpolate", "interp", "interpolated"}, BinaryClassID::InterpolatedMeltingCurve,
         &c_make_entry<c_MeltingCurveBase, c_InterpolatedMeltingCurve>},
    };
    return registry;
}

inline std::unique_ptr<c_MeltingCurveBase> c_find_melting_curve(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_melting_curve_registry(), model_name, params);
}

inline std::unique_ptr<c_MeltingCurveBase> c_melting_curve_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_melting_curve_registry(), in, force);
}

inline std::string c_melting_curve_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_melting_curve_registry(), model_name);
}

inline std::vector<std::string> c_melting_curve_model_names() {
    return c_model_names(c_melting_curve_registry());
}

}  // namespace tidalpy
