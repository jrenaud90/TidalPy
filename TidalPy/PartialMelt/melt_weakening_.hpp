#pragma once
/* Melt weakening: how a partially molten aggregate's shear modulus [Pa] and viscosity [Pa s] fall from the solid
 * phase's values toward the liquid phase's as the melt fraction rises. All MKS.
 *
 * A model sees the local temperature, the solidus and liquidus at the local pressure, the melt fraction, and both
 * phases' values, and returns the aggregate's pair, floored at the liquid's. Every model returns the solid pair with
 * no melt and the liquid pair when fully molten. The material decides, from the shear modulus a model returns,
 * where the radial solver treats the aggregate as a liquid.
 *
 * The weakening laws are continuous in the melt fraction, so an integration through the melting range sees no jump.
 * Each law describes the partially molten framework up to a rheological transition, the breakdown band
 * [crit, crit + width] of melt fraction, across which the framework's pair blends into the liquid's (the viscosity
 * log-linearly, the shear modulus linearly, since the liquid's may be zero); past the band the aggregate is the liquid.
 * A zero width makes the transition a step.
 *
 * Both temperature laws are anchored at the solidus, so they carry over to materials whose solidus is not the 1600 K
 * silicate value the published fits assume; the published forms are recovered at T_sol = 1600 K.
 *
 * References
 * ----------
 * - Fischer and Spohn (1990), Icarus 83, 39.
 * - Henning, O'Connell, and Sasselov (2009), ApJ 707, 1000; Renaud and Henning (2018), ApJ 857, 98.
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
#include "../Utilities/math/numerics_.hpp"  // c_safe_pow, c_safe_exp

namespace tidalpy {

// Everything a weakening model sees at a point.
struct c_MeltWeakeningInputs {
    double temperature      = TidalPyConstants::d_NAN;   // [K]
    double solidus          = TidalPyConstants::d_NAN;   // at the local pressure [K]
    double liquidus         = TidalPyConstants::d_NAN;   // at the local pressure [K]
    double melt_fraction    = 0.0;                       // [m3 m-3]
    double solid_shear      = TidalPyConstants::d_NAN;   // [Pa]
    double solid_viscosity  = TidalPyConstants::d_NAN;   // [Pa s]
    double liquid_shear     = 0.0;                       // [Pa]
    double liquid_viscosity = TidalPyConstants::d_NAN;   // [Pa s]
    // The solid phase's pair at the solidus and the local pressure, for a law anchored there
    // (get_uses_solidus_state); NaN otherwise.
    double solid_shear_at_solidus     = TidalPyConstants::d_NAN;   // [Pa]
    double solid_viscosity_at_solidus = TidalPyConstants::d_NAN;   // [Pa s]
};

struct c_MeltWeakeningResult {
    double shear_modulus = TidalPyConstants::d_NAN;   // [Pa]
    double viscosity     = TidalPyConstants::d_NAN;   // [Pa s]
};

class c_MeltWeakeningBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "melt weakening";

    explicit c_MeltWeakeningBase(const std::string& model_name) : c_PhysicsBase(model_name) {}
    ~c_MeltWeakeningBase() override = default;

    // The aggregate's shear modulus and viscosity. NaN for a non-finite melt fraction.
    c_MeltWeakeningResult calc_weakening(const c_MeltWeakeningInputs& inputs) const noexcept {
        c_MeltWeakeningResult result;
        const double phi = inputs.melt_fraction;
        if (!std::isfinite(phi)) { return result; }
        if (phi <= 0.0) {
            result.shear_modulus = inputs.solid_shear;
            result.viscosity     = inputs.solid_viscosity;
            return result;
        }
        const double band_start = this->p_crit_melt_frac;
        const double band_end = std::min(band_start + this->p_crit_melt_frac_width, 1.0);
        if (phi >= band_end) {
            result.shear_modulus = inputs.liquid_shear;
            result.viscosity     = inputs.liquid_viscosity;
            return result;
        }
        this->p_calc_partial(inputs, result);
        // Never weaker than the liquid.
        if (!(result.shear_modulus >= inputs.liquid_shear)) { result.shear_modulus = inputs.liquid_shear; }
        if (!(result.viscosity >= inputs.liquid_viscosity))  { result.viscosity     = inputs.liquid_viscosity; }
        // Across the breakdown band the framework's pair blends into the liquid's, reaching it at the band's end.
        if (phi > band_start) {
            const double blend = (phi - band_start) / (band_end - band_start);
            result.shear_modulus = (1.0 - blend) * result.shear_modulus + blend * inputs.liquid_shear;
            if ((result.viscosity > 0.0) && (inputs.liquid_viscosity > 0.0) && std::isfinite(result.viscosity)) {
                result.viscosity = std::exp(
                    (1.0 - blend) * std::log(result.viscosity) + blend * std::log(inputs.liquid_viscosity));
            }
        }
        return result;
    }

    // Whether the law reads the solid's pair at the solidus (c_MeltWeakeningInputs::solid_*_at_solidus).
    virtual bool get_uses_solidus_state() const noexcept { return false; }

protected:
    // The framework's pair inside the melting range, 0 < phi < crit + width.
    virtual void p_calc_partial(const c_MeltWeakeningInputs& inputs, c_MeltWeakeningResult& out) const noexcept = 0;

    // The breakdown band's parameters, as rows of a model's parameter table.
    template <class Model>
    static void p_append_breakdown_specs(std::vector<c_ParamSpec<Model>>& specs) {
        specs.push_back({"crit_melt_frac", "crit_melt_frac", &Model::p_crit_melt_frac, 0.5,
                         c_ParamBounds::UnitInterval, "Critical melt fraction, where the breakdown starts [m3 m-3]."});
        specs.push_back({"crit_melt_frac_width", "crit_melt_frac_width", &Model::p_crit_melt_frac_width, 0.05,
                         c_ParamBounds::UnitInterval,
                         "Width of the breakdown band, across which the aggregate becomes the liquid [m3 m-3]."});
    }

    // The breakdown band [crit, crit + width] of melt fraction; a law without one keeps its framework to phi = 1.
    double p_crit_melt_frac       = 1.0;
    double p_crit_melt_frac_width = 0.0;
};

// No weakening (alias "off"): the solid's values until the material is fully molten, then the liquid's.
class c_NoMeltWeakening final : public c_SpecModel<c_NoMeltWeakening, c_MeltWeakeningBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::NoMeltWeakening;
    static const std::vector<c_ParamSpec<c_NoMeltWeakening>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_NoMeltWeakening>> specs;
        return specs;
    }

    c_NoMeltWeakening() : c_NoMeltWeakening(c_ParamMap{}) {}
    explicit c_NoMeltWeakening(const c_ParamMap& params) : c_SpecModel("none") { this->p_initialize(params); }

protected:
    void p_calc_partial(const c_MeltWeakeningInputs& inputs, c_MeltWeakeningResult& out) const noexcept override {
        out.shear_modulus = inputs.solid_shear;
        out.viscosity     = inputs.solid_viscosity;
    }
};

// Fischer and Spohn (1990) temperature law (aliases "fischer", "fischer_spohn"): above the solidus both strengths are
// 10^(log10_at_solidus + slope (1 / T - 1 / T_sol)). By default log10_at_solidus is the solid's own value at the
// solidus and local pressure, so the law continues the solid's; a finite value gives the absolute fits,
// 10^(27000 / T - 1) Pa s and 10^(82000 / T - 40.6) Pa, at T_sol = 1600 K with 15.875 and 10.65 (and steps from the
// solid's value at the solidus to it).
class c_SpohnMeltWeakening final : public c_SpecModel<c_SpohnMeltWeakening, c_MeltWeakeningBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::SpohnMeltWeakening;

    static const std::vector<c_ParamSpec<c_SpohnMeltWeakening>>& parameter_specs() {
        using Self = c_SpohnMeltWeakening;
        static const std::vector<c_ParamSpec<Self>> specs = [] {
            std::vector<c_ParamSpec<Self>> rows = {
                {"visc_power_slope", "fs_visc_power_slope_k", &Self::p_visc_power_slope, 27000.0,
                 c_ParamBounds::Finite, "Viscosity slope s in 10^(L + s (1/T - 1/T_sol)) [K]."},
                {"visc_log10_at_solidus", "fs_visc_log10_at_solidus", &Self::p_visc_log10_at_solidus,
                 TidalPyConstants::d_NAN, c_ParamBounds::Any,
                 "log10 of the viscosity at the solidus [log10 Pa s]; unset takes the solid's own."},
                {"shear_power_slope", "fs_shear_power_slope_k", &Self::p_shear_power_slope, 82000.0,
                 c_ParamBounds::Finite, "Shear-modulus slope s in 10^(L + s (1/T - 1/T_sol)) [K]."},
                {"shear_log10_at_solidus", "fs_shear_log10_at_solidus", &Self::p_shear_log10_at_solidus,
                 TidalPyConstants::d_NAN, c_ParamBounds::Any,
                 "log10 of the shear modulus at the solidus [log10 Pa]; unset takes the solid's own."},
            };
            p_append_breakdown_specs(rows);
            return rows;
        }();
        return specs;
    }

    c_SpohnMeltWeakening() : c_SpohnMeltWeakening(c_ParamMap{}) {}
    explicit c_SpohnMeltWeakening(const c_ParamMap& params) : c_SpecModel("spohn") { this->p_initialize(params); }

    bool get_uses_solidus_state() const noexcept override {
        return !std::isfinite(this->p_visc_log10_at_solidus) || !std::isfinite(this->p_shear_log10_at_solidus);
    }

protected:
    void p_calc_partial(const c_MeltWeakeningInputs& inputs, c_MeltWeakeningResult& out) const noexcept override {
        const double inv_temp_shift = (1.0 / inputs.temperature) - (1.0 / inputs.solidus);
        const double visc_log10_at_solidus = std::isfinite(this->p_visc_log10_at_solidus)
            ? this->p_visc_log10_at_solidus : std::log10(inputs.solid_viscosity_at_solidus);
        const double shear_log10_at_solidus = std::isfinite(this->p_shear_log10_at_solidus)
            ? this->p_shear_log10_at_solidus : std::log10(inputs.solid_shear_at_solidus);
        out.viscosity     = c_safe_pow(10.0, visc_log10_at_solidus + this->p_visc_power_slope * inv_temp_shift);
        out.shear_modulus = c_safe_pow(10.0, shear_log10_at_solidus + this->p_shear_power_slope * inv_temp_shift);
    }

    double p_visc_power_slope       = 0.0;
    double p_visc_log10_at_solidus  = 0.0;
    double p_shear_power_slope      = 0.0;
    double p_shear_log10_at_solidus = 0.0;
};

// Henning (2009, 2010) three-regime weakening: exponential below the critical melt fraction, a steeper falloff across
// the breakdown band [crit, crit + width] that blends into the liquid's values by the band's end, then the liquid's.
// The sub-critical shear law mu_s exp[b1 (1/T - 1/T_sol)] is 1 at the solidus; Henning et al. (2009) Eq. 20,
// exp(40000 / T - 25), is this form for b1 = 40000 K at T_sol = 1600 K.
class c_HenningMeltWeakening final : public c_SpecModel<c_HenningMeltWeakening, c_MeltWeakeningBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::HenningMeltWeakening;

    static const std::vector<c_ParamSpec<c_HenningMeltWeakening>>& parameter_specs() {
        using Self = c_HenningMeltWeakening;
        static const std::vector<c_ParamSpec<Self>> specs = [] {
            std::vector<c_ParamSpec<Self>> rows;
            p_append_breakdown_specs(rows);
            rows.insert(rows.end(), {
                {"visc_slope_1", "hn_visc_slope_1", &Self::p_visc_slope_1, 13.5, c_ParamBounds::NonNegative,
                 "Sub-critical viscosity slope in exp(-slope phi) [dimensionless]."},
                {"visc_falloff_slope", "hn_visc_falloff_slope", &Self::p_visc_falloff_slope, 370.0,
                 c_ParamBounds::NonNegative, "Breakdown-band viscosity slope [dimensionless]."},
                {"shear_param_1", "hn_shear_param_1_k", &Self::p_shear_param_1, 40000.0, c_ParamBounds::NonNegative,
                 "b1 in the sub-critical shear law mu_s exp[b1 (1/T - 1/T_sol)] [K]."},
                {"shear_falloff_slope", "hn_shear_falloff_slope", &Self::p_shear_falloff_slope, 700.0,
                 c_ParamBounds::NonNegative, "Breakdown-band shear slope [dimensionless]."},
            });
            return rows;
        }();
        return specs;
    }

    c_HenningMeltWeakening() : c_HenningMeltWeakening(c_ParamMap{}) {}
    explicit c_HenningMeltWeakening(const c_ParamMap& params) : c_SpecModel("henning") { this->p_initialize(params); }

protected:
    void p_calc_partial(const c_MeltWeakeningInputs& inputs, c_MeltWeakeningResult& out) const noexcept override {
        const double phi        = inputs.melt_fraction;
        const double crit       = this->p_crit_melt_frac;
        const double inv_solidus = 1.0 / inputs.solidus;
        if (phi < crit) {
            out.viscosity     = inputs.solid_viscosity * c_safe_exp(-this->p_visc_slope_1 * phi);
            out.shear_modulus = inputs.solid_shear
                              * c_safe_exp(this->p_shear_param_1 * ((1.0 / inputs.temperature) - inv_solidus));
        } else {
            // The full sub-critical effect, then a steep falloff (blended into the liquid by the base).
            const double break_temperature = inputs.solidus + crit * (inputs.liquidus - inputs.solidus);
            out.viscosity     = inputs.solid_viscosity * c_safe_exp(-this->p_visc_slope_1 * crit)
                              * c_safe_exp(-this->p_visc_falloff_slope * (phi - crit));
            out.shear_modulus = inputs.solid_shear
                              * c_safe_exp(this->p_shear_param_1 * ((1.0 / break_temperature) - inv_solidus))
                              * c_safe_exp(-this->p_shear_falloff_slope * (phi - crit));
        }
    }

    double p_visc_slope_1         = 0.0;
    double p_visc_falloff_slope   = 0.0;
    double p_shear_param_1        = 0.0;
    double p_shear_falloff_slope  = 0.0;
};

inline const c_ModelRegistry<c_MeltWeakeningBase>& c_melt_weakening_registry() {
    static const c_ModelRegistry<c_MeltWeakeningBase> registry = {
        {{"none", "off"}, BinaryClassID::NoMeltWeakening, &c_make_entry<c_MeltWeakeningBase, c_NoMeltWeakening>},
        {{"spohn", "fischer", "fischer_spohn"}, BinaryClassID::SpohnMeltWeakening,
         &c_make_entry<c_MeltWeakeningBase, c_SpohnMeltWeakening>},
        {{"henning"}, BinaryClassID::HenningMeltWeakening,
         &c_make_entry<c_MeltWeakeningBase, c_HenningMeltWeakening>},
    };
    return registry;
}

inline std::unique_ptr<c_MeltWeakeningBase> c_find_melt_weakening(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_melt_weakening_registry(), model_name, params);
}

inline std::unique_ptr<c_MeltWeakeningBase> c_melt_weakening_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_melt_weakening_registry(), in, force);
}

inline std::string c_melt_weakening_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_melt_weakening_registry(), model_name);
}

inline std::vector<std::string> c_melt_weakening_model_names() {
    return c_model_names(c_melt_weakening_registry());
}

}  // namespace tidalpy
