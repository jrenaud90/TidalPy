#pragma once
/* TidalPy's radiogenic heating models. All MKS.
 *
 * Each model declares its parameters in one table (c_SpecModel, spec_model_.hpp), which gives its constructor,
 * validation, config entries, and binary record. The isotope and fixed models share one reference time, so times can
 * be measured from any fixed epoch such as solar-system formation.
 *
 * References
 * ----------
 * - Hussmann and Spohn (2004); Turcotte and Schubert (2001): chondritic isotope data.
 * - Castillo-Rogez et al. (2007): long- and short-lived radiogenic isotopes.
 * - McDonough and Sun (1995): bulk silicate Earth elemental abundances.
 */

#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "binary_.hpp"
#include "radiogenics_base_.hpp"
#include "../Utilities/math/numerics_.hpp"  // c_safe_exp
#include "constants_.hpp"
#include "model_names_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"

namespace tidalpy {

// One isotope of a built-in dataset: its label and the four values the isotope model takes for it (see
// c_IsotopeRadiogenics).
struct c_Isotope {
    std::string name;                  // isotope label (e.g. "U238")
    double heat_production = 0.0;      // specific heat production of the pure isotope [W/kg]
    double half_life       = 0.0;      // half life [s]
    double mass_frac       = 0.0;      // isotopic mass fraction within its element at the reference time [kg/kg]
    double concentration   = 0.0;      // element concentration in the layer material at the reference time [kg/kg]

    c_Isotope() = default;
    c_Isotope(std::string isotope_name,
              double hpr,
              double half_life_seconds,
              double isotopic_mass_frac,
              double element_concentration)
        : name(std::move(isotope_name)),
          heat_production(hpr),
          half_life(half_life_seconds),
          mass_frac(isotopic_mass_frac),
          concentration(element_concentration) {}
};

// Built-in catalogs of literature isotope sets, so a caller need not hand-enter abundances. The source
// literature quotes half lives and reference times in Myr; they are converted to seconds here.
struct c_IsotopeDataset {
    std::vector<c_Isotope> isotopes;
    double ref_time = 0.0;
};

inline std::vector<std::string> c_isotope_dataset_names() {
    return {"modern_day_chondritic", "llri", "slri", "llri_and_slri", "bulk_silicate_earth"};
}

// Named built-in isotope datasets (case-insensitive):
//
//   "modern_day_chondritic"
//       Present-day chondritic abundances of the four long-lived heat producers (U238, U235, Th232,
//       K40), for rocky or icy bodies of broadly chondritic composition near the present epoch.
//       Hussmann and Spohn (2004); Turcotte and Schubert (2001).
//
//   "llri"
//       The long-lived isotopes (U238, U235, Th232, K40) of ordinary chondritic rock at CAI formation, for a body
//       that formed too late for the short-lived isotopes to matter. Castillo-Rogez et al. (2007).
//
//   "slri"
//       The short-lived isotopes (Al26, Fe60, Mn53) of the same rock at CAI formation. Castillo-Rogez et al. (2007).
//
//   "llri_and_slri"
//       Both of the above, for early-solar-system thermal evolution where the short-lived isotopes dominate the
//       heat budget.
//
//   "bulk_silicate_earth"
//       Present-day bulk silicate Earth: the four long-lived heat producers at BSE concentrations
//       (U = 20.3 ppb, Th = 79.5 ppb, K = 240 ppm), for Earth-like silicate mantles.
//       McDonough and Sun (1995) for the concentrations; Turcotte and Schubert (2002) for the heat
//       production rates and half lives.
//
// The chondritic and BSE sets quote present-epoch concentrations (ref_time = 4600 Myr after
// solar-system formation); the three Castillo-Rogez sets quote formation (CAI) abundances (ref_time = 0), so
// time is measured from formation and the short-lived isotopes decay away over the first 10 Myr or so.
inline c_IsotopeDataset c_get_isotope_dataset(const std::string& name);

// Castillo-Rogez et al. (2007) long-lived isotopes at CAI formation: heat production and half lives from their
// Table 4, concentrations from their Table 3. Table 3 quotes each isotope's own concentration, already decayed back
// to formation, while the isotopic abundances of Table 4 are present-day values that do not hold at formation, so
// each isotope is entered with its own concentration and a mass fraction of 1. The Th232 half life is the middle
// of the quoted 14010 to 14050 Myr. The K40 half life is 1248 Myr, not the table's 1277 Myr, which its own decay
// constant (5.54e-10 per year) contradicts. Table 4 labels heat production per kg of element; the values are per kg
// of isotope, as used here.
inline std::vector<c_Isotope> c_castillo_rogez_2007_llri() {
    const double myr = TidalPyConstants::d_SECONDS_PER_MYR;
    return {
        c_Isotope("U238",  9.465e-5, 4468.0  * myr, 1.0, 26.2e-9),
        c_Isotope("U235",  5.687e-4, 703.81  * myr, 1.0, 8.2e-9),
        c_Isotope("Th232", 2.638e-5, 14030.0 * myr, 1.0, 53.8e-9),
        c_Isotope("K40",   2.917e-5, 1248.0  * myr, 1.0, 1104.0e-9),
    };
}

// Castillo-Rogez et al. (2007) short-lived isotopes at CAI formation: heat production, half lives, and initial
// isotopic ratios from their Table 5. Each element concentration is the Table 3 isotope concentration divided by
// that ratio, as the paper builds Table 3 (26Al: 5e-5 of 1.2 wt% aluminum is 600 ppb), so the product reproduces
// Table 3 exactly. The ratio to the stable isotope stands in for the mass fraction within the element, as it does
// in the paper; the Fe60 ratio is the 1e-6 of the paper's short-lived-isotope models. The decay data are corrected
// from Table 5, whose values are inconsistent with the decay energies: Al26 deposits 3.12 MeV per decay, 0.355 W/kg
// at a 0.717 Myr half life (Lebrun et al. 2013, after Sramek et al. 2012; Table 5's 0.146 W/kg implies 1.29 MeV);
// Fe60 has a 2.62 Myr half life (Rugel et al. 2009, PRL 103, 072502) and deposits 2.71 MeV with its Co60 daughter,
// 0.0366 W/kg; Mn53 decays by electron capture and deposits only its X-ray and Auger energy, about 5 keV, 5.8e-5 W/kg
// (Table 5's 0.027 W/kg exceeds its whole 0.597 MeV decay energy).
inline std::vector<c_Isotope> c_castillo_rogez_2007_slri() {
    const double myr = TidalPyConstants::d_SECONDS_PER_MYR;
    return {
        c_Isotope("Al26", 0.355,  0.717 * myr, 5.0e-5, 1.2e-2),
        c_Isotope("Fe60", 0.0366, 2.62  * myr, 1.0e-6, 0.225),
        c_Isotope("Mn53", 5.8e-5, 3.7   * myr, 1.0e-5, 2.57e-3),
    };
}

inline c_IsotopeDataset c_get_isotope_dataset(const std::string& name) {
    const std::string key = c_to_lower(name);
    // The datasets quote half lives and reference times in Myr; convert so the C++ API stays MKS.
    const double myr = TidalPyConstants::d_SECONDS_PER_MYR;
    c_IsotopeDataset dataset;
    dataset.ref_time = 4600.0 * myr;

    if (key == "modern_day_chondritic") {
        // Hussmann and Spohn (2004); Turcotte and Schubert (2001).
        dataset.isotopes = {
            c_Isotope("U238",  9.48e-5, 4470.0  * myr, 0.9928,   0.012e-6),
            c_Isotope("U235",  5.69e-4, 704.0   * myr, 0.0071,   0.012e-6),
            c_Isotope("Th232", 2.69e-5, 14000.0 * myr, 0.9998,   0.04e-6),
            c_Isotope("K40",   2.92e-5, 1250.0  * myr, 1.19e-4,  840.0e-6),
        };
        return dataset;
    }
    if ((key == "llri") || (key == "slri") || (key == "llri_and_slri")) {
        // Castillo-Rogez et al. (2007) quotes formation (CAI) abundances, so ref_time is formation.
        dataset.ref_time = 0.0;
        if (key != "slri") {
            dataset.isotopes = c_castillo_rogez_2007_llri();
        }
        if (key != "llri") {
            const std::vector<c_Isotope> short_lived = c_castillo_rogez_2007_slri();
            dataset.isotopes.insert(dataset.isotopes.end(), short_lived.begin(), short_lived.end());
        }
        return dataset;
    }
    if (key == "bulk_silicate_earth") {
        // McDonough and Sun (1995) concentrations (U 20.3 ppb, Th 79.5 ppb, K 240 ppm);
        // Turcotte and Schubert (2002) heat production rates and half lives.
        dataset.isotopes = {
            c_Isotope("U238",  9.48e-5, 4470.0  * myr, 0.9928,   20.3e-9),
            c_Isotope("U235",  5.69e-4, 704.0   * myr, 0.0071,   20.3e-9),
            c_Isotope("Th232", 2.69e-5, 14000.0 * myr, 0.9998,   79.5e-9),
            c_Isotope("K40",   2.92e-5, 1250.0  * myr, 1.19e-4,  240.0e-6),
        };
        return dataset;
    }

    throw std::invalid_argument("TidalPy: unknown isotope dataset '" + name + "'");
}

// Radiogenics disabled (alias "none").
class c_OffRadiogenics final : public c_SpecModel<c_OffRadiogenics, c_RadiogenicsBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::OffRadiogenics;

    static const std::vector<c_ParamSpec<c_OffRadiogenics>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_OffRadiogenics>> specs = {};
        return specs;
    }

    c_OffRadiogenics() : c_OffRadiogenics(c_ParamMap{}) {}
    explicit c_OffRadiogenics(const c_ParamMap& params) : c_SpecModel("off") { this->p_initialize(params); }

    double calc_heating(double /*time*/, double /*mass*/) const override { return 0.0; }
};

// The sum over individually decaying isotopes (alias "isotopes"). Isotope i heats a unit of layer mass at
//
//   q_i(t) = mass_frac_i * concentration_i * heat_production_i * exp(gamma_i * (t - ref_time))
//
// with gamma_i = ln(0.5) / half_life_i the (negative) decay constant. Both abundances are values at ref_time. A source
// that quotes the isotope's own concentration rather than its element's is entered with a mass fraction of 1. The four
// tables hold one value per isotope. The isotopes may carry labels (isotope_names), none or one per isotope; they name
// the isotopes in the config and the binary record and take no part in the heating.
class c_IsotopeRadiogenics final : public c_SpecModel<c_IsotopeRadiogenics, c_RadiogenicsBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::IsotopeRadiogenics;

    static const std::vector<c_ParamSpec<c_IsotopeRadiogenics>>& parameter_specs() {
        using Self = c_IsotopeRadiogenics;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"heat_production", "heat_production_w_kg", &Self::p_heat_production, 0.0, c_ParamBounds::NonNegative,
             "Specific heat production of each pure isotope [W/kg]."},
            {"half_lives", "half_lives_s", &Self::p_half_lives, 0.0, c_ParamBounds::PositiveOrInfinite,
             "Half life of each isotope [s]; infinite for a stable one."},
            {"mass_fracs", "mass_fracs", &Self::p_mass_fracs, 0.0, c_ParamBounds::UnitInterval,
             "Mass fraction of each isotope within its element at the reference time [kg/kg]."},
            {"concentrations", "concentrations", &Self::p_concentrations, 0.0, c_ParamBounds::NonNegative,
             "Concentration of each isotope's element in the layer material at the reference time [kg/kg]."},
            {"ref_time", "ref_time_s", &Self::p_ref_time, 0.0, c_ParamBounds::Finite,
             "Time the abundances apply at [s], on the same clock as the heating's time."},
        };
        return specs;
    }

    c_IsotopeRadiogenics() : c_IsotopeRadiogenics(c_ParamMap{}) {}
    explicit c_IsotopeRadiogenics(const c_ParamMap& params) : c_IsotopeRadiogenics(params, {}) {}
    c_IsotopeRadiogenics(const c_ParamMap& params, std::vector<std::string> isotope_names) :
            c_SpecModel("isotope"),
            p_isotope_names(std::move(isotope_names)) {
        this->p_initialize(params);
    }

    double get_ref_time() const noexcept override { return this->p_ref_time; }
    std::size_t get_num_isotopes() const noexcept { return this->p_heat_production.size(); }
    const std::vector<std::string>& get_isotope_names() const noexcept { return this->p_isotope_names; }

    // The four tables always appear, empty ones too: a model given no isotopes is then rebuilt with none, where the
    // Python factory would otherwise take its default dataset for the missing tables.
    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RadiogenicsBase::append_config_entries(out);
        out.push_back(c_config_doubles("heat_production_w_kg", this->p_heat_production));
        out.push_back(c_config_doubles("half_lives_s", this->p_half_lives));
        out.push_back(c_config_doubles("mass_fracs", this->p_mass_fracs));
        out.push_back(c_config_doubles("concentrations", this->p_concentrations));
        out.push_back(c_config_double("ref_time_s", this->p_ref_time));
        if (!this->p_isotope_names.empty()) {
            out.push_back(c_config_strings("isotope_names", this->p_isotope_names));
        }
    }

    // Each isotope's specific heating, summed, then scaled by the layer mass.
    double calc_heating(double time, double mass) const override {
        double specific_heating = 0.0;
        for (std::size_t i = 0; i < this->p_reference_heating.size(); ++i) {
            specific_heating += this->p_reference_heating[i]
                              * c_safe_exp(this->p_decay_constants[i] * (time - this->p_ref_time));
        }
        return specific_heating * mass;
    }

    // Each point sums the same terms in the same order as calc_heating, then scales by its mass, so the two agree
    // exactly.
    void calc_heating_vectorize(
            const std::vector<double>& time,
            const std::vector<double>& mass,
            std::vector<double>& out_heating) const override {
        const std::size_t num_points = c_broadcast_length({time.size(), mass.size()}, "calc_heating_vectorize");
        const std::size_t time_stride = c_broadcast_stride(time.size());
        const std::size_t mass_stride = c_broadcast_stride(mass.size());
        out_heating.assign(num_points, 0.0);
        for (std::size_t isotope_i = 0; isotope_i < this->p_reference_heating.size(); ++isotope_i) {
            const double reference_heating = this->p_reference_heating[isotope_i];
            const double decay_constant    = this->p_decay_constants[isotope_i];
            for (std::size_t i = 0; i < num_points; ++i) {
                const double elapsed_time = time[i * time_stride] - this->p_ref_time;
                out_heating[i] += reference_heating * c_safe_exp(decay_constant * elapsed_time);
            }
        }
        for (std::size_t i = 0; i < num_points; ++i) {
            out_heating[i] *= mass[i * mass_stride];
        }
    }

protected:
    // One value per isotope in every table, and no labels or one per isotope.
    void p_validate() const override {
        const std::size_t num_isotopes = this->p_heat_production.size();
        if ((this->p_half_lives.size() != num_isotopes) || (this->p_mass_fracs.size() != num_isotopes) ||
            (this->p_concentrations.size() != num_isotopes)) {
            throw std::invalid_argument(
                this->p_describe() + " needs one value per isotope in each table; got "
                + std::to_string(num_isotopes) + " heat production rates, " + std::to_string(this->p_half_lives.size())
                + " half lives, " + std::to_string(this->p_mass_fracs.size()) + " mass fractions, and "
                + std::to_string(this->p_concentrations.size()) + " concentrations.");
        }
        if (!this->p_isotope_names.empty() && (this->p_isotope_names.size() != num_isotopes)) {
            throw std::invalid_argument(
                this->p_describe() + " has " + std::to_string(this->p_isotope_names.size()) + " isotope names for "
                + std::to_string(num_isotopes) + " isotopes; give one name per isotope or none.");
        }
    }

    // Each isotope's specific heating at the reference time [W/kg] and its decay constant [1/s].
    void p_update_derived() noexcept override {
        const std::size_t num_isotopes = this->p_heat_production.size();
        this->p_reference_heating.resize(num_isotopes);
        this->p_decay_constants.resize(num_isotopes);
        for (std::size_t i = 0; i < num_isotopes; ++i) {
            this->p_reference_heating[i] =
                this->p_mass_fracs[i] * this->p_concentrations[i] * this->p_heat_production[i];
            this->p_decay_constants[i] = TidalPyConstants::d_LN_HALF / c_guard_denominator(this->p_half_lives[i]);
        }
    }

    // The label count, then each label.
    void p_write_extra(std::ostream& out) const override {
        const uint64_t num_names = static_cast<uint64_t>(this->p_isotope_names.size());
        out.write(reinterpret_cast<const char*>(&num_names), sizeof(uint64_t));
        for (const std::string& name : this->p_isotope_names) { write_binary_string(out, name); }
    }

    void p_read_extra(std::istream& in) override {
        uint64_t num_names = 0;
        in.read(reinterpret_cast<char*>(&num_names), sizeof(uint64_t));
        if (!in) { throw std::runtime_error("TidalPy: failed to read isotope names from binary data"); }
        check_binary_count(in, num_names, sizeof(uint64_t), "isotope name");
        std::vector<std::string> names;
        names.reserve(static_cast<std::size_t>(num_names));
        for (uint64_t i = 0; i < num_names; ++i) { names.push_back(read_binary_string(in)); }
        this->p_isotope_names = std::move(names);
    }

    std::vector<double> p_heat_production;
    std::vector<double> p_half_lives;
    std::vector<double> p_mass_fracs;
    std::vector<double> p_concentrations;
    double p_ref_time = 0.0;
    std::vector<std::string> p_isotope_names;

    // From p_update_derived.
    std::vector<double> p_reference_heating;
    std::vector<double> p_decay_constants;
};

// One lumped rate with optional decay (alias "constant"): heating = mass * rate * exp(gamma (t - ref_time)), with
// gamma = ln(0.5) / average_half_life; an average half life of 0 or below means no decay.
class c_FixedRadiogenics final : public c_SpecModel<c_FixedRadiogenics, c_RadiogenicsBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::FixedRadiogenics;

    static const std::vector<c_ParamSpec<c_FixedRadiogenics>>& parameter_specs() {
        using Self = c_FixedRadiogenics;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"fixed_heat_production", "fixed_heat_production_w_kg", &Self::p_fixed_heat_production, 0.0,
             c_ParamBounds::NonNegative, "Specific heat production at the reference time [W/kg]."},
            {"average_half_life", "average_half_life_s", &Self::p_average_half_life, 0.0, c_ParamBounds::Finite,
             "Half life of the lumped rate's decay [s]; 0 or below for no decay."},
            {"ref_time", "ref_time_s", &Self::p_ref_time, 0.0, c_ParamBounds::Finite,
             "Time the rate applies at [s], on the same clock as the heating's time."},
        };
        return specs;
    }

    c_FixedRadiogenics() : c_FixedRadiogenics(c_ParamMap{}) {}
    explicit c_FixedRadiogenics(const c_ParamMap& params) : c_SpecModel("fixed") { this->p_initialize(params); }

    double get_ref_time() const noexcept override { return this->p_ref_time; }

    double calc_heating(double time, double mass) const override {
        if (this->p_average_half_life <= 0.0) {
            return mass * this->p_fixed_heat_production;
        }
        // The half-life guard reads the configured floor at call time.
        const double decay_constant = TidalPyConstants::d_LN_HALF / c_guard_denominator(this->p_average_half_life);
        return mass * this->p_fixed_heat_production * c_safe_exp(decay_constant * (time - this->p_ref_time));
    }

protected:
    double p_fixed_heat_production = 0.0;
    double p_average_half_life     = 0.0;
    double p_ref_time              = 0.0;
};

inline const c_ModelRegistry<c_RadiogenicsBase>& c_radiogenics_registry() {
    static const c_ModelRegistry<c_RadiogenicsBase> registry = {
        {{"off", "none"},
         BinaryClassID::OffRadiogenics,     &c_make_entry<c_RadiogenicsBase, c_OffRadiogenics>},
        {{"isotope", "isotopes"},
         BinaryClassID::IsotopeRadiogenics, &c_make_entry<c_RadiogenicsBase, c_IsotopeRadiogenics>},
        {{"fixed", "constant"},
         BinaryClassID::FixedRadiogenics,   &c_make_entry<c_RadiogenicsBase, c_FixedRadiogenics>},
    };
    return registry;
}

// The family's entry points, each one line over the generic registry functions.
inline std::unique_ptr<c_RadiogenicsBase> c_find_radiogenics(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_radiogenics_registry(), model_name, params);
}

// An isotope model whose isotopes carry labels, one per isotope (or none).
inline std::unique_ptr<c_RadiogenicsBase> c_make_isotope_radiogenics(
        const c_ParamMap& params, const std::vector<std::string>& isotope_names) {
    return std::make_unique<c_IsotopeRadiogenics>(params, isotope_names);
}

inline std::unique_ptr<c_RadiogenicsBase> c_radiogenics_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_radiogenics_registry(), in, force);
}

inline std::string c_radiogenics_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_radiogenics_registry(), model_name);
}

inline std::vector<std::string> c_radiogenics_model_names() {
    return c_model_names(c_radiogenics_registry());
}

}  // namespace tidalpy
