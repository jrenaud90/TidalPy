#pragma once
/* TidalPy's radiogenic heating models. All MKS.
 *
 * The isotope and fixed models share one reference time, so times can be measured from any fixed epoch
 * such as solar-system formation.
 *
 * References
 * ----------
 * - Hussmann and Spohn (2004); Turcotte and Schubert (2001): chondritic isotope data.
 * - Castillo-Rogez et al. (2007): long- and short-lived radiogenic isotopes.
 * - McDonough and Sun (1995): bulk silicate Earth elemental abundances.
 *
 * Binary payload: model name then the model parameters. Off and Fixed use the shared c_PhysicsBase
 * helpers; Isotope writes its variable-length isotope list itself.
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

#include "radiogenics_base_.hpp"
#include "../Utilities_x/math_x/numerics_.hpp"  // c_safe_exp
#include "constants_.hpp"

namespace tidalpy {

// Guard a half-life denominator that may approach zero; clamps to a signed numerical floor.
inline double rad_guard(double value) noexcept {
    const double floor_value = tidalpy_config_ptr->d_NUMERICAL_FLOOR;
    if (std::abs(value) < floor_value) {
        return (value < 0.0) ? -floor_value : floor_value;
    }
    return value;
}

// One radioactive isotope and its decay heating. Specific heating per unit layer mass is
//
//   q(t) = mass_frac * concentration * heat_production * exp(gamma * (t - t_ref))
//
// with gamma = ln(0.5) / half_life the (negative) decay constant. Both abundances are values at t_ref. A source
// that quotes the isotope's own concentration rather than its element's is entered with mass_frac = 1.
struct c_Isotope {
    std::string name;                       // isotope label (e.g. "U238")
    double heat_production = 0.0;      // specific heat production of the pure isotope [W/kg]
    double half_life = 0.0;      // half life [s]
    double mass_frac = 0.0;      // isotopic mass fraction within its element at t_ref [kg/kg]
    double concentration = 0.0;      // element concentration in the layer material at t_ref [kg/kg]

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

    // Decay constant gamma = ln(0.5) / half_life [1/s]; negative, larger in magnitude for short half lives.
    double decay_constant() const noexcept {
        return TidalPyConstants::d_LN_HALF / rad_guard(this->half_life);
    }

    // Specific heating [W/kg] at the given time. The guarded exponential gives NaN rather than inf when
    // the time is so far before ref_time that the back-extrapolation overflows.
    double specific_heating(double time, double ref_time) const noexcept {
        const double q_ref = this->mass_frac * this->concentration * this->heat_production;
        return q_ref * c_safe_exp(this->decay_constant() * (time - ref_time));
    }
};

// Combined construction parameters; each model reads only the fields it needs.
struct c_RadiogenicsConfig {
    // Isotope model.
    std::vector<c_Isotope> isotopes;

    // Fixed model.
    double fixed_heat_production = 0.0;    // lumped specific rate    [W/kg]
    double average_half_life        = 0.0;    // decay half life (<=0 => no decay) [s]

    // Shared: the epoch the rates and concentrations were measured at.
    double ref_time = 0.0;                    // reference time          [s]
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
// of the quoted 14010 to 14050 Myr.
inline std::vector<c_Isotope> c_castillo_rogez_2007_llri() {
    const double myr = TidalPyConstants::d_SECONDS_PER_MYR;
    return {
        c_Isotope("U238",  9.465e-5, 4468.0  * myr, 1.0, 26.2e-9),
        c_Isotope("U235",  5.687e-4, 703.81  * myr, 1.0, 8.2e-9),
        c_Isotope("Th232", 2.638e-5, 14030.0 * myr, 1.0, 53.8e-9),
        c_Isotope("K40",   2.917e-5, 1277.0  * myr, 1.0, 1104.0e-9),
    };
}

// Castillo-Rogez et al. (2007) short-lived isotopes at CAI formation: heat production, half lives, and initial
// isotopic ratios from their Table 5. Each element concentration is the Table 3 isotope concentration divided by
// that ratio, as the paper builds Table 3 (26Al: 5e-5 of 1.2 wt% aluminum is 600 ppb), so the product reproduces
// Table 3 exactly. The ratio to the stable isotope stands in for the mass fraction within the element, as it does
// in the paper. Where Table 5 quotes a range, the Fe60 ratio is the 1e-6 of the paper's short-lived-isotope models
// and the Al26 half life and Fe60 heat production are the middle of the range.
inline std::vector<c_Isotope> c_castillo_rogez_2007_slri() {
    const double myr = TidalPyConstants::d_SECONDS_PER_MYR;
    return {
        c_Isotope("Al26", 0.146, 0.723 * myr, 5.0e-5, 1.2e-2),
        c_Isotope("Fe60", 0.071, 1.5   * myr, 1.0e-6, 0.225),
        c_Isotope("Mn53", 0.027, 3.7   * myr, 1.0e-5, 2.57e-3),
    };
}

// Radiogenics disabled.
inline double rad_heating_off(double /*time*/, double /*mass*/) noexcept {
    return 0.0;
}

// Sum each isotope's specific heating, then scale by the layer mass.
inline double rad_heating_isotope(
        double time,
        double mass,
        const std::vector<c_Isotope>& isotopes,
        double ref_time) noexcept {
    double specific_heating = 0.0;
    for (const c_Isotope& isotope : isotopes) {
        specific_heating += isotope.specific_heating(time, ref_time);
    }
    return specific_heating * mass;
}

// One lumped rate with optional exponential decay; average_half_life <= 0 disables the decay.
inline double rad_heating_fixed(
        double time,
        double mass,
        double fixed_heat_production,
        double average_half_life,
        double ref_time) noexcept {
    if (average_half_life <= 0.0) {
        return mass * fixed_heat_production;
    }
    const double gamma = TidalPyConstants::d_LN_HALF / rad_guard(average_half_life);
    return mass * fixed_heat_production * c_safe_exp(gamma * (time - ref_time));
}

inline std::string rad_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

inline c_IsotopeDataset c_get_isotope_dataset(const std::string& name) {
    const std::string key = rad_to_lower(name);
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

// Off and Fixed serialize through the shared c_PhysicsBase helpers. Isotope carries a variable-length
// isotope list, so it writes its own payload after the shared header and model name.

// Radiogenics disabled (alias "none").
class c_OffRadiogenics : public c_RadiogenicsBase {
public:
    c_OffRadiogenics() : c_RadiogenicsBase("off") {}
    explicit c_OffRadiogenics(const c_RadiogenicsConfig& /*cfg*/) : c_RadiogenicsBase("off") {}
    ~c_OffRadiogenics() override = default;

    double calc_heating(double time, double mass) const override {
        return rad_heating_off(time, mass);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::OffRadiogenics));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// Sum over individually decaying isotopes.
class c_IsotopeRadiogenics : public c_RadiogenicsBase {
public:
    c_IsotopeRadiogenics() : c_RadiogenicsBase("isotope") {}
    explicit c_IsotopeRadiogenics(const c_RadiogenicsConfig& cfg)
        : c_RadiogenicsBase("isotope"),
          p_isotopes(cfg.isotopes),
          p_ref_time(cfg.ref_time) {}
    ~c_IsotopeRadiogenics() override = default;

    const std::vector<c_Isotope>& get_isotopes() const noexcept { return this->p_isotopes; }
    double get_ref_time()          const noexcept override { return this->p_ref_time; }
    std::size_t get_num_isotopes() const noexcept { return this->p_isotopes.size(); }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RadiogenicsBase::append_config_entries(out);
        std::vector<double> heat_production, half_lives, mass_fracs, concentrations;
        std::vector<std::string> names;
        for (const c_Isotope& isotope : this->p_isotopes) {
            heat_production.push_back(isotope.heat_production);
            half_lives.push_back(isotope.half_life);
            mass_fracs.push_back(isotope.mass_frac);
            concentrations.push_back(isotope.concentration);
            names.push_back(isotope.name);
        }
        out.push_back(c_config_doubles("heat_production_w_kg", heat_production));
        out.push_back(c_config_doubles("half_lives_s", half_lives));
        out.push_back(c_config_doubles("mass_fracs", mass_fracs));
        out.push_back(c_config_doubles("concentrations", concentrations));
        out.push_back(c_config_strings("isotope_names", names));
        out.push_back(c_config_double("ref_time_s", this->p_ref_time));
    }

    double calc_heating(double time, double mass) const override {
        return rad_heating_isotope(time, mass, this->p_isotopes, this->p_ref_time);
    }

    void write_binary(std::ostream& out) const override {
        const auto n = static_cast<uint64_t>(this->p_isotopes.size());
        uint64_t payload =
            binary_string_bytes(this->p_model_name)
            + sizeof(double)            // ref_time
            + sizeof(uint64_t);         // isotope count
        for (const c_Isotope& iso : this->p_isotopes) {
            payload += binary_string_bytes(iso.name) + 4 * sizeof(double);
        }
        write_binary_header(out, static_cast<uint32_t>(BinaryClassID::IsotopeRadiogenics), payload);
        write_binary_string(out, this->p_model_name);
        out.write(reinterpret_cast<const char*>(&this->p_ref_time), sizeof(double));
        out.write(reinterpret_cast<const char*>(&n), sizeof(uint64_t));
        for (const c_Isotope& iso : this->p_isotopes) {
            write_binary_string(out, iso.name);
            out.write(reinterpret_cast<const char*>(&iso.heat_production), sizeof(double));
            out.write(reinterpret_cast<const char*>(&iso.half_life),          sizeof(double));
            out.write(reinterpret_cast<const char*>(&iso.mass_frac),            sizeof(double));
            out.write(reinterpret_cast<const char*>(&iso.concentration),        sizeof(double));
        }
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write isotope radiogenics binary data");
        }
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->p_model_name = read_binary_string(in);
        in.read(reinterpret_cast<char*>(&this->p_ref_time), sizeof(double));
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(uint64_t));
        if (!in) { throw std::runtime_error("TidalPy: failed to read isotope radiogenics binary data"); }
        check_binary_count(in, n, sizeof(uint32_t) + 4 * sizeof(double), "isotope");
        this->p_isotopes.clear();
        this->p_isotopes.reserve(n);
        for (uint64_t i = 0; i < n; ++i) {
            c_Isotope iso;
            iso.name = read_binary_string(in);
            in.read(reinterpret_cast<char*>(&iso.heat_production), sizeof(double));
            in.read(reinterpret_cast<char*>(&iso.half_life),          sizeof(double));
            in.read(reinterpret_cast<char*>(&iso.mass_frac),            sizeof(double));
            in.read(reinterpret_cast<char*>(&iso.concentration),        sizeof(double));
            this->p_isotopes.push_back(std::move(iso));
        }
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read isotope radiogenics binary data");
        }
    }

protected:
    std::vector<c_Isotope> p_isotopes;
    double p_ref_time = 0.0;
};

// One lumped rate with optional decay (alias "constant").
class c_FixedRadiogenics : public c_RadiogenicsBase {
public:
    c_FixedRadiogenics() : c_RadiogenicsBase("fixed") {}
    explicit c_FixedRadiogenics(const c_RadiogenicsConfig& cfg)
        : c_RadiogenicsBase("fixed"),
          p_fixed_heat_production(cfg.fixed_heat_production),
          p_average_half_life(cfg.average_half_life),
          p_ref_time(cfg.ref_time) {}
    ~c_FixedRadiogenics() override = default;

    double get_fixed_heat_production() const noexcept { return this->p_fixed_heat_production; }
    double get_average_half_life()     const noexcept { return this->p_average_half_life; }
    double get_ref_time()              const noexcept override { return this->p_ref_time; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RadiogenicsBase::append_config_entries(out);
        out.push_back(c_config_double("fixed_heat_production_w_kg", this->p_fixed_heat_production));
        out.push_back(c_config_double("average_half_life_s", this->p_average_half_life));
        out.push_back(c_config_double("ref_time_s", this->p_ref_time));
    }

    double calc_heating(double time, double mass) const override {
        return rad_heating_fixed(
                time,
                mass,
                this->p_fixed_heat_production,
                this->p_average_half_life,
                this->p_ref_time);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::FixedRadiogenics),
            {this->p_fixed_heat_production, this->p_average_half_life, this->p_ref_time});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 3);
        this->p_fixed_heat_production = params[0];
        this->p_average_half_life = params[1];
        this->p_ref_time          = params[2];
    }

protected:
    double p_fixed_heat_production = 0.0;
    double p_average_half_life = 0.0;
    double p_ref_time          = 0.0;
};

// One value per model, so c_find_radiogenics dispatches without string comparisons.
enum class c_RadiogenicsModel : uint8_t {
    Off     = 0,
    Isotope = 1,
    Fixed   = 2,
};

// Model names are matched case-insensitively.
inline c_RadiogenicsModel c_radiogenics_model_from_name(const std::string& model_name) {
    const std::string name = rad_to_lower(model_name);

    if (name == "off"     || name == "none")     { return c_RadiogenicsModel::Off; }
    if (name == "isotope" || name == "isotopes") { return c_RadiogenicsModel::Isotope; }
    if (name == "fixed"   || name == "constant") { return c_RadiogenicsModel::Fixed; }

    throw std::invalid_argument("TidalPy: unknown radiogenics model name '" + model_name + "'");
}

// The canonical C++ factory: layers, binary reconstruction, and the Cython wrapper all route here.
inline std::unique_ptr<c_RadiogenicsBase> c_find_radiogenics(
        c_RadiogenicsModel model, const c_RadiogenicsConfig& cfg) {
    switch (model) {
        case c_RadiogenicsModel::Off:     return std::make_unique<c_OffRadiogenics>(cfg);
        case c_RadiogenicsModel::Isotope: return std::make_unique<c_IsotopeRadiogenics>(cfg);
        case c_RadiogenicsModel::Fixed:   return std::make_unique<c_FixedRadiogenics>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_RadiogenicsModel enum value");
}

inline std::unique_ptr<c_RadiogenicsBase> c_find_radiogenics(
        const std::string& model_name, const c_RadiogenicsConfig& cfg) {
    return c_find_radiogenics(c_radiogenics_model_from_name(model_name), cfg);
}

// The class id is peeked without consuming the header so the default-constructed model restores itself.
// Used by the layer recursive deserialization in structures_x/layers.
inline std::unique_ptr<c_RadiogenicsBase> c_radiogenics_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_RadiogenicsBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::OffRadiogenics:     model = std::make_unique<c_OffRadiogenics>();     break;
        case BinaryClassID::IsotopeRadiogenics: model = std::make_unique<c_IsotopeRadiogenics>(); break;
        case BinaryClassID::FixedRadiogenics:   model = std::make_unique<c_FixedRadiogenics>();   break;
        default:
            throw std::runtime_error("TidalPy: unknown radiogenics class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
