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

#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "radiogenics_base_.hpp"
#include "../Utilities_x/math_x/numerics_.hpp"  // c_safe_exp
#include "constants_.hpp"
#include "model_names_.hpp"

namespace tidalpy {

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
        return TidalPyConstants::d_LN_HALF / c_guard_denominator(this->half_life);
    }

    // Specific heating at t_ref, q(t_ref) [W/kg].
    double reference_heating() const noexcept {
        return this->mass_frac * this->concentration * this->heat_production;
    }

    // Specific heating [W/kg] at the given time. The guarded exponential gives NaN rather than inf when
    // the time is so far before ref_time that the back-extrapolation overflows.
    double specific_heating(double time, double ref_time) const noexcept {
        return this->reference_heating() * c_safe_exp(this->decay_constant() * (time - ref_time));
    }
};

// The time-independent factors of an isotope's q(t), formed once when a model's isotope list is set rather than on
// every evaluation. The decay constant is formed per call, because its half-life guard reads the configured floor at
// call time. specific_heating repeats c_Isotope::specific_heating's operations, so the two agree exactly.
struct c_IsotopeDecayTerms {
    double reference_heating = 0.0;  // q(t_ref) [W/kg]
    double half_life         = 0.0;  // [s]

    c_IsotopeDecayTerms() = default;
    explicit c_IsotopeDecayTerms(const c_Isotope& isotope) noexcept :
        reference_heating(isotope.reference_heating()),
        half_life(isotope.half_life) {}

    // Decay constant gamma = ln(0.5) / half_life [1/s], as c_Isotope::decay_constant.
    double decay_constant() const noexcept {
        return TidalPyConstants::d_LN_HALF / c_guard_denominator(this->half_life);
    }

    // Specific heating [W/kg] at elapsed_time = t - t_ref [s].
    double specific_heating(double elapsed_time) const noexcept {
        return this->reference_heating * c_safe_exp(this->decay_constant() * elapsed_time);
    }
};

// Combined construction parameters; each model reads only the fields it needs. Its defaults are the models' defaults.
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

// Radiogenics disabled.
inline double rad_heating_off(double /*time*/, double /*mass*/) noexcept {
    return 0.0;
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

// Off and Fixed serialize through the shared c_PhysicsBase helpers. Isotope carries a variable-length
// isotope list, so it writes its own payload after the shared header and model name.

// Radiogenics disabled (alias "none").
class c_OffRadiogenics final : public c_RadiogenicsBase {
public:
    c_OffRadiogenics() : c_OffRadiogenics(c_RadiogenicsConfig{}) {}
    explicit c_OffRadiogenics(const c_RadiogenicsConfig& /*cfg*/) : c_RadiogenicsBase("off") {}
    ~c_OffRadiogenics() override = default;

    double calc_heating(double time, double mass) const override {
        return rad_heating_off(time, mass);
    }

    uint32_t get_binary_class_id() const override {
        return static_cast<uint32_t>(BinaryClassID::OffRadiogenics);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id());
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// Sum over individually decaying isotopes.
class c_IsotopeRadiogenics final : public c_RadiogenicsBase {
public:
    c_IsotopeRadiogenics() : c_IsotopeRadiogenics(c_RadiogenicsConfig{}) {}
    explicit c_IsotopeRadiogenics(const c_RadiogenicsConfig& cfg)
        : c_RadiogenicsBase("isotope"),
          p_isotopes(cfg.isotopes),
          p_ref_time(cfg.ref_time) {
        this->p_cache_decay_terms();
    }
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

    // Sum each isotope's specific heating, then scale by the layer mass.
    double calc_heating(double time, double mass) const override {
        const double elapsed_time = time - this->p_ref_time;
        double specific_heating = 0.0;
        for (const c_IsotopeDecayTerms& decay_terms : this->p_decay_terms) {
            specific_heating += decay_terms.specific_heating(elapsed_time);
        }
        return specific_heating * mass;
    }

    // A sweep forms each isotope's decay constant once rather than once per point. Each point sums the same terms in
    // the same order as calc_heating, then scales by its mass, so the two agree exactly.
    void calc_heating_vectorize(
            const std::vector<double>& time,
            const std::vector<double>& mass,
            std::vector<double>& out_heating) const override {
        const std::size_t num_points = c_broadcast_length({time.size(), mass.size()}, "calc_heating_vectorize");
        const std::size_t time_stride = c_broadcast_stride(time.size());
        const std::size_t mass_stride = c_broadcast_stride(mass.size());
        out_heating.assign(num_points, 0.0);
        for (const c_IsotopeDecayTerms& decay_terms : this->p_decay_terms) {
            const double decay_constant = decay_terms.decay_constant();
            for (std::size_t i = 0; i < num_points; ++i) {
                const double elapsed_time = time[i * time_stride] - this->p_ref_time;
                out_heating[i] += decay_terms.reference_heating * c_safe_exp(decay_constant * elapsed_time);
            }
        }
        for (std::size_t i = 0; i < num_points; ++i) {
            out_heating[i] *= mass[i * mass_stride];
        }
    }

    uint32_t get_binary_class_id() const override {
        return static_cast<uint32_t>(BinaryClassID::IsotopeRadiogenics);
    }

    void write_binary(std::ostream& out) const override {
        const auto n = static_cast<uint64_t>(this->p_isotopes.size());
        write_binary_header(out, this->get_binary_class_id(), p_payload_bytes(this->p_model_name, this->p_isotopes));
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

    // The header's payload size must be exactly what the isotope list read from it occupies; any other size means the
    // record was written with a different layout (or is corrupt), so it raises even with force, which relaxes only
    // the schema-version check. The model changes only once the whole record is read and checked.
    void read_binary(std::istream& in, bool force = false) override {
        const c_BinaryHeader header = c_read_binary_record_header(in, force);
        std::string model_name = read_binary_string(in);
        double ref_time = 0.0;
        in.read(reinterpret_cast<char*>(&ref_time), sizeof(double));
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(uint64_t));
        if (!in) { throw std::runtime_error("TidalPy: failed to read isotope radiogenics binary data"); }
        check_binary_count(in, n, sizeof(uint32_t) + 4 * sizeof(double), "isotope");
        std::vector<c_Isotope> isotopes;
        isotopes.reserve(n);
        for (uint64_t i = 0; i < n; ++i) {
            c_Isotope iso;
            iso.name = read_binary_string(in);
            in.read(reinterpret_cast<char*>(&iso.heat_production), sizeof(double));
            in.read(reinterpret_cast<char*>(&iso.half_life),          sizeof(double));
            in.read(reinterpret_cast<char*>(&iso.mass_frac),            sizeof(double));
            in.read(reinterpret_cast<char*>(&iso.concentration),        sizeof(double));
            isotopes.push_back(std::move(iso));
        }
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read isotope radiogenics binary data");
        }
        const uint64_t expected_payload = p_payload_bytes(model_name, isotopes);
        if (header.payload_size != expected_payload) {
            throw std::runtime_error(
                "TidalPy: corrupt binary data: the isotope radiogenics record holds "
                + std::to_string(header.payload_size) + " payload bytes, but its " + std::to_string(n)
                + " isotopes occupy " + std::to_string(expected_payload)
                + ", so it was written with a different layout or is corrupt");
        }
        this->p_model_name = std::move(model_name);
        this->p_ref_time   = ref_time;
        this->p_isotopes   = std::move(isotopes);
        this->p_cache_decay_terms();
    }

protected:
    // Payload bytes of a record: the model name, the reference time, the isotope count, then each isotope's name
    // and four doubles.
    static uint64_t p_payload_bytes(const std::string& model_name, const std::vector<c_Isotope>& isotopes) {
        uint64_t payload = binary_string_bytes(model_name) + sizeof(double) + sizeof(uint64_t);
        for (const c_Isotope& iso : isotopes) {
            payload += binary_string_bytes(iso.name) + 4 * sizeof(double);
        }
        return payload;
    }

    void p_cache_decay_terms() {
        this->p_decay_terms.clear();
        this->p_decay_terms.reserve(this->p_isotopes.size());
        for (const c_Isotope& isotope : this->p_isotopes) {
            this->p_decay_terms.emplace_back(isotope);
        }
    }

    std::vector<c_Isotope> p_isotopes;
    double p_ref_time;
    // One entry per isotope, in p_isotopes order, rebuilt whenever p_isotopes is set.
    std::vector<c_IsotopeDecayTerms> p_decay_terms;
};

// One lumped rate with optional decay (alias "constant"); average_half_life <= 0 disables the decay.
class c_FixedRadiogenics final : public c_RadiogenicsBase {
public:
    c_FixedRadiogenics() : c_FixedRadiogenics(c_RadiogenicsConfig{}) {}
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
        if (this->p_average_half_life <= 0.0) {
            return mass * this->p_fixed_heat_production;
        }
        // The half-life guard reads the configured floor at call time.
        const double decay_constant = TidalPyConstants::d_LN_HALF / c_guard_denominator(this->p_average_half_life);
        return mass * this->p_fixed_heat_production * c_safe_exp(decay_constant * (time - this->p_ref_time));
    }

    uint32_t get_binary_class_id() const override {
        return static_cast<uint32_t>(BinaryClassID::FixedRadiogenics);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, this->get_binary_class_id(),
            {this->p_fixed_heat_production, this->p_average_half_life, this->p_ref_time});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 3);
        this->p_fixed_heat_production = params[0];
        this->p_average_half_life = params[1];
        this->p_ref_time          = params[2];
    }

protected:
    double p_fixed_heat_production;
    double p_average_half_life;
    double p_ref_time;
};

// One value per model, so c_find_radiogenics dispatches without string comparisons.
enum class c_RadiogenicsModel : uint8_t {
    Off     = 0,
    Isotope = 1,
    Fixed   = 2,
};

// Model names are matched case-insensitively.
inline c_RadiogenicsModel c_radiogenics_model_from_name(const std::string& model_name) {
    const std::string name = c_to_lower(model_name);

    if (name == "off"     || name == "none")     { return c_RadiogenicsModel::Off; }
    if (name == "isotope" || name == "isotopes") { return c_RadiogenicsModel::Isotope; }
    if (name == "fixed"   || name == "constant") { return c_RadiogenicsModel::Fixed; }

    throw std::invalid_argument("TidalPy: unknown radiogenics model name '" + model_name + "'");
}

// Builds a model from its enum value and parameters; the Cython wrappers construct through it. A saved record is
// restored by c_radiogenics_from_binary instead.
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
    const c_BinaryHeader header = c_peek_binary_header(in);

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
