#pragma once
/*
 * radiogenics_.hpp - Implements TidalPy's radiogenic heating models.
 *
 * Inherits c_RadiogenicsBase (radiogenics_base_.hpp), which itself inherits
 * c_PhysicsBase. Each model implements calc_heating(time, mass), returning the
 * total radiogenic heating [W] produced by the given mass at the given time.
 *
 * Models (with config aliases handled by the factory):
 *   c_OffRadiogenics     (alias "none")     — radiogenics disabled, heating == 0
 *   c_IsotopeRadiogenics                    — sum of individual decaying isotopes
 *   c_FixedRadiogenics   (alias "constant") — single lumped rate with optional decay
 *
 * A single radioactive isotope is described by the lightweight c_Isotope struct
 * (heat production, half life, isotopic mass fraction, element concentration).
 * The isotope model carries a std::vector<c_Isotope>.
 *
 * All quantities are MKS: time and half-lives in seconds [s], mass in [kg],
 * heat production rates in [W/kg], and heating in [W]. The isotope and fixed
 * models share a single reference time so that times can be expressed relative
 * to any fixed epoch (e.g. solar-system formation).
 *
 * References
 * ----------
 * - Hussmann and Spohn (2004); Turcotte and Schubert (2001) — chondritic isotope data.
 * - Castillo-Rogez et al. (2007) — long- and short-lived radiogenic isotopes.
 * - McDonough and Sun (1995) — bulk silicate Earth elemental abundances.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::<Model> (501-503)
 *   payload: model_name length (uint32_t) | model_name bytes | model params
 *   Off / Fixed write scalar doubles via the shared c_PhysicsBase helpers.
 *   Isotope writes its variable-length isotope list directly (see below).
 *   The layer observer pointer (p_layer_ptr) is NOT serialized.
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

namespace tidalpy {

// -------------------------------------------------------------------------------
// Module-level constants.
// -------------------------------------------------------------------------------

// Natural log of one half. std::log is not constexpr under C++20, so the value
// is given as a literal (ln(0.5) = -0.6931471805599453). The decay constant for
// a half life t_half is gamma = d_LN_HALF / t_half.
inline constexpr double d_LN_HALF = -0.6931471805599453094172321214581765680755001344;

// Seconds in one mega-year (Julian year, 365.25 days). Built-in isotope datasets
// quote half lives and reference times in Myr; they are converted to seconds at
// construction so the C++ API stays MKS.
inline constexpr double d_SECONDS_PER_MYR = 1.0e6 * 365.25 * 24.0 * 3600.0;

// Numerical floor used to guard half-life denominators that may approach zero.
inline constexpr double d_RADIOGENICS_FLOOR = 1.0e-100;

// Replace a magnitude smaller than the floor with a signed floor value.
inline double rad_guard(double value) noexcept {
    if (std::abs(value) < d_RADIOGENICS_FLOOR) {
        return (value < 0.0) ? -d_RADIOGENICS_FLOOR : d_RADIOGENICS_FLOOR;
    }
    return value;
}

// =====================================================================================================================
// c_Isotope — a single radioactive isotope and its decay heating
//
// A lightweight value type (no base class, no virtuals) describing one isotope's
// contribution to radiogenic heating. The specific heating per unit layer mass is
//
//   q(t) = mass_frac * concentration * heat_production * exp(gamma * (t - t_ref))
//
// where gamma = ln(0.5) / half_life is the (negative) decay constant.
// =====================================================================================================================
struct c_Isotope {
    std::string name;                       // isotope label (e.g. "U238")
    double heat_production_w_kg = 0.0;      // specific heat production of the pure isotope [W/kg]
    double half_life_s          = 0.0;      // half life [s]
    double mass_frac            = 0.0;      // isotopic mass fraction within its element [kg/kg]
    double concentration        = 0.0;      // element concentration in the layer material [kg/kg]

    c_Isotope() = default;
    c_Isotope(std::string isotope_name,
              double hpr_w_kg,
              double half_life_seconds,
              double isotopic_mass_frac,
              double element_concentration)
        : name(std::move(isotope_name)),
          heat_production_w_kg(hpr_w_kg),
          half_life_s(half_life_seconds),
          mass_frac(isotopic_mass_frac),
          concentration(element_concentration) {}

    // Decay constant gamma = ln(0.5) / half_life [1/s] (negative; magnitude grows
    // as the half life shortens).
    double decay_constant() const noexcept {
        return d_LN_HALF / rad_guard(this->half_life_s);
    }

    // Specific radiogenic heating per unit layer mass [W/kg] at the given time.
    // The guarded exponential returns NaN (rather than inf) if the requested time is so far
    // before the reference time that the back-extrapolated heating overflows.
    double specific_heating(double time_s, double ref_time_s) const noexcept {
        const double q_ref = this->mass_frac * this->concentration * this->heat_production_w_kg;
        return q_ref * c_safe_exp(this->decay_constant() * (time_s - ref_time_s));
    }
};

// -------------------------------------------------------------------------------
// c_RadiogenicsConfig — combined construction parameters for all models.
// Each model reads only the fields it needs.
// -------------------------------------------------------------------------------
struct c_RadiogenicsConfig {
    // Isotope model — one c_Isotope per radioactive isotope.
    std::vector<c_Isotope> isotopes;

    // Fixed model.
    double fixed_heat_production_w_kg = 0.0;    // lumped specific rate    [W/kg]
    double average_half_life_s        = 0.0;    // decay half life (<=0 => no decay) [s]

    // Shared reference time at which the rates/concentrations were measured.
    double ref_time_s = 0.0;                    // reference time          [s]
};

// =====================================================================================================================
// Built-in isotope datasets
//
// Convenience catalogs of well-characterized radiogenic isotope sets from the
// literature, so a caller can build a realistic c_IsotopeRadiogenics without
// hand-entering abundances. Half lives and reference times are quoted in Myr in
// the source literature and converted to seconds here (MKS).
// =====================================================================================================================
struct c_IsotopeDataset {
    std::vector<c_Isotope> isotopes;
    double ref_time_s = 0.0;
};

// Names of the available built-in datasets (see c_get_isotope_dataset).
inline std::vector<std::string> c_isotope_dataset_names() {
    return {"modern_day_chondritic", "llri_and_slri", "bulk_silicate_earth"};
}

// Build a named built-in isotope dataset. Throws std::invalid_argument on an
// unknown name. Recognized names (case-insensitive):
//
//   "modern_day_chondritic"
//       Present-day chondritic abundances of the four long-lived heat producers
//       (U238, U235, Th232, K40). Applicable to rocky/icy bodies of broadly
//       chondritic composition evaluated near the present epoch.
//       Source: Hussmann and Spohn (2004); Turcotte and Schubert (2001).
//
//   "llri_and_slri"
//       Long-lived (U238, U235, Th232, K40) plus short-lived (Mn53, Fe60, Al26)
//       radiogenic isotopes. Applicable to early-solar-system thermal evolution
//       (the first few Myr) where short-lived isotopes dominate the heat budget.
//       Source: Castillo-Rogez et al. (2007).
//
//   "bulk_silicate_earth"
//       Present-day bulk silicate Earth: the four long-lived heat producers with
//       BSE elemental concentrations (U = 20.3 ppb, Th = 79.5 ppb, K = 240 ppm).
//       Applicable to Earth-like silicate mantles at the present epoch.
//       Source: McDonough and Sun (1995) (concentrations); heat production rates
//       and half lives from Turcotte and Schubert (2002).
//
// All reference times are 4600 Myr (i.e. concentrations are quoted at the present
// epoch, 4.6 Gyr after solar-system formation).
inline c_IsotopeDataset c_get_isotope_dataset(const std::string& name);

// =====================================================================================================================
// Radiogenic heating functions [W]
//
// Each mirrors the validated legacy implementation in
// TidalPy/radiogenics/radiogenic_models.py.
// =====================================================================================================================

// Off: radiogenics disabled, heating == 0.
inline double rad_heating_off(double /*time_s*/, double /*mass_kg*/) noexcept {
    return 0.0;
}

// Isotope: sum each isotope's specific heating, then scale by the layer mass.
inline double rad_heating_isotope(
        double time_s,
        double mass_kg,
        const std::vector<c_Isotope>& isotopes,
        double ref_time_s) noexcept {
    double specific_heating = 0.0;
    for (const c_Isotope& isotope : isotopes) {
        specific_heating += isotope.specific_heating(time_s, ref_time_s);
    }
    return specific_heating * mass_kg;
}

// Fixed: single lumped rate with optional exponential decay.
// average_half_life_s <= 0 disables decay (constant heating rate).
inline double rad_heating_fixed(
        double time_s,
        double mass_kg,
        double fixed_heat_production_w_kg,
        double average_half_life_s,
        double ref_time_s) noexcept {
    if (average_half_life_s <= 0.0) {
        return mass_kg * fixed_heat_production_w_kg;
    }
    const double gamma = d_LN_HALF / rad_guard(average_half_life_s);
    return mass_kg * fixed_heat_production_w_kg * std::exp(gamma * (time_s - ref_time_s));
}

// -------------------------------------------------------------------------------
// Lower-case a model name for case-insensitive factory lookup.
// -------------------------------------------------------------------------------
inline std::string rad_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// -------------------------------------------------------------------------------
// c_get_isotope_dataset implementation (declared above).
// -------------------------------------------------------------------------------
inline c_IsotopeDataset c_get_isotope_dataset(const std::string& name) {
    const std::string key = rad_to_lower(name);
    const double myr = d_SECONDS_PER_MYR;
    c_IsotopeDataset dataset;
    dataset.ref_time_s = 4600.0 * myr;

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
    if (key == "llri_and_slri") {
        // Castillo-Rogez et al. (2007).
        dataset.isotopes = {
            c_Isotope("U238",  9.465e-5, 4468.0   * myr, 0.9928,   0.026e-6),
            c_Isotope("U235",  5.687e-4, 703.81   * myr, 0.0071,   0.0082e-6),
            c_Isotope("Th232", 2.638e-5, 14025.0  * myr, 1.0,      0.0538e-6),
            c_Isotope("K40",   2.917e-5, 1277.0   * myr, 1.176e-4, 1.104e-6),
            c_Isotope("Mn53",  0.027,    3.7      * myr, 2.0e-5,   0.0257e-6),
            c_Isotope("Fe60",  0.07,     1.5      * myr, 1.0e-6,   0.1e-6),
            c_Isotope("Al26",  0.146,    0.72     * myr, 5.0e-5,   0.6e-6),
        };
        return dataset;
    }
    if (key == "bulk_silicate_earth") {
        // McDonough and Sun (1995) concentrations; Turcotte and Schubert (2002)
        // heat production rates and half lives. BSE: U 20.3 ppb, Th 79.5 ppb, K 240 ppm.
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

// =====================================================================================================================
// Radiogenics models
//
// Binary serialization uses the shared c_PhysicsBase helpers
// (write_physics_binary / read_physics_binary) for the scalar-only models (Off,
// Fixed). The Isotope model has a variable-length isotope list, so it writes its
// own payload after the shared header + model name (a justified exception to the
// "single helper call" rule for fixed-parameter models).
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_OffRadiogenics — radiogenics disabled (heating == 0).
// -------------------------------------------------------------------------------
class c_OffRadiogenics : public c_RadiogenicsBase {
public:
    c_OffRadiogenics() : c_RadiogenicsBase("off") {}
    explicit c_OffRadiogenics(const c_RadiogenicsConfig& /*cfg*/) : c_RadiogenicsBase("off") {}
    ~c_OffRadiogenics() override = default;

    double calc_heating(double time_s, double mass_kg) const override {
        return rad_heating_off(time_s, mass_kg);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::OffRadiogenics));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// -------------------------------------------------------------------------------
// c_IsotopeRadiogenics — sum of individual decaying isotopes.
// -------------------------------------------------------------------------------
class c_IsotopeRadiogenics : public c_RadiogenicsBase {
public:
    c_IsotopeRadiogenics() : c_RadiogenicsBase("isotope") {}
    explicit c_IsotopeRadiogenics(const c_RadiogenicsConfig& cfg)
        : c_RadiogenicsBase("isotope"),
          p_isotopes(cfg.isotopes),
          p_ref_time_s(cfg.ref_time_s) {}
    ~c_IsotopeRadiogenics() override = default;

    const std::vector<c_Isotope>& get_isotopes() const noexcept { return this->p_isotopes; }
    double get_ref_time()          const noexcept { return this->p_ref_time_s; }
    std::size_t get_num_isotopes() const noexcept { return this->p_isotopes.size(); }

    double calc_heating(double time_s, double mass_kg) const override {
        return rad_heating_isotope(time_s, mass_kg, this->p_isotopes, this->p_ref_time_s);
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
        out.write(reinterpret_cast<const char*>(&this->p_ref_time_s), sizeof(double));
        out.write(reinterpret_cast<const char*>(&n), sizeof(uint64_t));
        for (const c_Isotope& iso : this->p_isotopes) {
            write_binary_string(out, iso.name);
            out.write(reinterpret_cast<const char*>(&iso.heat_production_w_kg), sizeof(double));
            out.write(reinterpret_cast<const char*>(&iso.half_life_s),          sizeof(double));
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
        in.read(reinterpret_cast<char*>(&this->p_ref_time_s), sizeof(double));
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(uint64_t));
        this->p_isotopes.clear();
        this->p_isotopes.reserve(n);
        for (uint64_t i = 0; i < n; ++i) {
            c_Isotope iso;
            iso.name = read_binary_string(in);
            in.read(reinterpret_cast<char*>(&iso.heat_production_w_kg), sizeof(double));
            in.read(reinterpret_cast<char*>(&iso.half_life_s),          sizeof(double));
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
    double p_ref_time_s = 0.0;
};

// -------------------------------------------------------------------------------
// c_FixedRadiogenics — single lumped rate with optional decay (alias "constant").
// -------------------------------------------------------------------------------
class c_FixedRadiogenics : public c_RadiogenicsBase {
public:
    c_FixedRadiogenics() : c_RadiogenicsBase("fixed") {}
    explicit c_FixedRadiogenics(const c_RadiogenicsConfig& cfg)
        : c_RadiogenicsBase("fixed"),
          p_fixed_heat_production_w_kg(cfg.fixed_heat_production_w_kg),
          p_average_half_life_s(cfg.average_half_life_s),
          p_ref_time_s(cfg.ref_time_s) {}
    ~c_FixedRadiogenics() override = default;

    double get_fixed_heat_production() const noexcept { return this->p_fixed_heat_production_w_kg; }
    double get_average_half_life()     const noexcept { return this->p_average_half_life_s; }
    double get_ref_time()              const noexcept { return this->p_ref_time_s; }

    double calc_heating(double time_s, double mass_kg) const override {
        return rad_heating_fixed(
            time_s, mass_kg,
            this->p_fixed_heat_production_w_kg,
            this->p_average_half_life_s,
            this->p_ref_time_s);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::FixedRadiogenics),
            {this->p_fixed_heat_production_w_kg, this->p_average_half_life_s, this->p_ref_time_s});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 3);
        this->p_fixed_heat_production_w_kg = params[0];
        this->p_average_half_life_s        = params[1];
        this->p_ref_time_s                 = params[2];
    }

protected:
    double p_fixed_heat_production_w_kg = 0.0;
    double p_average_half_life_s        = 0.0;
    double p_ref_time_s                 = 0.0;
};

// =====================================================================================================================
// Factory
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_RadiogenicsModel — one named value per radiogenics model.
// Used by c_find_radiogenics to dispatch to the correct class without string
// comparisons. The Cython layer maps Python strings (and aliases) to these
// values via c_radiogenics_model_from_name.
// -------------------------------------------------------------------------------
enum class c_RadiogenicsModel : uint8_t {
    Off     = 0,
    Isotope = 1,
    Fixed   = 2,
};

// -------------------------------------------------------------------------------
// c_radiogenics_model_from_name — map a (case-insensitive) model name or alias to
// a c_RadiogenicsModel enum value.
//
// Recognized names and aliases:
//   "off" / "none"
//   "isotope" / "isotopes"
//   "fixed" / "constant"
//
// Throws std::invalid_argument if the model name is unknown.
// -------------------------------------------------------------------------------
inline c_RadiogenicsModel c_radiogenics_model_from_name(const std::string& model_name) {
    const std::string name = rad_to_lower(model_name);

    if (name == "off"     || name == "none")     { return c_RadiogenicsModel::Off; }
    if (name == "isotope" || name == "isotopes") { return c_RadiogenicsModel::Isotope; }
    if (name == "fixed"   || name == "constant") { return c_RadiogenicsModel::Fixed; }

    throw std::invalid_argument("TidalPy: unknown radiogenics model name '" + model_name + "'");
}

// -------------------------------------------------------------------------------
// c_find_radiogenics — build the radiogenics model named by a c_RadiogenicsModel.
//
// Returns a unique_ptr to a newly heap-allocated concrete model constructed from
// the supplied config. This is the canonical C++ factory; C++ consumers (layers
// attaching radiogenics, binary reconstruction, the Cython wrapper) all route
// through it. Throws std::invalid_argument for an unrecognised enum value.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_RadiogenicsBase> c_find_radiogenics(
        c_RadiogenicsModel model, const c_RadiogenicsConfig& cfg) {
    switch (model) {
        case c_RadiogenicsModel::Off:     return std::make_unique<c_OffRadiogenics>(cfg);
        case c_RadiogenicsModel::Isotope: return std::make_unique<c_IsotopeRadiogenics>(cfg);
        case c_RadiogenicsModel::Fixed:   return std::make_unique<c_FixedRadiogenics>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_RadiogenicsModel enum value");
}

// -------------------------------------------------------------------------------
// c_find_radiogenics (name overload) — maps a name/alias to the enum and builds
// the model. Throws std::invalid_argument on unknown names.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_RadiogenicsBase> c_find_radiogenics(
        const std::string& model_name, const c_RadiogenicsConfig& cfg) {
    return c_find_radiogenics(c_radiogenics_model_from_name(model_name), cfg);
}

// -------------------------------------------------------------------------------
// c_radiogenics_from_binary — reconstruct a radiogenics model from a binary stream.
//
// Peeks the upcoming record's BinaryClassID (without consuming the header),
// constructs the matching default-initialized concrete model, then delegates to
// its read_binary to restore the model name and parameters. Used by the layer
// recursive deserialization (see structures_x/layers). Throws std::runtime_error
// if the class id is not a known radiogenics model.
// -------------------------------------------------------------------------------
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
