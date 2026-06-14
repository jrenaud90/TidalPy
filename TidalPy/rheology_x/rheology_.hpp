#pragma once
/*
 * rheology_.hpp - Implements TidalPy's rheology models.
 *
 * Inherits c_RheologyBase (rheology_base_.hpp), which itself inherits
 * c_PhysicsBase.  Each model implements calc_complex_modulus(modulus, viscosity, frequency),
 * returning the complex (shear/bulk) modulus [Pa] directly.
 *
 * Models (with config aliases handled by the factory):
 *   c_Elastic   (alias "off")             — purely elastic, no dissipation
 *   c_Viscous   (alias "newton")          — purely viscous (Newtonian fluid)
 *   c_Voigt     (alias "voigt-kelvin")    — Voigt-Kelvin element
 *   c_Maxwell                             — standard Maxwell body
 *   c_Burgers                             — Maxwell + Voigt in series
 *   c_Andrade                             — Maxwell + Andrade transient term
 *   c_Sundberg  (alias "sundberg-cooper") — Andrade + Voigt
 *
 * All complex moduli are in [Pa]; the internal element
 * Inputs are reference (background) values at the layer mid-point in MKS units.
 *
 * References
 * ----------
 * - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000
 *   (Maxwell, Voigt-Kelvin, Burgers).
 * - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
 *   (complex compliances and Love numbers).
 * - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784
 *   (Andrade and Sundberg-Cooper).
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::<Model> (301-307)
 *   payload: model_name length (uint32_t) | model_name bytes | model params (doubles)
 *   The layer observer pointer (p_layer_ptr) is NOT serialized.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants_.hpp"
#include "rheology_base_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// Numerical floor used to guard denominators that may approach zero (e.g. at
// zero forcing frequency).
// -------------------------------------------------------------------------------
inline constexpr double d_RHEOLOGY_FLOOR = 1.0e-100;

// Replace a magnitude smaller than the floor with a signed floor value.
inline double rheo_guard(double value) noexcept {
    if (std::abs(value) < d_RHEOLOGY_FLOOR) {
        return (value < 0.0) ? -d_RHEOLOGY_FLOOR : d_RHEOLOGY_FLOOR;
    }
    return value;
}

// -------------------------------------------------------------------------------
// c_RheologyConfig — combined construction parameters for all rheology models.
// Each model reads only the fields it needs.
// -------------------------------------------------------------------------------
struct c_RheologyConfig {
    double alpha                = 0.3;     // Andrade exponent           [dimensionless]
    double zeta                 = 1.0;     // Andrade timescale ratio    [dimensionless]
    double voigt_modulus_frac   = 5.0;     // Voigt modulus fraction     [dimensionless]
    double voigt_viscosity_frac = 0.02;    // Voigt viscosity fraction   [dimensionless]
};

// =====================================================================================================================
// Internal element compliances [Pa^-1]
//
// The public interface of every rheology model is the complex (shear/bulk)
// MODULUS (see the rheo_modulus_* functions below).  However, the standard
// composite rheologies (Burgers, Andrade, Sundberg) combine their constituent
// viscoelastic elements *in series*, which means the element COMPLIANCES add
// and the resulting modulus is the reciprocal of that sum.  These element
// compliances are therefore required intermediate quantities; they are kept in
// an internal `detail` namespace and are never exposed to Python.
//
// Each takes the (shear or bulk) viscosity [Pa·s], the unrelaxed modulus [Pa]
// (from which the static compliance J = 1/modulus is derived), and the forcing
// frequency [rad/s].  The math mirrors the validated legacy implementations in
// TidalPy/rheology/complex_compliance/compliance_models.py.
//
// Assumptions
// -----------
// - Linear viscoelastic regime, single forcing frequency.
// - Andrade-family models assume a positive forcing frequency.
// =====================================================================================================================
namespace detail {

// Maxwell element compliance: J* = J - i / (viscosity * frequency).
inline c_ComplexCompliance element_compliance_maxwell(
        double modulus_pa,
        double viscosity_pas,
        double frequency_rad_s) noexcept {
    const double static_compliance = 1.0 / rheo_guard(modulus_pa);
    const double denom = rheo_guard(viscosity_pas * frequency_rad_s);
    return c_ComplexCompliance(static_compliance, -1.0 / denom);
}

// Voigt-Kelvin element compliance using fractional modulus/viscosity offsets.
// The Voigt arm's compliance is the layer compliance divided by the modulus
// fraction: J_voigt = (1 / modulus) / voigt_modulus_frac.
inline c_ComplexCompliance element_compliance_voigt(
        double modulus_pa,
        double viscosity_pas,
        double frequency_rad_s,
        double voigt_modulus_frac,
        double voigt_viscosity_frac) noexcept {
    const double static_compliance = 1.0 / rheo_guard(modulus_pa);
    const double voigt_compliance  = static_compliance / rheo_guard(voigt_modulus_frac);
    const double voigt_viscosity   = voigt_viscosity_frac * viscosity_pas;

    const double scaled  = voigt_compliance * voigt_viscosity * frequency_rad_s;
    const double denom   = scaled * scaled + 1.0;
    const double real_j  = voigt_compliance / denom;
    const double imag_j  = -(voigt_compliance * voigt_compliance) * voigt_viscosity
                           * frequency_rad_s / denom;
    return c_ComplexCompliance(real_j, imag_j);
}

// Andrade element compliance: Maxwell compliance plus a transient term ~ omega^{-alpha}.
inline c_ComplexCompliance element_compliance_andrade(
        double modulus_pa,
        double viscosity_pas,
        double frequency_rad_s,
        double alpha,
        double zeta) noexcept {
    const double static_compliance = 1.0 / rheo_guard(modulus_pa);
    const double andrade_term =
        rheo_guard(static_compliance * viscosity_pas * frequency_rad_s * zeta);

    const double const_term =
        static_compliance * std::pow(andrade_term, -alpha) * std::tgamma(1.0 + alpha);
    const double half_pi_alpha = alpha * TidalPyConstants::d_PI / 2.0;
    const c_ComplexCompliance andrade_transient(
        std::cos(half_pi_alpha) * const_term,
        -std::sin(half_pi_alpha) * const_term
    );

    return element_compliance_maxwell(modulus_pa, viscosity_pas, frequency_rad_s)
         + andrade_transient;
}

}  // namespace detail

// =====================================================================================================================
// Complex (shear/bulk) modulus functions [Pa]
//
// Model constitutive laws return the complex modulus mu* directly:
// real = storage (in-phase), imag = loss (out-of-phase, positive
// = energy loss). Simple models are evaluated analytically; series composites
// invert the sum of their element compliances (see the detail namespace above).
// =====================================================================================================================

// Elastic: mu* = modulus (real).  No dissipation, frequency-independent.
inline c_ComplexModulus rheo_modulus_elastic(
        double modulus_pa,
        double /*viscosity_pas*/,
        double /*frequency_rad_s*/) noexcept {
    return c_ComplexModulus(modulus_pa, 0.0);
}

// Viscous (Newton): mu* = i * viscosity * frequency (purely dissipative).
inline c_ComplexModulus rheo_modulus_viscous(
        double /*modulus_pa*/,
        double viscosity_pas,
        double frequency_rad_s) noexcept {
    return c_ComplexModulus(0.0, viscosity_pas * frequency_rad_s);
}

// Maxwell: mu* = 1 / J_maxwell.
inline c_ComplexModulus rheo_modulus_maxwell(
        double modulus_pa,
        double viscosity_pas,
        double frequency_rad_s) noexcept {
    return c_ComplexModulus(1.0, 0.0)
         / detail::element_compliance_maxwell(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s);
}

// Voigt-Kelvin: mu* = 1 / J_voigt.
inline c_ComplexModulus rheo_modulus_voigt(
        double modulus_pa,
        double viscosity_pas,
        double frequency_rad_s,
        double voigt_modulus_frac,
        double voigt_viscosity_frac) noexcept {
    return c_ComplexModulus(1.0, 0.0)
         / detail::element_compliance_voigt(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s,
            voigt_modulus_frac,
            voigt_viscosity_frac);
}

// Burgers: Maxwell and Voigt elements in series; mu* = 1 / (J_maxwell + J_voigt).
inline c_ComplexModulus rheo_modulus_burgers(
        double modulus_pa,
        double viscosity_pas,
        double frequency_rad_s,
        double voigt_modulus_frac,
        double voigt_viscosity_frac) noexcept {
    const c_ComplexCompliance total =
        detail::element_compliance_maxwell(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s)
      + detail::element_compliance_voigt(
        modulus_pa,
        viscosity_pas,
        frequency_rad_s,
        voigt_modulus_frac,
        voigt_viscosity_frac);
    return c_ComplexModulus(1.0, 0.0) / total;
}

// Andrade: mu* = 1 / J_andrade.
inline c_ComplexModulus rheo_modulus_andrade(
        double modulus_pa,
        double viscosity_pas,
        double frequency_rad_s,
        double alpha,
        double zeta) noexcept {
    return c_ComplexModulus(1.0, 0.0)
         / detail::element_compliance_andrade(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s,
            alpha,
            zeta);
}

// Sundberg-Cooper: Andrade and Voigt elements in series; mu* = 1 / (J_andrade + J_voigt).
inline c_ComplexModulus rheo_modulus_sundberg(
        double modulus_pa,
        double viscosity_pas,
        double frequency_rad_s,
        double alpha,
        double zeta,
        double voigt_modulus_frac,
        double voigt_viscosity_frac) noexcept {
    const c_ComplexCompliance total =
        detail::element_compliance_andrade(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s,
            alpha,
            zeta)
      + detail::element_compliance_voigt(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s,
            voigt_modulus_frac,
            voigt_viscosity_frac);
    return c_ComplexModulus(1.0, 0.0) / total;
}

// -------------------------------------------------------------------------------
// Lower-case a model name for case-insensitive factory lookup.
// -------------------------------------------------------------------------------
inline std::string rheo_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// =====================================================================================================================
// Rheology models
//
// Binary serialization uses the shared c_PhysicsBase helpers
// (write_physics_binary / read_physics_binary), so each model only supplies its
// BinaryClassID and its scalar params — the header, model name, and byte layout
// are handled in the base class (see Utilities_x/classes_x/physics_base_.hpp).
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_Elastic — purely elastic response (alias "off").
// -------------------------------------------------------------------------------
class c_Elastic : public c_RheologyBase {
public:
    c_Elastic() : c_RheologyBase("elastic") {}
    explicit c_Elastic(const c_RheologyConfig& /*cfg*/) : c_RheologyBase("elastic") {}
    ~c_Elastic() override = default;

    c_ComplexModulus calc_complex_modulus(
            double modulus_pa,
            double viscosity_pas,
            double frequency_rad_s) const override {
        return rheo_modulus_elastic(
            modulus_pa, 
            viscosity_pas,
            frequency_rad_s);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::Elastic));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// -------------------------------------------------------------------------------
// c_Viscous — purely viscous response (alias "newton").
// -------------------------------------------------------------------------------
class c_Viscous : public c_RheologyBase {
public:
    c_Viscous() : c_RheologyBase("viscous") {}
    explicit c_Viscous(const c_RheologyConfig& /*cfg*/) : c_RheologyBase("viscous") {}
    ~c_Viscous() override = default;

    c_ComplexModulus calc_complex_modulus(
            double modulus_pa,
            double viscosity_pas,
            double frequency_rad_s) const override {
        return rheo_modulus_viscous(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::Viscous));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// -------------------------------------------------------------------------------
// c_Maxwell — standard Maxwell body.
// -------------------------------------------------------------------------------
class c_Maxwell : public c_RheologyBase {
public:
    c_Maxwell() : c_RheologyBase("maxwell") {}
    explicit c_Maxwell(const c_RheologyConfig& /*cfg*/) : c_RheologyBase("maxwell") {}
    ~c_Maxwell() override = default;

    c_ComplexModulus calc_complex_modulus(
            double modulus_pa,
            double viscosity_pas,
            double frequency_rad_s) const override {
        return rheo_modulus_maxwell(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::Maxwell));
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// -------------------------------------------------------------------------------
// c_Voigt — Voigt-Kelvin element (alias "voigt-kelvin").
// -------------------------------------------------------------------------------
class c_Voigt : public c_RheologyBase {
public:
    c_Voigt() : c_RheologyBase("voigt") {}
    explicit c_Voigt(const c_RheologyConfig& cfg)
        : c_RheologyBase("voigt"),
          p_voigt_modulus_frac(cfg.voigt_modulus_frac),
          p_voigt_viscosity_frac(cfg.voigt_viscosity_frac) {}
    ~c_Voigt() override = default;

    double get_voigt_modulus_frac() const noexcept { return this->p_voigt_modulus_frac; }
    double get_voigt_viscosity_frac() const noexcept { return this->p_voigt_viscosity_frac; }

    c_ComplexModulus calc_complex_modulus(
            double modulus_pa,
            double viscosity_pas,
            double frequency_rad_s) const override {
        return rheo_modulus_voigt(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s,
            this->p_voigt_modulus_frac,
            this->p_voigt_viscosity_frac);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::Voigt),
                                   {this->p_voigt_modulus_frac, this->p_voigt_viscosity_frac});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 2);
        this->p_voigt_modulus_frac   = params[0];
        this->p_voigt_viscosity_frac = params[1];
    }

protected:
    double p_voigt_modulus_frac = 5.0;
    double p_voigt_viscosity_frac  = 0.02;
};

// -------------------------------------------------------------------------------
// c_Burgers — Maxwell + Voigt in series.
// -------------------------------------------------------------------------------
class c_Burgers : public c_RheologyBase {
public:
    c_Burgers() : c_RheologyBase("burgers") {}
    explicit c_Burgers(const c_RheologyConfig& cfg)
        : c_RheologyBase("burgers"),
          p_voigt_modulus_frac(cfg.voigt_modulus_frac),
          p_voigt_viscosity_frac(cfg.voigt_viscosity_frac) {}
    ~c_Burgers() override = default;

    double get_voigt_modulus_frac()   const noexcept { return this->p_voigt_modulus_frac; }
    double get_voigt_viscosity_frac() const noexcept { return this->p_voigt_viscosity_frac; }

    c_ComplexModulus calc_complex_modulus(
            double modulus_pa, 
            double viscosity_pas,
            double frequency_rad_s) const override {
        return rheo_modulus_burgers(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s,
            this->p_voigt_modulus_frac,
            this->p_voigt_viscosity_frac);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::Burgers),
                                   {this->p_voigt_modulus_frac, this->p_voigt_viscosity_frac});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 2);
        this->p_voigt_modulus_frac   = params[0];
        this->p_voigt_viscosity_frac = params[1];
    }

protected:
    double p_voigt_modulus_frac = 5.0;
    double p_voigt_viscosity_frac  = 0.02;
};

// -------------------------------------------------------------------------------
// c_Andrade — Maxwell + Andrade transient term.
// -------------------------------------------------------------------------------
class c_Andrade : public c_RheologyBase {
public:
    c_Andrade() : c_RheologyBase("andrade") {}
    explicit c_Andrade(const c_RheologyConfig& cfg)
        : c_RheologyBase("andrade"),
          p_alpha(cfg.alpha),
          p_zeta(cfg.zeta) {}
    ~c_Andrade() override = default;

    double get_alpha() const noexcept { return this->p_alpha; }
    double get_zeta()  const noexcept { return this->p_zeta; }

    c_ComplexModulus calc_complex_modulus(
            double modulus_pa,
            double viscosity_pas,
            double frequency_rad_s) const override {
        return rheo_modulus_andrade(
            modulus_pa,
            viscosity_pas,
            frequency_rad_s,
            this->p_alpha,
            this->p_zeta);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::Andrade),
                                   {this->p_alpha, this->p_zeta});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 2);
        this->p_alpha = params[0];
        this->p_zeta  = params[1];
    }

protected:
    double p_alpha = 0.3;
    double p_zeta  = 1.0;
};

// -------------------------------------------------------------------------------
// c_Sundberg — Andrade + Voigt (alias "sundberg-cooper").
// -------------------------------------------------------------------------------
class c_Sundberg : public c_RheologyBase {
public:
    c_Sundberg() : c_RheologyBase("sundberg") {}
    explicit c_Sundberg(const c_RheologyConfig& cfg)
        : c_RheologyBase("sundberg"),
          p_alpha(cfg.alpha),
          p_zeta(cfg.zeta),
          p_voigt_modulus_frac(cfg.voigt_modulus_frac),
          p_voigt_viscosity_frac(cfg.voigt_viscosity_frac) {}
    ~c_Sundberg() override = default;

    double get_alpha()                 const noexcept { return this->p_alpha; }
    double get_zeta()                  const noexcept { return this->p_zeta; }
    double get_voigt_modulus_frac() const noexcept { return this->p_voigt_modulus_frac; }
    double get_voigt_viscosity_frac()  const noexcept { return this->p_voigt_viscosity_frac; }

    c_ComplexModulus calc_complex_modulus(
            double modulus_pa,
            double viscosity_pas,
            double frequency_rad_s) const override {
        return rheo_modulus_sundberg(
            modulus_pa, 
            viscosity_pas,
            frequency_rad_s,
            this->p_alpha,
            this->p_zeta,
            this->p_voigt_modulus_frac,
            this->p_voigt_viscosity_frac);
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, static_cast<uint32_t>(BinaryClassID::Sundberg),
            {this->p_alpha, this->p_zeta, this->p_voigt_modulus_frac, this->p_voigt_viscosity_frac});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 4);
        this->p_alpha                = params[0];
        this->p_zeta                 = params[1];
        this->p_voigt_modulus_frac   = params[2];
        this->p_voigt_viscosity_frac = params[3];
    }

protected:
    double p_alpha                 = 0.3;
    double p_zeta                  = 1.0;
    double p_voigt_modulus_frac = 5.0;
    double p_voigt_viscosity_frac  = 0.02;
};

// =====================================================================================================================
// Factory
// =====================================================================================================================

// -------------------------------------------------------------------------------
// c_RheologyModel — one named value per rheology model.
// Used by c_find_rheology to dispatch to the correct class without string
// comparisons. The Cython layer maps Python strings (and aliases) to these
// values via c_rheology_model_from_name.
// -------------------------------------------------------------------------------
enum class c_RheologyModel : uint8_t {
    Elastic  = 0,
    Viscous  = 1,
    Voigt    = 2,
    Maxwell  = 3,
    Burgers  = 4,
    Andrade  = 5,
    Sundberg = 6,
};

// -------------------------------------------------------------------------------
// c_rheology_model_from_name — map a (case-insensitive) model name or alias to a
// c_RheologyModel enum value.
//
// Recognised names and aliases:
//   "elastic" / "off"
//   "viscous" / "newton"
//   "voigt"   / "voigt-kelvin" / "voigt_kelvin"
//   "maxwell"
//   "burgers"
//   "andrade"
//   "sundberg" / "sundberg-cooper" / "sundberg_cooper"
//
// Throws std::invalid_argument if the model name is unknown.
// -------------------------------------------------------------------------------
inline c_RheologyModel c_rheology_model_from_name(const std::string& model_name) {
    const std::string name = rheo_to_lower(model_name);

    if (name == "elastic" || name == "off")          { return c_RheologyModel::Elastic; }
    if (name == "viscous" || name == "newton")       { return c_RheologyModel::Viscous; }
    if (name == "voigt"   || name == "voigt-kelvin"
                          || name == "voigt_kelvin") { return c_RheologyModel::Voigt; }
    if (name == "maxwell")                           { return c_RheologyModel::Maxwell; }
    if (name == "burgers")                           { return c_RheologyModel::Burgers; }
    if (name == "andrade")                           { return c_RheologyModel::Andrade; }
    if (name == "sundberg"
             || name == "sundberg-cooper"
             || name == "sundberg_cooper")           { return c_RheologyModel::Sundberg; }

    throw std::invalid_argument("TidalPy: unknown rheology model name '" + model_name + "'");
}

// -------------------------------------------------------------------------------
// c_find_rheology — build the rheology model named by a c_RheologyModel enum.
//
// Returns a unique_ptr to a newly heap-allocated concrete model constructed from
// the supplied config.  This is the canonical C++ factory; C++ consumers (layers
// attaching rheology, binary reconstruction, the Cython wrapper) all route
// through it.  Throws std::invalid_argument for an unrecognised enum value.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_RheologyBase> c_find_rheology(
        c_RheologyModel model, const c_RheologyConfig& cfg) {
    switch (model) {
        case c_RheologyModel::Elastic:  return std::make_unique<c_Elastic>(cfg);
        case c_RheologyModel::Viscous:  return std::make_unique<c_Viscous>(cfg);
        case c_RheologyModel::Voigt:    return std::make_unique<c_Voigt>(cfg);
        case c_RheologyModel::Maxwell:  return std::make_unique<c_Maxwell>(cfg);
        case c_RheologyModel::Burgers:  return std::make_unique<c_Burgers>(cfg);
        case c_RheologyModel::Andrade:  return std::make_unique<c_Andrade>(cfg);
        case c_RheologyModel::Sundberg: return std::make_unique<c_Sundberg>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_RheologyModel enum value");
}

// -------------------------------------------------------------------------------
// c_find_rheology (name overload) — convenience wrapper that maps a name/alias to
// the enum and builds the model.  Throws std::invalid_argument on unknown names.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_RheologyBase> c_find_rheology(
        const std::string& model_name, const c_RheologyConfig& cfg) {
    return c_find_rheology(c_rheology_model_from_name(model_name), cfg);
}

// -------------------------------------------------------------------------------
// c_rheology_from_binary — reconstruct a rheology model from a binary stream.
//
// Peeks the upcoming record's BinaryClassID (without consuming the header),
// constructs the matching default-initialised concrete model, then delegates to
// its read_binary to restore the model name and parameters. Used by the layer
// recursive deserialization (see structures_x/layers). Throws std::runtime_error
// if the class id is not a known rheology model.
// -------------------------------------------------------------------------------
inline std::unique_ptr<c_RheologyBase> c_rheology_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_RheologyBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::Elastic:  model = std::make_unique<c_Elastic>();  break;
        case BinaryClassID::Viscous:  model = std::make_unique<c_Viscous>();  break;
        case BinaryClassID::Voigt:    model = std::make_unique<c_Voigt>();    break;
        case BinaryClassID::Maxwell:  model = std::make_unique<c_Maxwell>();  break;
        case BinaryClassID::Burgers:  model = std::make_unique<c_Burgers>();  break;
        case BinaryClassID::Andrade:  model = std::make_unique<c_Andrade>();  break;
        case BinaryClassID::Sundberg: model = std::make_unique<c_Sundberg>(); break;
        default:
            throw std::runtime_error("TidalPy: unknown rheology class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
