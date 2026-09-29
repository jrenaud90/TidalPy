#pragma once
/* TidalPy's rheology models. Inputs are reference (background) MKS values at the layer mid-point.
 *
 * References
 * ----------
 * - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000
 *   (Maxwell, Voigt-Kelvin, Burgers).
 * - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150 (complex compliances and Love numbers).
 * - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784 (Andrade and Sundberg-Cooper).
 * - Zener (1948), Elasticity and Anelasticity of Metals; Nowick and Berry (1972), Anelastic Relaxation in
 *   Crystalline Solids (the standard linear solid).
 * - Kanamori and Anderson (1977), Rev. Geophys., DOI: 10.1029/RG015i001p00105; Wahr and Bergen (1986), GJRAS,
 *   DOI: 10.1111/j.1365-246X.1986.tb06642.x (a seismic Q and the dispersion it carries to tidal frequencies).
 *
 * Binary payload: model name then the model's doubles. The layer observer pointer is not serialized.
 */

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
#include "model_names_.hpp"
#include "rheology_base_.hpp"

namespace tidalpy {

// Combined construction parameters; each model reads only the fields it needs. Its defaults are the models' defaults.
struct c_RheologyConfig {
    double alpha                = 0.3;     // Andrade exponent           [dimensionless]
    double zeta                 = 1.0;     // Andrade timescale ratio    [dimensionless]
    double voigt_modulus_frac   = 5.0;     // Voigt modulus fraction     [dimensionless]
    double voigt_viscosity_frac = 0.02;    // Voigt viscosity fraction   [dimensionless]
    double relaxed_modulus_frac = 0.5;     // Zener relaxed modulus as a fraction of the unrelaxed one
    double reference_frequency  = 2.0 * TidalPyConstants::d_PI;  // Seismic Q reference frequency [rad s-1]; 1 s
    double q_frequency_exponent = 0.0;     // Seismic Q grows as frequency^exponent [dimensionless]
};

// The factors of the Andrade transient that depend on alpha alone: Gamma(1 + alpha) and the cosine and sine of
// alpha pi / 2. A model computes them once when its alpha is set, not on every call.
struct c_AndradeFactors {
    double gamma_term = 1.0;
    double cos_term   = 1.0;
    double sin_term   = 0.0;

    c_AndradeFactors() = default;
    explicit c_AndradeFactors(double alpha) noexcept :
        gamma_term(std::tgamma(1.0 + alpha)),
        cos_term(std::cos(0.5 * alpha * TidalPyConstants::d_PI)),
        sin_term(std::sin(0.5 * alpha * TidalPyConstants::d_PI)) {}
};

// Internal element compliances [Pa^-1]. The composite rheologies (Burgers, Andrade, Sundberg) put their
// elements in series, so the element compliances add and the modulus is the reciprocal of that sum.
// The Andrade family additionally assumes a positive forcing frequency.
//
// An infinite viscosity is the viscosity models' cold limit (rigid): its dashpots never move, at any frequency
// including zero, so the viscous and transient terms vanish rather than forming inf / inf or inf * 0.
namespace detail {

// Maxwell element compliance: J* = J - i / (viscosity * frequency).
inline c_ComplexCompliance element_compliance_maxwell(
        double modulus,
        double viscosity,
        double frequency) noexcept {
    const double static_compliance = 1.0 / c_guard_denominator(modulus);
    if (std::isinf(viscosity)) {
        return c_ComplexCompliance(static_compliance, 0.0);
    }
    const double denom = c_guard_denominator(viscosity * frequency);
    return c_ComplexCompliance(static_compliance, -1.0 / denom);
}

// Voigt-Kelvin element. The Voigt arm's compliance is the layer compliance divided by the modulus
// fraction: J_voigt = (1 / modulus) / voigt_modulus_frac.
inline c_ComplexCompliance element_compliance_voigt(
        double modulus,
        double viscosity,
        double frequency,
        double voigt_modulus_frac,
        double voigt_viscosity_frac) noexcept {
    const double static_compliance = 1.0 / c_guard_denominator(modulus);
    const double voigt_compliance = static_compliance / c_guard_denominator(voigt_modulus_frac);
    const double voigt_viscosity  = voigt_viscosity_frac * viscosity;

    const double scaled = voigt_compliance * voigt_viscosity * frequency;
    if (std::isinf(viscosity) || std::isinf(scaled)) {
        // A locked (or effectively locked) dashpot: the Voigt arm does not deform.
        return c_ComplexCompliance(0.0, 0.0);
    }
    const double denom  = scaled * scaled + 1.0;
    const double real_j = voigt_compliance / denom;
    const double imag_j = -(voigt_compliance * voigt_compliance) * voigt_viscosity
                           * frequency / denom;
    return c_ComplexCompliance(real_j, imag_j);
}

// Andrade element compliance: Maxwell compliance plus a transient term ~ omega^{-alpha}. factors must be
// c_AndradeFactors(alpha).
inline c_ComplexCompliance element_compliance_andrade(
        double modulus,
        double viscosity,
        double frequency,
        double alpha,
        double zeta,
        const c_AndradeFactors& factors) noexcept {
    if (std::isinf(viscosity)) {
        // The transient creep scales with the Maxwell time, so it vanishes with the viscous flow.
        return element_compliance_maxwell(modulus, viscosity, frequency);
    }
    const double static_compliance = 1.0 / c_guard_denominator(modulus);
    const double andrade_term =
        c_guard_denominator(static_compliance * viscosity * frequency * zeta);

    const double const_term =
        static_compliance * std::pow(andrade_term, -alpha) * factors.gamma_term;
    const c_ComplexCompliance andrade_transient(
        factors.cos_term * const_term,
        -factors.sin_term * const_term
    );

    return element_compliance_maxwell(modulus, viscosity, frequency)
         + andrade_transient;
}

}  // namespace detail

// Complex (shear or bulk) modulus functions [Pa]. Simple models are analytic; the series composites
// invert the sum of their element compliances.

// Elastic: mu* = modulus. No dissipation, frequency independent.
inline c_ComplexModulus rheo_modulus_elastic(
        double modulus,
        double /*viscosity*/,
        double /*frequency*/) noexcept {
    return c_ComplexModulus(modulus, 0.0);
}

// Viscous (Newton): mu* = i * viscosity * frequency (purely dissipative).
inline c_ComplexModulus rheo_modulus_viscous(
        double /*modulus*/,
        double viscosity,
        double frequency) noexcept {
    return c_ComplexModulus(0.0, viscosity * frequency);
}

// Maxwell: mu* = 1 / J_maxwell.
inline c_ComplexModulus rheo_modulus_maxwell(
        double modulus,
        double viscosity,
        double frequency) noexcept {
    return c_ComplexModulus(1.0, 0.0)
         / detail::element_compliance_maxwell(
            modulus,
            viscosity,
            frequency);
}

// Voigt-Kelvin: mu* = 1 / J_voigt = voigt_modulus_frac * modulus + i * voigt_viscosity_frac * viscosity * frequency,
// the spring and dashpot in parallel. An infinite viscosity is infinitely stiff at any nonzero frequency; at zero
// frequency only the spring is loaded.
inline c_ComplexModulus rheo_modulus_voigt(
        double modulus,
        double viscosity,
        double frequency,
        double voigt_modulus_frac,
        double voigt_viscosity_frac) noexcept {
    const double spring = c_guard_denominator(modulus) * c_guard_denominator(voigt_modulus_frac);
    const double dashpot = (frequency == 0.0) ? 0.0 : voigt_viscosity_frac * viscosity * frequency;
    return c_ComplexModulus(spring, dashpot);
}

// Burgers: Maxwell and Voigt elements in series; mu* = 1 / (J_maxwell + J_voigt).
inline c_ComplexModulus rheo_modulus_burgers(
        double modulus,
        double viscosity,
        double frequency,
        double voigt_modulus_frac,
        double voigt_viscosity_frac) noexcept {
    const c_ComplexCompliance total =
        detail::element_compliance_maxwell(
            modulus,
            viscosity,
            frequency)
      + detail::element_compliance_voigt(
        modulus,
        viscosity,
        frequency,
        voigt_modulus_frac,
        voigt_viscosity_frac);
    return c_ComplexModulus(1.0, 0.0) / total;
}

// Andrade: mu* = 1 / J_andrade. factors must be c_AndradeFactors(alpha); the form without them builds them.
inline c_ComplexModulus rheo_modulus_andrade(
        double modulus,
        double viscosity,
        double frequency,
        double alpha,
        double zeta,
        const c_AndradeFactors& factors) noexcept {
    return c_ComplexModulus(1.0, 0.0)
         / detail::element_compliance_andrade(
            modulus,
            viscosity,
            frequency,
            alpha,
            zeta,
            factors);
}
inline c_ComplexModulus rheo_modulus_andrade(
        double modulus,
        double viscosity,
        double frequency,
        double alpha,
        double zeta) noexcept {
    return rheo_modulus_andrade(modulus, viscosity, frequency, alpha, zeta, c_AndradeFactors(alpha));
}

// Sundberg-Cooper: Andrade and Voigt elements in series; mu* = 1 / (J_andrade + J_voigt). factors must be
// c_AndradeFactors(alpha); the form without them builds them.
inline c_ComplexModulus rheo_modulus_sundberg(
        double modulus,
        double viscosity,
        double frequency,
        double alpha,
        double zeta,
        double voigt_modulus_frac,
        double voigt_viscosity_frac,
        const c_AndradeFactors& factors) noexcept {
    const c_ComplexCompliance total =
        detail::element_compliance_andrade(
            modulus,
            viscosity,
            frequency,
            alpha,
            zeta,
            factors)
      + detail::element_compliance_voigt(
            modulus,
            viscosity,
            frequency,
            voigt_modulus_frac,
            voigt_viscosity_frac);
    return c_ComplexModulus(1.0, 0.0) / total;
}
inline c_ComplexModulus rheo_modulus_sundberg(
        double modulus,
        double viscosity,
        double frequency,
        double alpha,
        double zeta,
        double voigt_modulus_frac,
        double voigt_viscosity_frac) noexcept {
    return rheo_modulus_sundberg(
        modulus, viscosity, frequency, alpha, zeta, voigt_modulus_frac, voigt_viscosity_frac,
        c_AndradeFactors(alpha));
}

// Zener (standard linear solid): a spring of the relaxed modulus r * modulus in parallel with a Maxwell arm whose
// spring is the rest, (1 - r) * modulus, and whose dashpot is the viscosity. With tau = viscosity / ((1 - r) modulus),
//     mu* = r modulus + (1 - r) modulus * i omega tau / (1 + i omega tau).
// The unrelaxed modulus (high frequency, or a locked dashpot) is the modulus; the relaxed one (zero frequency) is
// r * modulus, where a Maxwell body relaxes to zero. r = 0 is exactly Maxwell and r = 1 is elastic. The loss peaks at
// omega tau = 1. The forms in x = omega tau and in 1 / x keep either extreme free of inf / inf.
inline c_ComplexModulus rheo_modulus_zener(
        double modulus,
        double viscosity,
        double frequency,
        double relaxed_modulus_frac) noexcept {
    const double relaxed = relaxed_modulus_frac * modulus;
    const double arm     = modulus - relaxed;
    if (arm == 0.0 || std::isinf(viscosity)) {
        return c_ComplexModulus(modulus, 0.0);
    }
    if (frequency == 0.0) {
        return c_ComplexModulus(relaxed, 0.0);
    }
    const double x = viscosity * frequency / arm;
    if (std::abs(x) >= 1.0) {
        const double s = 1.0 / x;
        const double denom = 1.0 + s * s;
        return c_ComplexModulus(relaxed + arm / denom, arm * s / denom);
    }
    const double denom = 1.0 + x * x;
    return c_ComplexModulus(relaxed + arm * x * x / denom, arm * x / denom);
}

// Seismic Q: the complex modulus from a quality factor measured at a reference frequency, with no viscosity. This
// model reads its viscosity input as that quality factor, Q_ref. With s = reference_frequency / |frequency| and the
// exponent a in [0, 1),
//     Q(omega) = Q_ref s^-a                                   (a = 0: one Q at every frequency)
//     Re mu*   = modulus / (1 + D(s) / Q_ref),   D(s) = cot(a pi / 2) (s^a - 1)   (a = 0: D = (2 / pi) ln s)
//     Im mu*   = Re mu* / Q(omega)
// D is the dispersion causality ties to that loss (Kramers-Kronig, first order in 1 / Q): even a constant Q softens
// the modulus logarithmically toward low frequency. The modulus is the one at the reference frequency, so a seismic
// profile's moduli go in as given. Written as a compliance the storage modulus stays positive; to first order in
// 1 / Q it is the Kanamori-Anderson (a = 0) and Wahr-Bergen (a > 0) form. dispersion_coefficient must be
// cot(a pi / 2), or 2 / pi at a = 0. Zero frequency or an infinite Q is unforced: the modulus, with no loss. A Q that
// is not positive, or a frequency so far above the reference that the dispersion drives the modulus through zero,
// gives NaN.
inline c_ComplexModulus rheo_modulus_seismic_q(
        double modulus,
        double quality_factor,
        double frequency,
        double reference_frequency,
        double q_frequency_exponent,
        double dispersion_coefficient) noexcept {
    if (frequency == 0.0 || std::isinf(quality_factor)) {
        return c_ComplexModulus(modulus, 0.0);
    }
    if (!(quality_factor > 0.0)) {
        return c_ComplexModulus(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
    }
    const double log_s = std::log(reference_frequency / std::abs(frequency));
    // expm1 keeps D continuous as the exponent goes to zero, where cot(a pi / 2) (s^a - 1) -> (2 / pi) ln s.
    const double dispersion = (q_frequency_exponent == 0.0)
        ? dispersion_coefficient * log_s
        : dispersion_coefficient * std::expm1(q_frequency_exponent * log_s);
    const double denom = 1.0 + dispersion / quality_factor;
    if (!(denom > 0.0)) {
        return c_ComplexModulus(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
    }
    const double storage = modulus / denom;
    const double loss_over_q = std::exp(q_frequency_exponent * log_s) / quality_factor;
    // Odd in frequency, like the viscous models: the loss changes sign with the forcing.
    return c_ComplexModulus(storage, std::copysign(storage * loss_over_q, frequency));
}

// Each model supplies only its BinaryClassID (get_binary_class_id) and its scalar params; c_PhysicsBase handles the
// header, the model name, and the byte layout.

// Purely elastic response (alias "off").
class c_Elastic final : public c_RheologyBase {
public:
    c_Elastic() : c_Elastic(c_RheologyConfig{}) {}
    explicit c_Elastic(const c_RheologyConfig& /*cfg*/) : c_RheologyBase("elastic") {}
    ~c_Elastic() override = default;

    c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const override {
        return rheo_modulus_elastic(
            modulus, 
            viscosity,
            frequency);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Elastic); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id());
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// Purely viscous response (alias "newton").
class c_Viscous final : public c_RheologyBase {
public:
    c_Viscous() : c_Viscous(c_RheologyConfig{}) {}
    explicit c_Viscous(const c_RheologyConfig& /*cfg*/) : c_RheologyBase("viscous") {}
    ~c_Viscous() override = default;

    c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const override {
        return rheo_modulus_viscous(
            modulus,
            viscosity,
            frequency);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Viscous); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id());
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

class c_Maxwell final : public c_RheologyBase {
public:
    c_Maxwell() : c_Maxwell(c_RheologyConfig{}) {}
    explicit c_Maxwell(const c_RheologyConfig& /*cfg*/) : c_RheologyBase("maxwell") {}
    ~c_Maxwell() override = default;

    c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const override {
        return rheo_modulus_maxwell(
            modulus,
            viscosity,
            frequency);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Maxwell); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id());
    }
    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }
};

// Voigt-Kelvin element (alias "voigt-kelvin").
class c_Voigt final : public c_RheologyBase {
public:
    c_Voigt() : c_Voigt(c_RheologyConfig{}) {}
    explicit c_Voigt(const c_RheologyConfig& cfg)
        : c_RheologyBase("voigt"),
          p_voigt_modulus_frac(cfg.voigt_modulus_frac),
          p_voigt_viscosity_frac(cfg.voigt_viscosity_frac) {}
    ~c_Voigt() override = default;

    double get_voigt_modulus_frac() const noexcept { return this->p_voigt_modulus_frac; }
    double get_voigt_viscosity_frac() const noexcept { return this->p_voigt_viscosity_frac; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RheologyBase::append_config_entries(out);
        out.push_back(c_config_double("voigt_modulus_frac", this->p_voigt_modulus_frac));
        out.push_back(c_config_double("voigt_viscosity_frac", this->p_voigt_viscosity_frac));
    }

    c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const override {
        return rheo_modulus_voigt(
            modulus,
            viscosity,
            frequency,
            this->p_voigt_modulus_frac,
            this->p_voigt_viscosity_frac);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Voigt); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id(),
                                   {this->p_voigt_modulus_frac, this->p_voigt_viscosity_frac});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 2);
        this->p_voigt_modulus_frac   = params[0];
        this->p_voigt_viscosity_frac = params[1];
    }

protected:
    double p_voigt_modulus_frac;
    double p_voigt_viscosity_frac;
};

// Maxwell and Voigt in series.
class c_Burgers final : public c_RheologyBase {
public:
    c_Burgers() : c_Burgers(c_RheologyConfig{}) {}
    explicit c_Burgers(const c_RheologyConfig& cfg)
        : c_RheologyBase("burgers"),
          p_voigt_modulus_frac(cfg.voigt_modulus_frac),
          p_voigt_viscosity_frac(cfg.voigt_viscosity_frac) {}
    ~c_Burgers() override = default;

    double get_voigt_modulus_frac()   const noexcept { return this->p_voigt_modulus_frac; }
    double get_voigt_viscosity_frac() const noexcept { return this->p_voigt_viscosity_frac; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RheologyBase::append_config_entries(out);
        out.push_back(c_config_double("voigt_modulus_frac", this->p_voigt_modulus_frac));
        out.push_back(c_config_double("voigt_viscosity_frac", this->p_voigt_viscosity_frac));
    }

    c_ComplexModulus calc_complex_modulus(
            double modulus, 
            double viscosity,
            double frequency) const override {
        return rheo_modulus_burgers(
            modulus,
            viscosity,
            frequency,
            this->p_voigt_modulus_frac,
            this->p_voigt_viscosity_frac);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Burgers); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id(),
                                   {this->p_voigt_modulus_frac, this->p_voigt_viscosity_frac});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 2);
        this->p_voigt_modulus_frac   = params[0];
        this->p_voigt_viscosity_frac = params[1];
    }

protected:
    double p_voigt_modulus_frac;
    double p_voigt_viscosity_frac;
};

// Maxwell plus an Andrade transient term.
class c_Andrade final : public c_RheologyBase {
public:
    c_Andrade() : c_Andrade(c_RheologyConfig{}) {}
    explicit c_Andrade(const c_RheologyConfig& cfg)
        : c_RheologyBase("andrade"),
          p_alpha(cfg.alpha),
          p_zeta(cfg.zeta) {}
    ~c_Andrade() override = default;

    double get_alpha() const noexcept { return this->p_alpha; }
    double get_zeta()  const noexcept { return this->p_zeta; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RheologyBase::append_config_entries(out);
        out.push_back(c_config_double("alpha", this->p_alpha));
        out.push_back(c_config_double("zeta", this->p_zeta));
    }

    c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const override {
        return rheo_modulus_andrade(
            modulus,
            viscosity,
            frequency,
            this->p_alpha,
            this->p_zeta,
            this->p_andrade_factors);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Andrade); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id(),
                                   {this->p_alpha, this->p_zeta});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 2);
        this->p_alpha = params[0];
        this->p_zeta  = params[1];
        this->p_andrade_factors = c_AndradeFactors(this->p_alpha);
    }

protected:
    double p_alpha;
    double p_zeta;
    // Declared after p_alpha, so every constructor builds it from the alpha it set.
    c_AndradeFactors p_andrade_factors{this->p_alpha};
};

// Andrade and Voigt (alias "sundberg-cooper").
class c_Sundberg final : public c_RheologyBase {
public:
    c_Sundberg() : c_Sundberg(c_RheologyConfig{}) {}
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

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RheologyBase::append_config_entries(out);
        out.push_back(c_config_double("alpha", this->p_alpha));
        out.push_back(c_config_double("zeta", this->p_zeta));
        out.push_back(c_config_double("voigt_modulus_frac", this->p_voigt_modulus_frac));
        out.push_back(c_config_double("voigt_viscosity_frac", this->p_voigt_viscosity_frac));
    }

    c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const override {
        return rheo_modulus_sundberg(
            modulus, 
            viscosity,
            frequency,
            this->p_alpha,
            this->p_zeta,
            this->p_voigt_modulus_frac,
            this->p_voigt_viscosity_frac,
            this->p_andrade_factors);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Sundberg); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(
            out, this->get_binary_class_id(),
            {this->p_alpha, this->p_zeta, this->p_voigt_modulus_frac, this->p_voigt_viscosity_frac});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 4);
        this->p_alpha                = params[0];
        this->p_zeta                 = params[1];
        this->p_voigt_modulus_frac   = params[2];
        this->p_voigt_viscosity_frac = params[3];
        this->p_andrade_factors      = c_AndradeFactors(this->p_alpha);
    }

protected:
    double p_alpha;
    double p_zeta;
    double p_voigt_modulus_frac;
    double p_voigt_viscosity_frac;
    // Declared after p_alpha, so every constructor builds it from the alpha it set.
    c_AndradeFactors p_andrade_factors{this->p_alpha};
};

// Zener, the standard linear solid (aliases "sls", "standard_linear_solid"). Unlike Maxwell it relaxes to a finite
// modulus, which suits a bulk response: melt-driven compaction relaxes a partially molten rock's bulk modulus toward
// its drained value, not to zero.
class c_Zener final : public c_RheologyBase {
public:
    c_Zener() : c_Zener(c_RheologyConfig{}) {}
    explicit c_Zener(const c_RheologyConfig& cfg)
        : c_RheologyBase("zener"),
          p_relaxed_modulus_frac(cfg.relaxed_modulus_frac) {
        c_check_relaxed_modulus_frac(this->p_relaxed_modulus_frac);
    }
    ~c_Zener() override = default;

    double get_relaxed_modulus_frac() const noexcept { return this->p_relaxed_modulus_frac; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RheologyBase::append_config_entries(out);
        out.push_back(c_config_double("relaxed_modulus_frac", this->p_relaxed_modulus_frac));
    }

    c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const override {
        return rheo_modulus_zener(
            modulus,
            viscosity,
            frequency,
            this->p_relaxed_modulus_frac);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::Zener); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id(), {this->p_relaxed_modulus_frac});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 1);
        c_check_relaxed_modulus_frac(params[0]);
        this->p_relaxed_modulus_frac = params[0];
    }

protected:
    double p_relaxed_modulus_frac;

    // A fraction outside [0, 1] would give a negative spring.
    static void c_check_relaxed_modulus_frac(double value) {
        if (!(value >= 0.0 && value <= 1.0)) {
            throw std::invalid_argument(
                "TidalPy: the Zener relaxed_modulus_frac must lie in [0, 1], got " + std::to_string(value));
        }
    }
};

// Seismic Q (aliases "constant_q", "power_law_q"): the loss comes from a quality factor rather than a viscosity. Its
// viscosity input is read as Q at the reference frequency, which is how a seismic profile's Q(r) reaches the Love
// solve; see rheo_modulus_seismic_q.
class c_SeismicQ final : public c_RheologyBase {
public:
    c_SeismicQ() : c_SeismicQ(c_RheologyConfig{}) {}
    explicit c_SeismicQ(const c_RheologyConfig& cfg)
        : c_RheologyBase("seismic_q"),
          p_reference_frequency(cfg.reference_frequency),
          p_q_frequency_exponent(cfg.q_frequency_exponent) {
        c_check_params(this->p_reference_frequency, this->p_q_frequency_exponent);
        this->p_set_dispersion_coefficient();
    }
    ~c_SeismicQ() override = default;

    double get_reference_frequency()  const noexcept { return this->p_reference_frequency; }
    double get_q_frequency_exponent() const noexcept { return this->p_q_frequency_exponent; }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        c_RheologyBase::append_config_entries(out);
        out.push_back(c_config_double("reference_frequency_rad_s", this->p_reference_frequency));
        out.push_back(c_config_double("q_frequency_exponent", this->p_q_frequency_exponent));
    }

    c_ComplexModulus calc_complex_modulus(
            double modulus,
            double viscosity,
            double frequency) const override {
        return rheo_modulus_seismic_q(
            modulus,
            viscosity,
            frequency,
            this->p_reference_frequency,
            this->p_q_frequency_exponent,
            this->p_dispersion_coefficient);
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::SeismicQ); }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, this->get_binary_class_id(),
                                   {this->p_reference_frequency, this->p_q_frequency_exponent});
    }
    void read_binary(std::istream& in, bool force = false) override {
        const std::vector<double> params = this->read_physics_binary(in, force, 2);
        c_check_params(params[0], params[1]);
        this->p_reference_frequency  = params[0];
        this->p_q_frequency_exponent = params[1];
        this->p_set_dispersion_coefficient();
    }

protected:
    double p_reference_frequency;
    double p_q_frequency_exponent;
    double p_dispersion_coefficient = 0.0;

    // cot(a pi / 2), or its a -> 0 limit taken with ln s in place of s^a - 1.
    void p_set_dispersion_coefficient() noexcept {
        this->p_dispersion_coefficient = (this->p_q_frequency_exponent == 0.0)
            ? 2.0 / TidalPyConstants::d_PI
            : 1.0 / std::tan(0.5 * this->p_q_frequency_exponent * TidalPyConstants::d_PI);
    }

    // An exponent of 1 or more has no finite dispersion (cot(a pi / 2) <= 0), and a negative one a Q that falls with
    // frequency, which no absorption band gives.
    static void c_check_params(double reference_frequency, double q_frequency_exponent) {
        if (!(reference_frequency > 0.0) || std::isinf(reference_frequency)) {
            throw std::invalid_argument(
                "TidalPy: the seismic Q reference_frequency_rad_s must be positive and finite, got "
                + std::to_string(reference_frequency));
        }
        if (!(q_frequency_exponent >= 0.0 && q_frequency_exponent < 1.0)) {
            throw std::invalid_argument(
                "TidalPy: the seismic Q q_frequency_exponent must lie in [0, 1), got "
                + std::to_string(q_frequency_exponent));
        }
    }
};

// One value per model, so c_find_rheology dispatches without string comparisons.
enum class c_RheologyModel : uint8_t {
    Elastic  = 0,
    Viscous  = 1,
    Voigt    = 2,
    Maxwell  = 3,
    Burgers  = 4,
    Andrade  = 5,
    Sundberg = 6,
    Zener    = 7,
    SeismicQ = 8,
};

// Model names are matched case-insensitively.
inline c_RheologyModel c_rheology_model_from_name(const std::string& model_name) {
    const std::string name = c_to_lower(model_name);

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
    if (name == "zener"    || name == "sls"
             || name == "standard_linear_solid")     { return c_RheologyModel::Zener; }
    if (name == "seismic_q" || name == "constant_q"
             || name == "power_law_q")               { return c_RheologyModel::SeismicQ; }

    throw std::invalid_argument("TidalPy: unknown rheology model name '" + model_name + "'");
}

// Builds a model from its enum value and parameters; the Cython wrappers construct through it. A saved record is
// restored by c_rheology_from_binary instead.
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
        case c_RheologyModel::Zener:    return std::make_unique<c_Zener>(cfg);
        case c_RheologyModel::SeismicQ: return std::make_unique<c_SeismicQ>(cfg);
    }
    throw std::invalid_argument("TidalPy: unrecognised c_RheologyModel enum value");
}

inline std::unique_ptr<c_RheologyBase> c_find_rheology(
        const std::string& model_name, const c_RheologyConfig& cfg) {
    return c_find_rheology(c_rheology_model_from_name(model_name), cfg);
}

// The class id is peeked without consuming the header so the default-constructed model restores itself.
// Used by the layer recursive deserialization in Structures/layers.
inline std::unique_ptr<c_RheologyBase> c_rheology_from_binary(std::istream& in, bool force = false) {
    const c_BinaryHeader header = c_peek_binary_header(in);

    std::unique_ptr<c_RheologyBase> model;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::Elastic:  model = std::make_unique<c_Elastic>();  break;
        case BinaryClassID::Viscous:  model = std::make_unique<c_Viscous>();  break;
        case BinaryClassID::Voigt:    model = std::make_unique<c_Voigt>();    break;
        case BinaryClassID::Maxwell:  model = std::make_unique<c_Maxwell>();  break;
        case BinaryClassID::Burgers:  model = std::make_unique<c_Burgers>();  break;
        case BinaryClassID::Andrade:  model = std::make_unique<c_Andrade>();  break;
        case BinaryClassID::Sundberg: model = std::make_unique<c_Sundberg>(); break;
        case BinaryClassID::Zener:    model = std::make_unique<c_Zener>();    break;
        case BinaryClassID::SeismicQ: model = std::make_unique<c_SeismicQ>(); break;
        default:
            throw std::runtime_error("TidalPy: unknown rheology class id in binary stream");
    }
    model->read_binary(in, force);
    return model;
}

}  // namespace tidalpy
