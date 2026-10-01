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
 * Each model's parameters are declared once in its table (c_SpecModel), which gives its config entries and
 * binary record.
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
#include "registry_.hpp"
#include "rheology_base_.hpp"
#include "spec_model_.hpp"

namespace tidalpy {

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

// Each model declares its parameters in one table (c_SpecModel, spec_model_.hpp), which gives its construction,
// validation, config entries, and binary record. The physics is in the rheo_modulus_* functions above.

// The parameter rows several models share, so each is described once.
template <class Model>
inline c_ParamSpec<Model> c_voigt_modulus_frac_spec(double Model::* member) {
    return {"voigt_modulus_frac", "voigt_modulus_frac", member, 5.0, c_ParamBounds::Positive,
            "Voigt spring as a multiple of the unrelaxed modulus [dimensionless]."};
}
template <class Model>
inline c_ParamSpec<Model> c_voigt_viscosity_frac_spec(double Model::* member) {
    return {"voigt_viscosity_frac", "voigt_viscosity_frac", member, 0.02, c_ParamBounds::Positive,
            "Voigt dashpot as a fraction of the viscosity [dimensionless]."};
}
template <class Model>
inline c_ParamSpec<Model> c_andrade_alpha_spec(double Model::* member) {
    return {"alpha", "alpha", member, 0.3, c_ParamBounds::UnitInterval,
            "Andrade exponent, in (0, 1] [dimensionless]."};
}
template <class Model>
inline c_ParamSpec<Model> c_andrade_zeta_spec(double Model::* member) {
    return {"zeta", "zeta", member, 1.0, c_ParamBounds::Positive,
            "Andrade timescale over the Maxwell time [dimensionless]."};
}

// The empty table of a model with no parameters.
template <class Model>
inline const std::vector<c_ParamSpec<Model>>& c_no_param_specs() {
    static const std::vector<c_ParamSpec<Model>> specs;
    return specs;
}

// Purely elastic response (alias "off").
class c_Elastic final : public c_SpecModel<c_Elastic, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Elastic;
    static const std::vector<c_ParamSpec<c_Elastic>>& parameter_specs() { return c_no_param_specs<c_Elastic>(); }

    c_Elastic() : c_Elastic(c_ParamMap{}) {}
    explicit c_Elastic(const c_ParamMap& params) : c_SpecModel("elastic") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_elastic(modulus, viscosity, frequency);
    }
};

// Purely viscous response (alias "newton").
class c_Viscous final : public c_SpecModel<c_Viscous, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Viscous;
    static const std::vector<c_ParamSpec<c_Viscous>>& parameter_specs() { return c_no_param_specs<c_Viscous>(); }

    c_Viscous() : c_Viscous(c_ParamMap{}) {}
    explicit c_Viscous(const c_ParamMap& params) : c_SpecModel("viscous") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_viscous(modulus, viscosity, frequency);
    }
};

class c_Maxwell final : public c_SpecModel<c_Maxwell, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Maxwell;
    static const std::vector<c_ParamSpec<c_Maxwell>>& parameter_specs() { return c_no_param_specs<c_Maxwell>(); }

    c_Maxwell() : c_Maxwell(c_ParamMap{}) {}
    explicit c_Maxwell(const c_ParamMap& params) : c_SpecModel("maxwell") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_maxwell(modulus, viscosity, frequency);
    }
};

// Voigt-Kelvin element (alias "voigt-kelvin").
class c_Voigt final : public c_SpecModel<c_Voigt, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Voigt;
    static const std::vector<c_ParamSpec<c_Voigt>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_Voigt>> specs = {
            c_voigt_modulus_frac_spec<c_Voigt>(&c_Voigt::p_voigt_modulus_frac),
            c_voigt_viscosity_frac_spec<c_Voigt>(&c_Voigt::p_voigt_viscosity_frac),
        };
        return specs;
    }

    c_Voigt() : c_Voigt(c_ParamMap{}) {}
    explicit c_Voigt(const c_ParamMap& params) : c_SpecModel("voigt") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_voigt(
            modulus, viscosity, frequency, this->p_voigt_modulus_frac, this->p_voigt_viscosity_frac);
    }

protected:
    double p_voigt_modulus_frac   = 0.0;
    double p_voigt_viscosity_frac = 0.0;
};

// Maxwell and Voigt in series.
class c_Burgers final : public c_SpecModel<c_Burgers, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Burgers;
    static const std::vector<c_ParamSpec<c_Burgers>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_Burgers>> specs = {
            c_voigt_modulus_frac_spec<c_Burgers>(&c_Burgers::p_voigt_modulus_frac),
            c_voigt_viscosity_frac_spec<c_Burgers>(&c_Burgers::p_voigt_viscosity_frac),
        };
        return specs;
    }

    c_Burgers() : c_Burgers(c_ParamMap{}) {}
    explicit c_Burgers(const c_ParamMap& params) : c_SpecModel("burgers") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_burgers(
            modulus, viscosity, frequency, this->p_voigt_modulus_frac, this->p_voigt_viscosity_frac);
    }

protected:
    double p_voigt_modulus_frac   = 0.0;
    double p_voigt_viscosity_frac = 0.0;
};

// Maxwell plus an Andrade transient term.
class c_Andrade final : public c_SpecModel<c_Andrade, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Andrade;
    static const std::vector<c_ParamSpec<c_Andrade>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_Andrade>> specs = {
            c_andrade_alpha_spec<c_Andrade>(&c_Andrade::p_alpha),
            c_andrade_zeta_spec<c_Andrade>(&c_Andrade::p_zeta),
        };
        return specs;
    }

    c_Andrade() : c_Andrade(c_ParamMap{}) {}
    explicit c_Andrade(const c_ParamMap& params) : c_SpecModel("andrade") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_andrade(modulus, viscosity, frequency, this->p_alpha, this->p_zeta, this->p_factors);
    }

protected:
    // A zero exponent has no transient creep, which is the Maxwell model.
    void p_validate() const override {
        if (!(this->p_alpha > 0.0)) {
            throw std::invalid_argument(this->p_describe() + " needs an 'alpha' above 0; maxwell has none.");
        }
    }
    void p_update_derived() noexcept override { this->p_factors = c_AndradeFactors(this->p_alpha); }

    double p_alpha = 0.0;
    double p_zeta  = 0.0;
    c_AndradeFactors p_factors;
};

// Andrade and Voigt in series (alias "sundberg-cooper").
class c_Sundberg final : public c_SpecModel<c_Sundberg, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Sundberg;
    static const std::vector<c_ParamSpec<c_Sundberg>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_Sundberg>> specs = {
            c_andrade_alpha_spec<c_Sundberg>(&c_Sundberg::p_alpha),
            c_andrade_zeta_spec<c_Sundberg>(&c_Sundberg::p_zeta),
            c_voigt_modulus_frac_spec<c_Sundberg>(&c_Sundberg::p_voigt_modulus_frac),
            c_voigt_viscosity_frac_spec<c_Sundberg>(&c_Sundberg::p_voigt_viscosity_frac),
        };
        return specs;
    }

    c_Sundberg() : c_Sundberg(c_ParamMap{}) {}
    explicit c_Sundberg(const c_ParamMap& params) : c_SpecModel("sundberg") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_sundberg(
            modulus, viscosity, frequency, this->p_alpha, this->p_zeta, this->p_voigt_modulus_frac,
            this->p_voigt_viscosity_frac, this->p_factors);
    }

protected:
    void p_validate() const override {
        if (!(this->p_alpha > 0.0)) {
            throw std::invalid_argument(this->p_describe() + " needs an 'alpha' above 0; burgers has none.");
        }
    }
    void p_update_derived() noexcept override { this->p_factors = c_AndradeFactors(this->p_alpha); }

    double p_alpha                = 0.0;
    double p_zeta                 = 0.0;
    double p_voigt_modulus_frac   = 0.0;
    double p_voigt_viscosity_frac = 0.0;
    c_AndradeFactors p_factors;
};

// Zener, the standard linear solid (aliases "sls", "standard_linear_solid"). Unlike Maxwell it relaxes to a finite
// modulus, which suits a bulk response: melt-driven compaction relaxes a partially molten rock's bulk modulus toward
// its drained value, not to zero.
class c_Zener final : public c_SpecModel<c_Zener, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::Zener;
    static const std::vector<c_ParamSpec<c_Zener>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_Zener>> specs = {
            {"relaxed_modulus_frac", "relaxed_modulus_frac", &c_Zener::p_relaxed_modulus_frac, 0.5,
             c_ParamBounds::UnitInterval,
             "Relaxed (zero-frequency) modulus as a fraction of the unrelaxed one [dimensionless]."},
        };
        return specs;
    }

    c_Zener() : c_Zener(c_ParamMap{}) {}
    explicit c_Zener(const c_ParamMap& params) : c_SpecModel("zener") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_zener(modulus, viscosity, frequency, this->p_relaxed_modulus_frac);
    }

protected:
    double p_relaxed_modulus_frac = 0.0;
};

// Seismic Q (aliases "constant_q", "power_law_q"): the loss comes from a quality factor rather than a viscosity. Its
// viscosity input is read as Q at the reference frequency, which is how a seismic profile's Q(r) reaches the Love
// solve; see rheo_modulus_seismic_q.
class c_SeismicQ final : public c_SpecModel<c_SeismicQ, c_RheologyBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::SeismicQ;
    static const std::vector<c_ParamSpec<c_SeismicQ>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_SeismicQ>> specs = {
            {"reference_frequency", "reference_frequency_rad_s", &c_SeismicQ::p_reference_frequency,
             2.0 * TidalPyConstants::d_PI, c_ParamBounds::Positive,
             "Frequency at which the quality factor is given [rad s-1]; 2 pi is a 1 s period."},
            {"q_frequency_exponent", "q_frequency_exponent", &c_SeismicQ::p_q_frequency_exponent, 0.0,
             c_ParamBounds::UnitInterval, "Exponent a of Q ~ omega^a, in [0, 1) [dimensionless]."},
        };
        return specs;
    }

    c_SeismicQ() : c_SeismicQ(c_ParamMap{}) {}
    explicit c_SeismicQ(const c_ParamMap& params) : c_SpecModel("seismic_q") { this->p_initialize(params); }

    c_ComplexModulus calc_complex_modulus(double modulus, double viscosity, double frequency) const override {
        return rheo_modulus_seismic_q(
            modulus, viscosity, frequency, this->p_reference_frequency, this->p_q_frequency_exponent,
            this->p_dispersion_coefficient);
    }

protected:
    // An exponent of 1 has no finite dispersion (cot(a pi / 2) = 0).
    void p_validate() const override {
        if (!(this->p_q_frequency_exponent < 1.0)) {
            throw std::invalid_argument(this->p_describe() + " needs a 'q_frequency_exponent' below 1.");
        }
    }
    // cot(a pi / 2), or its a -> 0 limit taken with ln s in place of s^a - 1.
    void p_update_derived() noexcept override {
        this->p_dispersion_coefficient = (this->p_q_frequency_exponent == 0.0)
            ? 2.0 / TidalPyConstants::d_PI
            : 1.0 / std::tan(0.5 * this->p_q_frequency_exponent * TidalPyConstants::d_PI);
    }

    double p_reference_frequency    = 0.0;
    double p_q_frequency_exponent   = 0.0;
    double p_dispersion_coefficient = 0.0;
};

inline const c_ModelRegistry<c_RheologyBase>& c_rheology_registry() {
    static const c_ModelRegistry<c_RheologyBase> registry = {
        {{"elastic", "off"},                                 BinaryClassID::Elastic,  &c_make_entry<c_RheologyBase, c_Elastic>},
        {{"viscous", "newton"},                              BinaryClassID::Viscous,  &c_make_entry<c_RheologyBase, c_Viscous>},
        {{"voigt", "voigt-kelvin", "voigt_kelvin"},          BinaryClassID::Voigt,    &c_make_entry<c_RheologyBase, c_Voigt>},
        {{"maxwell"},                                        BinaryClassID::Maxwell,  &c_make_entry<c_RheologyBase, c_Maxwell>},
        {{"burgers"},                                        BinaryClassID::Burgers,  &c_make_entry<c_RheologyBase, c_Burgers>},
        {{"andrade"},                                        BinaryClassID::Andrade,  &c_make_entry<c_RheologyBase, c_Andrade>},
        {{"sundberg", "sundberg-cooper", "sundberg_cooper"}, BinaryClassID::Sundberg, &c_make_entry<c_RheologyBase, c_Sundberg>},
        {{"zener", "sls", "standard_linear_solid"},          BinaryClassID::Zener,    &c_make_entry<c_RheologyBase, c_Zener>},
        {{"seismic_q", "constant_q", "power_law_q"},         BinaryClassID::SeismicQ, &c_make_entry<c_RheologyBase, c_SeismicQ>},
    };
    return registry;
}

// The family's entry points, each one line over the generic registry functions.
inline std::unique_ptr<c_RheologyBase> c_find_rheology(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_rheology_registry(), model_name, params);
}

inline std::unique_ptr<c_RheologyBase> c_rheology_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_rheology_registry(), in, force);
}

inline std::string c_rheology_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_rheology_registry(), model_name);
}

inline std::vector<std::string> c_rheology_model_names() {
    return c_model_names(c_rheology_registry());
}

}  // namespace tidalpy
