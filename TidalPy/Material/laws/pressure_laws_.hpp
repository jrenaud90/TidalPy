#pragma once
/* Analytic pressure laws P(eta), eta = rho / rho0, and the density-from-pressure inversion the equation-of-state laws
 * share: the 3rd-order Birch-Murnaghan and Vinet laws, the compression range over which each rises, and a
 * safeguarded Newton inversion. All MKS.
 *
 * References
 * ----------
 * Birch (1947), Phys. Rev. 71, 809. Vinet et al. (1987), J. Geophys. Res. 92, 9319.
 */

#include <cmath>
#include <cstdint>
#include <limits>

#include "../../constants_.hpp"
#include "../../Utilities/math/numerics_.hpp"  // c_safe_pow, c_safe_exp

namespace tidalpy {

// Safeguarded Newton/bisection density-from-pressure inversion. Its tolerance and iteration cap come from
// TidalPy.config['numerical'] ('eos_invert_rtol', 'eos_invert_max_iters') when a model is built without its
// own; the fallbacks cover a model built before the config is loaded. The cap only guarantees termination;
// convergence normally takes well under 10 iterations.
inline constexpr double d_EOS_INVERT_RTOL_FALLBACK      = 1.0e-13;
inline constexpr int    d_EOS_INVERT_MAX_ITERS_FALLBACK = 60;

// The inversion tolerance a model uses: its own when positive and finite, otherwise the configured default.
inline double c_resolve_eos_invert_rtol(double invert_rtol) noexcept
{
    if (std::isfinite(invert_rtol) && invert_rtol > 0.0)
    {
        return invert_rtol;
    }
    if (tidalpy_config_ptr != nullptr && std::isfinite(tidalpy_config_ptr->d_EOS_INVERT_RTOL))
    {
        return tidalpy_config_ptr->d_EOS_INVERT_RTOL;
    }
    return d_EOS_INVERT_RTOL_FALLBACK;
}

// The inversion iteration cap a model uses: its own when positive, otherwise the configured default.
inline int c_resolve_eos_invert_max_iters(int invert_max_iters) noexcept
{
    if (invert_max_iters > 0)
    {
        return invert_max_iters;
    }
    if (tidalpy_config_ptr != nullptr && tidalpy_config_ptr->d_EOS_INVERT_MAX_ITERS > 0)
    {
        return tidalpy_config_ptr->d_EOS_INVERT_MAX_ITERS;
    }
    return d_EOS_INVERT_MAX_ITERS_FALLBACK;
}

// Ambient, where mineral-physics rho0 and K0 are quoted.
inline constexpr double d_EOS_REFERENCE_TEMPERATURE = 300.0;

// The analytic pressure laws.
enum class c_PressureLaw : uint8_t {
    BirchMurnaghan = 0,
    Vinet          = 1,
};

// Analytic pressure laws and the density-from-pressure inversion.
//
// All laws are written in the compression ratio eta = rho / rho0 = V0 / V. They rise monotonically in eta
// only near eta = 1, since the finite-strain corrections turn them over at extreme eta, so the inversion
// brackets its root within the monotonic range.

// 3rd-order Birch-Murnaghan pressure [Pa] and isothermal bulk modulus K = eta dP/deta [Pa]. One cube root
// serves every fractional power, and the inversion wants both values at once.
inline void eos_bm_pressure_and_bulk_modulus(
        double eta,
        double K0,
        double K0_prime,
        double& pressure,
        double& bulk_modulus) noexcept {
    const double cbrt_eta     = std::cbrt(eta);
    const double eta_23       = cbrt_eta * cbrt_eta;
    const double eta_53       = eta * eta_23;
    const double eta_73       = eta_53 * eta_23;
    const double strain_coeff = 0.75 * (K0_prime - 4.0);
    const double strain_term  = 1.0 + strain_coeff * (eta_23 - 1.0);
    pressure     = 1.5 * K0 * (eta_73 - eta_53) * strain_term;
    bulk_modulus = 1.5 * K0 * (
        ((7.0 / 3.0) * eta_73 - (5.0 / 3.0) * eta_53) * strain_term
        + (eta_73 - eta_53) * strain_coeff * (2.0 / 3.0) * eta_23);
}

// Vinet pressure [Pa] and isothermal bulk modulus K = eta dP/deta [Pa]; inv_cbrt_eta = (V/V0)^(1/3).
inline void eos_vinet_pressure_and_bulk_modulus(
        double eta,
        double K0,
        double K0_prime,
        double& pressure,
        double& bulk_modulus) noexcept {
    const double inv_cbrt_eta   = 1.0 / std::cbrt(eta);
    const double exponent_coeff = 1.5 * (K0_prime - 1.0);
    const double exponential    = c_safe_exp(exponent_coeff * (1.0 - inv_cbrt_eta));
    const double inv_square     = 1.0 / (inv_cbrt_eta * inv_cbrt_eta);
    pressure     = 3.0 * K0 * (1.0 - inv_cbrt_eta) * inv_square * exponential;
    bulk_modulus = K0 * exponential
        * (2.0 - inv_cbrt_eta + exponent_coeff * inv_cbrt_eta * (1.0 - inv_cbrt_eta)) * inv_square;
}

inline double eos_bm_pressure(double eta, double K0, double K0_prime) noexcept {
    double pressure;
    double bulk_modulus;
    eos_bm_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus);
    return pressure;
}

inline double eos_vinet_pressure(double eta, double K0, double K0_prime) noexcept {
    double pressure;
    double bulk_modulus;
    eos_vinet_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus);
    return pressure;
}

inline double eos_bm_bulk_modulus(double eta, double K0, double K0_prime) noexcept {
    double pressure;
    double bulk_modulus;
    eos_bm_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus);
    return bulk_modulus;
}

inline double eos_vinet_bulk_modulus(double eta, double K0, double K0_prime) noexcept {
    double pressure;
    double bulk_modulus;
    eos_vinet_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus);
    return bulk_modulus;
}

// The compressions over which a pressure law rises, and the pressures at the two ends. Every law turns over
// in tension, and the 3rd-order Birch-Murnaghan factor 1 + (3/4)(K0'-4)(eta^(2/3)-1) changes sign at large
// eta when K0' < 4, so it turns over in compression too. An end the search never reaches stays unbounded.
// The range depends on the law's constants alone, so a model finds it once rather than per inversion.
struct c_PressureLawRange {
    double compression_min = 0.0;
    double compression_max = TidalPyConstants::d_INF;
    double pressure_min    = -TidalPyConstants::d_INF;
    double pressure_max    = TidalPyConstants::d_INF;
};

// Step outward from eta = 1 until K = eta dP/deta stops being positive, then bisect that sign change.
template <typename LawFn>
inline c_PressureLawRange eos_find_monotonic_range(double K0, double K0_prime, LawFn law_fn, double rtol) noexcept {
    c_PressureLawRange range;
    double pressure = 0.0;
    double bulk     = 0.0;
    const auto rising = [&](double eta) {
        law_fn(eta, K0, K0_prime, pressure, bulk);
        return (bulk > 0.0) && std::isfinite(pressure);
    };
    for (int side = 0; side < 2; ++side) {
        const double growth = (side == 0) ? 0.8 : 1.25;
        double inside  = 1.0;
        double outside = 1.0;
        bool   bounded = false;
        for (int k = 0; k < 200; ++k) {
            const double candidate = inside * growth;
            if (!rising(candidate)) { outside = candidate; bounded = true; break; }
            inside = candidate;
        }
        if (!bounded) { continue; }
        for (int k = 0; (k < 200) && (std::abs(outside - inside) > rtol * inside); ++k) {
            const double middle = 0.5 * (inside + outside);
            if (rising(middle)) { inside = middle; } else { outside = middle; }
        }
        law_fn(inside, K0, K0_prime, pressure, bulk);
        if (side == 0) {
            range.compression_min = inside;
            range.pressure_min    = pressure;
        } else {
            range.compression_max = inside;
            range.pressure_max    = pressure;
        }
    }
    return range;
}

// Invert a pressure law for the compression eta at a target pressure. A target past either end of the
// monotonic range has no compression to find and takes that end, keeping the answer continuous in the
// pressure; the structure solve depends on that while its central pressure is still a guess and its outer
// radii sit in tension. Inside the range this is Newton's method on the exact slope K/eta, started from
// the Murnaghan law (closed-form invertible and close to both laws over planetary compressions) and kept
// inside a bracket that every evaluation tightens; a step that leaves the bracket takes its midpoint.
template <typename LawFn>
inline double eos_invert_eta(
        double pressure_target,
        double K0,
        double K0_prime,
        LawFn law_fn,
        const c_PressureLawRange& range,
        double rtol,
        int max_iters) noexcept {
    if (!std::isfinite(pressure_target)) { return TidalPyConstants::d_NAN; }
    if (std::abs(pressure_target) <= TidalPyConstants::d_EPS) { return 1.0; }
    if (pressure_target <= range.pressure_min) { return range.compression_min; }
    if (pressure_target >= range.pressure_max) { return range.compression_max; }

    double lo = range.compression_min;
    double hi = range.compression_max;

    // Murnaghan: P = (K0 / K0') (eta^K0' - 1).
    const double murnaghan_base = 1.0 + K0_prime * pressure_target / K0;
    double eta = (murnaghan_base > 0.0 && K0_prime > 0.0) ? c_safe_pow(murnaghan_base, 1.0 / K0_prime) : 1.0;
    if (!(eta > lo && eta < hi)) { eta = 1.0; }

    double pressure = 0.0;
    double bulk     = 0.0;
    for (int i = 0; i < max_iters; ++i) {
        law_fn(eta, K0, K0_prime, pressure, bulk);
        if (pressure < pressure_target) { lo = eta; } else { hi = eta; }

        // A step smaller than the spacing of doubles lands on the bracket's edge, which is convergence
        // rather than an escape, so the bounds are inclusive.
        double next = eta + (pressure_target - pressure) * eta / bulk;
        if (!(next >= lo && next <= hi)) { next = 0.5 * (lo + hi); }

        if (std::abs(next - eta) <= rtol * eta) { return next; }
        eta = next;
    }
    return eta;  // cap reached without full convergence; best estimate
}

}  // namespace tidalpy
