#pragma once
/*
 * potential_3d_.hpp - dynamic Kaula tidal-potential engine for the 3D stress, strain, and heating path.
 *
 * Given an orbital and spin state and truncation levels (max degree l, eccentricity truncation, obliquity
 * truncation), this returns every active tidal mode (l, m, p, q) with its signed forcing frequency and the
 * potential angular factor U(theta, phi, t) together with its first and second colatitude and longitude
 * derivatives (a complex c_PotentialPointC phasor). The active modes and their coefficients are built from
 * the same eccentricity (G_lpq) and obliquity (F_lmp) functions the global (1D) path uses (see
 * potential/global_.hpp), plus the associated Legendre functions P_lm (Utilities_x/legendre).
 *
 * Governing equation: Kaula's tide-raising potential, Efroimsky and Williams (2009) Eq. 18 (= Kaula 1964):
 *
 *   W(R, theta, phi, t) =  (G M_host / a) sum_l (R/a)^l
 *                            sum_m (l-m)!/(l+m)! (2 - d_m0) P_lm(cos theta)
 *                              sum_p F_lmp(i) sum_q G_lpq(e) * Trig_lm( omega_lmpq t - m phi )
 *
 * with Trig_lm = cos for (l - m) even, sin for (l - m) odd, and the tidal mode
 *   omega_lmpq = (l - 2p + q) n - m * spin        (n = orbital mean motion, spin = rotation rate),
 * dropping periapse and node precession. phi is the body-fixed east longitude (increasing in the direction of
 * rotation), zero at the host's ascending node at t = 0, where the host sits when the mean anomaly is zero. Kaula's
 * P_lm carries no Condon-Shortley phase while c_legendre does, so each amplitude carries (-1)^m to cancel it.
 * This is linear in F_lmp, G_lpq, and P_lm; the global path squares F
 * and G because global heating goes as the potential squared, while here the heating bilinearity is applied
 * downstream, after the mode stress and strain tensors are summed. The r^2 (R/a)^l coefficient is taken at
 * the surface radius R; the depth dependence is carried by the radial-solver y-functions in the kernel.
 *
 * All quantities MKS; frequencies rad s-1; angles radians. The (l-m)!/(l+m)!(2-d_m0) factor comes from
 * c_get_lm_coeff_map() (potential_common_.hpp).
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

#include "obliquity_driver_.hpp"        // c_obliquity_func, ObliquityFuncOutput
#include "eccentricity_driver_.hpp"     // c_eccentricity_func, EccentricityFuncOutput
#include "potential_common_.hpp"        // c_get_lm_coeff_map, keys
#include "potential_point_.hpp"         // tidalpy::c_PotentialPointC
#include "legendre_driver_.hpp"         // tidalpy::c_legendre
#include "constants_.hpp"               // TidalPyConstants

namespace tidalpy {

// One active tidal mode: degree, signed forcing frequency, and complex potential angular-factor amplitude,
// the time factor pulled out so U(t) = Re[U_c e^{i omega t}]. The 3D stress, strain, and heating paths
// consume the complex amplitude, which keeps the cycle average exact.
struct c_TidalPotential3DMode {
    int degree_l = 0;
    double mode_frequency = 0.0;      // signed omega_lmpq [rad s-1]
    c_PotentialPointC potential;      // complex amplitude of U and its theta/phi derivatives
};

// Everything about a mode the engine can determine without a colatitude or longitude; the angular factor at
// a point is evaluated on demand by c_eval_potential_point_3d. Splitting the engine this way lets a batch
// or map path build the mode list once and amortize the radial solve across every query point.
struct c_TidalPotential3DModeCoeff {
    int degree_l = 0;
    int order_m = 0;
    int p = 0;
    int q = 0;
    int parity = 0;                   // (l - m) & 1: 0 -> cos (even), 1 -> sin (odd)
    double mode_frequency = 0.0;      // signed omega_lmpq [rad s-1]
    double amplitude = 0.0;           // (-1)^m G_lpq F_lmp (R/a)^l (G M_host / a) (l-m)!/(l+m)! (2-d_m0)
    double amplitude_factor = 0.0;    // the amplitude without G_lpq, so products of G can be cut at e^N
    c_EccentricitySeriesTable eccentricity_series;   // the degree's eccentricity table, for those products
};

// Evaluate a mode's complex potential angular factor U_c and its theta and phi derivatives at a point:
//   even parity ((l-m) even): U_c = amplitude * P_lm(cos theta) * e^{-i m phi}
//   odd parity  ((l-m) odd):  U_c = -i * amplitude * P_lm(cos theta) * e^{-i m phi}
// Theta derivatives act on P_lm; phi derivatives bring a factor -i*m from d/dphi of e^{-i m phi}.
inline c_PotentialPointC c_eval_potential_point_3d(
        const c_TidalPotential3DModeCoeff& coeff,
        double colatitude,
        double longitude)
{
    const c_LegendreValue legendre = c_legendre(coeff.degree_l, coeff.order_m, colatitude);
    const double m_d = static_cast<double>(coeff.order_m);

    const std::complex<double> base = (coeff.parity == 0)
        ? std::complex<double>(1.0, 0.0)
        : std::complex<double>(0.0, -1.0);
    const std::complex<double> e_imphi(std::cos(m_d * longitude), -std::sin(m_d * longitude));
    const std::complex<double> phasor = coeff.amplitude * base * e_imphi;
    const std::complex<double> im(0.0, -m_d);   // d/dphi -> factor -i*m

    return c_PotentialPointC {
        phasor * legendre.p,                     // U_c
        phasor * legendre.dp_dtheta,             // dU/dtheta
        phasor * (im * legendre.p),              // dU/dphi
        phasor * legendre.d2p_dtheta2,           // d2U/dtheta2
        phasor * (-m_d * m_d * legendre.p),      // d2U/dphi2
        phasor * (im * legendre.dp_dtheta)       // d2U/dtheta_dphi
    };
}

// The same mode-discovery loop as the global (1D) engine, but linear in F_lmp and G_lpq, since the
// amplitude carries one power of each, and carrying the (l, m) labels so P_lm can be evaluated per point
// later. Nothing here depends on position, so a batch or map path builds it once.
inline std::vector<c_TidalPotential3DModeCoeff> c_tidal_potential_3d_mode_coeffs(
        double planet_radius,
        double semi_major_axis,
        double orbital_frequency,
        double spin_frequency,
        double obliquity,
        double eccentricity,
        double host_mass,
        double G_to_use,
        int min_degree_l,
        int max_degree_l,
        int obliquity_truncation,
        int eccentricity_truncation,
        int* error_code)
{
    error_code[0] = 0;
    std::vector<c_TidalPotential3DModeCoeff> coeffs;

    const double R_a  = planet_radius / semi_major_axis;
    double ra_l_coeff = std::pow(R_a, static_cast<double>(min_degree_l)) * (G_to_use * host_mass / semi_major_axis);

    auto& lm_coeff_map = c_get_lm_coeff_map();
    c_Key2 lm_key;

    for (int degree_l = min_degree_l; degree_l <= max_degree_l; ++degree_l)
    {
        if (degree_l > min_degree_l)
        {
            ra_l_coeff *= R_a;
        }

        ObliquityFuncOutput obliquity_funcs =
            c_obliquity_func(error_code, obliquity, degree_l, obliquity_truncation);
        if (error_code[0] != 0) { return coeffs; }

        EccentricityFuncOutput eccentricity_funcs =
            c_eccentricity_func(error_code, eccentricity, degree_l, eccentricity_truncation);
        if (error_code[0] != 0) { return coeffs; }
        const c_EccentricitySeriesTable eccentricity_series =
            c_eccentricity_series_table(error_code, degree_l, eccentricity_truncation);
        if (error_code[0] != 0) { return coeffs; }

        lm_key.a = degree_l;
        lm_key.b = -1;
        double lm_coeff = TidalPyConstants::d_NAN;

        for (const auto& [lmp_key, F_lmp] : obliquity_funcs.first)
        {
            if (F_lmp == 0.0) { continue; }

            const int order_m = lmp_key.b;
            if (order_m != lm_key.b)
            {
                lm_key.b = order_m;
                lm_key.rebuild_reference();
                bool found = false;
                lm_coeff = lm_coeff_map.get(found, lm_key);
                if (!found) { error_code[0] = -20; return coeffs; }
            }

            const int parity = (degree_l - order_m) & 1;
            const double lmp_coeff = F_lmp * ra_l_coeff * lm_coeff;

            bool found = false;
            const c_IntMap<c_Key1, double>* ecc_by_q =
                eccentricity_funcs.second.get_ptr(found, c_Key2(lmp_key.a, lmp_key.c));
            if (!found) { continue; }

            for (const auto& [q_key, G_lpq] : *ecc_by_q)
            {
                if (G_lpq == 0.0) { continue; }

                const int q = q_key.a;
                const double mode =
                    static_cast<double>(degree_l - 2 * lmp_key.c + q) * orbital_frequency
                    - static_cast<double>(order_m) * spin_frequency;

                c_TidalPotential3DModeCoeff out;
                out.degree_l = degree_l;
                out.order_m = order_m;
                out.p = lmp_key.c;
                out.q = q;
                out.parity = parity;
                out.mode_frequency = mode;
                out.amplitude_factor = ((order_m & 1) ? -1.0 : 1.0) * lmp_coeff;
                out.amplitude = out.amplitude_factor * G_lpq;
                out.eccentricity_series = eccentricity_series;
                coeffs.push_back(out);
            }
        }
    }
    return coeffs;
}

// Per-mode complex-phasor angular factors at a single point, the raw per-(l, m, p, q) view. The heating
// paths do not consume it directly: they first merge the modes into coherent waves, because modes sharing a
// real spatial function must be summed before the heating bilinear form.
inline std::vector<c_TidalPotential3DMode> c_tidal_potential_3d_modes(
        double planet_radius,
        double semi_major_axis,
        double orbital_frequency,
        double spin_frequency,
        double obliquity,
        double eccentricity,
        double host_mass,
        double G_to_use,
        int min_degree_l,
        int max_degree_l,
        int obliquity_truncation,
        int eccentricity_truncation,
        double colatitude,
        double longitude,
        int* error_code)
{
    const std::vector<c_TidalPotential3DModeCoeff> coeffs = c_tidal_potential_3d_mode_coeffs(
        planet_radius,
        semi_major_axis,
        orbital_frequency,
        spin_frequency,
        obliquity,
        eccentricity,
        host_mass,
        G_to_use,
        min_degree_l,
        max_degree_l,
        obliquity_truncation,
        eccentricity_truncation,
        error_code);

    std::vector<c_TidalPotential3DMode> modes;
    if (error_code[0] != 0) { return modes; }
    modes.reserve(coeffs.size());

    for (const c_TidalPotential3DModeCoeff& coeff : coeffs)
    {
        c_TidalPotential3DMode out;
        out.degree_l = coeff.degree_l;
        out.mode_frequency = coeff.mode_frequency;
        out.potential = c_eval_potential_point_3d(coeff, colatitude, longitude);
        modes.push_back(out);
    }
    return modes;
}

// ---------------------------------------------------------------------------------------------------------------
// Coherent waves
// ---------------------------------------------------------------------------------------------------------------
// A coherent tidal wave for the 3D heating paths. Every active mode is mapped onto a non-negative frequency,
// a mode at omega < 0 contributing the conjugate of its phasor at +|omega| since Re[U_c e^{i omega t}] =
// Re[conj(U_c) e^{-i omega t}], and merged with every other mode sharing the same real spatial function:
// same l, m, |omega|, and azimuthal sign. The complex amplitude carries the parity phase (-i for odd l - m),
// the conjugation, and the coherent sum over the merged modes.
//
// The merge matters because the secular heating at a frequency is (|omega|/2) Im(sigma_c : conj(eps_c)) of
// the total complex amplitude there. Summing the cycle-averaged powers of two modes that are the same real
// sinusoid halves their contribution: (a + a)^2 = 4 a^2, not 2 a^2. The m = 0 modes always come in such
// pairs, (l, 0, p, q) at +omega and (l, 0, l-p, -q) at -omega carrying equal amplitudes and the same
// function of time. For a homogeneous body at zero obliquity the zonal terms are 9/84 of the degree-2
// heating, so summing them incoherently loses 5.36% of the total at synchronous rotation. At nonzero
// obliquity, modes of one (l, m) with different (p, q) can share a signed frequency too; their relative
// phase is set by the argument of periapse, which this engine takes as zero, so they also merge coherently.
//
// A wave also keeps its member modes, each with its amplitude less G_lpq, because the secular heating is quadratic
// in the amplitudes: products of two eccentricity functions are cut at e^N (c_wave_pair_power_3d), as in the global
// (1D) path, while the instantaneous fields use the unsquared amplitude.
struct c_WaveMember3D {
    c_EccentricitySeriesTable eccentricity_series;
    int p = 0;
    int q = 0;
    std::complex<double> factor {0.0, 0.0};   // amplitude without G_lpq (parity phase and conjugation applied)
};

struct c_TidalWave3D {
    int degree_l = 0;
    int order_m = 0;
    int azimuthal_sign = 1;                 // -1: e^{-i m phi}; +1: e^{+i m phi} (conjugated omega < 0 mode, or m = 0)
    double frequency = 0.0;                 // |omega| [rad s-1], > 0
    std::complex<double> amplitude {0.0, 0.0};   // coherent complex amplitude (parity phase, conjugation, merge applied)
    std::vector<c_WaveMember3D> members;
};

// The product amplitude_a * conj(amplitude_b) of two waves with every product of two eccentricity functions cut at
// e^N, which is what the secular heating takes in place of the plain product.
inline std::complex<double> c_wave_pair_power_3d(
        const c_TidalWave3D& wave_a,
        const c_TidalWave3D& wave_b,
        double eccentricity) noexcept
{
    std::complex<double> total(0.0, 0.0);
    for (const c_WaveMember3D& member_a : wave_a.members)
    {
        for (const c_WaveMember3D& member_b : wave_b.members)
        {
            const double product = c_eccentricity_cut_product(
                member_a.eccentricity_series,
                member_a.p,
                member_a.q,
                member_b.eccentricity_series,
                member_b.p,
                member_b.q,
                eccentricity);
            if (product != 0.0)
            {
                total += member_a.factor * std::conj(member_b.factor) * product;
            }
        }
    }
    return total;
}

// Integer combinations of n and the spin rate that agree mathematically can still differ at the last bit; the
// tolerance is the 1D path's (c_frequency_match_rtol, [numerical] frequency_match_rtol).
inline bool c_tidal_wave_same_frequency(double frequency_a, double frequency_b) {
    return std::abs(frequency_a - frequency_b)
        <= c_frequency_match_rtol() * std::max(std::abs(frequency_a), std::abs(frequency_b));
}

// Drops modes at or below min_frequency, which dissipate nothing, and waves whose merged amplitude cancels.
inline std::vector<c_TidalWave3D> c_coherent_tidal_waves_3d(
        const std::vector<c_TidalPotential3DModeCoeff>& modes,
        double min_frequency)
{
    std::vector<c_TidalWave3D> waves;
    waves.reserve(modes.size());
    for (const c_TidalPotential3DModeCoeff& mode : modes)
    {
        const double frequency = std::abs(mode.mode_frequency);
        if (frequency <= min_frequency) {
            continue;
        }

        const bool negative = mode.mode_frequency < 0.0;
        // Parity phase of the phasor: 1 for cos (even l - m), -i for sin (odd l - m); conjugated for omega < 0.
        const std::complex<double> phase = (mode.parity == 0)
            ? std::complex<double>(1.0, 0.0)
            : std::complex<double>(0.0, negative ? 1.0 : -1.0);
        const int azimuthal_sign = (mode.order_m == 0 || negative) ? 1 : -1;
        const std::complex<double> amplitude = mode.amplitude * phase;
        c_WaveMember3D member;
        member.eccentricity_series = mode.eccentricity_series;
        member.p = mode.p;
        member.q = mode.q;
        member.factor = mode.amplitude_factor * phase;

        bool merged = false;
        for (c_TidalWave3D& wave : waves)
        {
            if (wave.degree_l == mode.degree_l && wave.order_m == mode.order_m
                && wave.azimuthal_sign == azimuthal_sign
                && c_tidal_wave_same_frequency(wave.frequency, frequency))
            {
                wave.amplitude += amplitude;
                wave.members.push_back(member);
                merged = true;
                break;
            }
        }
        if (!merged)
        {
            c_TidalWave3D wave;
            wave.degree_l       = mode.degree_l;
            wave.order_m        = mode.order_m;
            wave.azimuthal_sign = azimuthal_sign;
            wave.frequency      = frequency;
            wave.amplitude      = amplitude;
            wave.members.push_back(member);
            waves.push_back(wave);
        }
    }
    std::vector<c_TidalWave3D> active;
    active.reserve(waves.size());
    for (const c_TidalWave3D& wave : waves)
    {
        if (std::abs(wave.amplitude) > 0.0) { active.push_back(wave); }
    }
    return active;
}

// e^{i mu phi}, the longitude phasor of a wave with signed azimuthal wavenumber mu.
inline std::complex<double> c_azimuthal_phasor(double mu, double longitude)
{
    return std::complex<double>(std::cos(mu * longitude), std::sin(mu * longitude));
}

// A coherent wave's U_c and its theta/phi derivatives from its Legendre values and phasor, for the given amplitude (the
// wave's own, or 1 for a per-unit-amplitude field):
//   U_c = amplitude * P_lm(cos theta) * e^{i mu phi},   mu = azimuthal_sign * m,
// with U(t) = Re[U_c e^{i |omega| t}]. Theta derivatives act on P_lm, phi derivatives bring a factor i*mu.
// A grid shares the Legendre values across the waves of one (l, m) and the phasor across those of one mu.
inline c_PotentialPointC c_wave_point_from_parts(
        const c_TidalWave3D& wave,
        const c_LegendreValue& legendre,
        const std::complex<double>& e_imuphi,
        const std::complex<double>& amplitude)
{
    const double mu = static_cast<double>(wave.azimuthal_sign) * static_cast<double>(wave.order_m);
    const std::complex<double> phasor = amplitude * e_imuphi;
    const std::complex<double> i_mu(0.0, mu);

    return c_PotentialPointC {
        phasor * legendre.p,                     // U_c
        phasor * legendre.dp_dtheta,             // dU/dtheta
        phasor * (i_mu * legendre.p),            // dU/dphi
        phasor * legendre.d2p_dtheta2,           // d2U/dtheta2
        phasor * (-mu * mu * legendre.p),        // d2U/dphi2
        phasor * (i_mu * legendre.dp_dtheta)     // d2U/dtheta_dphi
    };
}

// Evaluate a coherent wave's complex potential angular factor U_c and its theta/phi derivatives at a point.
inline c_PotentialPointC c_eval_wave_point_3d(
        const c_TidalWave3D& wave,
        double colatitude,
        double longitude)
{
    const double mu = static_cast<double>(wave.azimuthal_sign) * static_cast<double>(wave.order_m);
    return c_wave_point_from_parts(
        wave,
        c_legendre(wave.degree_l, wave.order_m, colatitude),
        c_azimuthal_phasor(mu, longitude),
        wave.amplitude);
}

} // namespace tidalpy
