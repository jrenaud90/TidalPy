#pragma once
/*
 * potential_3d_.hpp - Dynamic Kaula tidal-potential engine for the 3D stress/strain/heating path.
 *
 * Given an orbital/spin state and truncation levels (max degree l, eccentricity truncation, obliquity
 * truncation), this returns every active tidal mode (l, m, p, q) with its signed forcing frequency and
 * the potential angular factor U(theta, phi, t) together with its first/second colatitude/longitude
 * derivatives (a c_PotentialPoint). It is the class-free replacement for the old per-scenario
 * c_TidalPotentialBase models: the active modes and their coefficients are built dynamically from the
 * same eccentricity (G_lpq) and inclination/obliquity (F_lmp) functions the global (1D) path uses
 * (see potential/global_.hpp), plus the associated Legendre functions P_lm (Utilities_x/legendre).
 *
 * Governing equation: Kaula's tide-raising potential, Efroimsky & Williams (2009) Eq. 18 (= Kaula 1964):
 *
 *   W(R, theta, phi, t) = -(G M_host / a) sum_l (R/a)^l
 *                            sum_m (l-m)!/(l+m)! (2 - d_m0) P_lm(cos theta)
 *                              sum_p F_lmp(i) sum_q G_lpq(e) * Trig_lm( omega_lmpq t + m phi )
 *
 * with Trig_lm = cos for (l - m) even, sin for (l - m) odd, and the tidal mode
 *   omega_lmpq = (l - 2p + q) n - m * spin        (n = orbital mean motion, spin = rotation rate),
 * dropping periapse/node precession. This is LINEAR in F_lmp, G_lpq, and P_lm (the global path squares
 * F and G because global heating goes as the potential squared; here the heating bilinearity is applied
 * downstream, after the mode stress/strain tensors are summed). The r^2 (R/a)^l coefficient is taken at
 * the surface radius R; the depth dependence is carried by the radial-solver y-functions in the kernel.
 *
 * All quantities MKS; frequencies rad s-1; angles radians. The (l-m)!/(l+m)!(2-d_m0) factor is provided
 * by c_get_lm_coeff_map() (potential_common_.hpp).
 */

#include <cmath>
#include <cstdint>
#include <vector>

#include "obliquity_driver_.hpp"        // c_obliquity_func, ObliquityFuncOutput
#include "eccentricity_driver_.hpp"     // c_eccentricity_func, EccentricityFuncOutput
#include "potential_common_.hpp"        // c_get_lm_coeff_map, keys
#include "potential_point_.hpp"         // tidalpy::c_PotentialPoint
#include "legendre_driver_.hpp"         // tidalpy::c_legendre
#include "constants_.hpp"               // TidalPyConstants

namespace tidalpy {

// One active tidal mode: its degree, signed forcing frequency, and COMPLEX potential angular-factor
// amplitude (the mode's time factor e^{i omega t} pulled out, U(t) = Re[U_c e^{i omega t}]). The 3D
// stress/strain/heating consumes these complex amplitudes so the cycle-average is exact (no time grid).
struct c_TidalPotential3DMode {
    int degree_l = 0;
    double mode_frequency = 0.0;      // signed omega_lmpq [rad s-1]
    c_PotentialPointC potential;      // complex amplitude of U and its theta/phi derivatives
};

// Position-INDEPENDENT description of one active tidal mode: everything the engine can determine
// without a colatitude/longitude. The angular factor at a point (P_lm and the e^{i m phi} phasor) is
// evaluated on demand by c_eval_potential_point_3d below. Splitting the engine this way lets a batch /
// map path build the mode list once and reuse it across every (radius, colatitude) query, and lets the
// radial (Love-number) solve be amortized across points, since neither depends on m or on position.
struct c_TidalPotential3DModeCoeff {
    int degree_l = 0;
    int order_m = 0;
    int parity = 0;                   // (l - m) & 1: 0 -> cos (even), 1 -> sin (odd)
    double mode_frequency = 0.0;      // signed omega_lmpq [rad s-1]
    double amplitude = 0.0;           // G_lpq * F_lmp * (R/a)^l * (G M_host / a) * (l-m)!/(l+m)!(2-d_m0)
};

// Evaluate a mode's complex potential angular factor U_c and its theta/phi derivatives at a point.
// This is the position-DEPENDENT part of the Kaula engine, factored out of c_tidal_potential_3d_modes:
//   even parity ((l-m) even): U_c = amplitude * P_lm(cos theta) * e^{i m phi}
//   odd parity  ((l-m) odd):  U_c = -i * amplitude * P_lm(cos theta) * e^{i m phi}
// theta derivatives act on P_lm; phi derivatives bring a factor i*m from d/dphi of e^{i m phi}.
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
    const std::complex<double> e_imphi(std::cos(m_d * longitude), std::sin(m_d * longitude));
    const std::complex<double> phasor = coeff.amplitude * base * e_imphi;
    const std::complex<double> im(0.0, m_d);   // d/dphi -> factor i*m

    return c_PotentialPointC {
        phasor * legendre.p,                     // U_c
        phasor * legendre.dp_dtheta,             // dU/dtheta
        phasor * (im * legendre.p),              // dU/dphi
        phasor * legendre.d2p_dtheta2,           // d2U/dtheta2
        phasor * (-m_d * m_d * legendre.p),      // d2U/dphi2
        phasor * (im * legendre.dp_dtheta)       // d2U/dtheta_dphi
    };
}

// Build the position-INDEPENDENT list of active tidal modes (degree, order, parity, signed frequency,
// scalar amplitude) from the orbital/spin state and truncation levels. This is the mode-discovery loop
// shared with the global (1D) engine, but LINEAR in F_lmp, G_lpq (the amplitude carries one power of
// each) and carrying the (l, m) labels so the P_lm angular factor can be evaluated per point later.
// Neither the mode list nor the downstream radial (Love-number) solve depends on colatitude/longitude,
// so a batch / map path builds this once and reuses it across all query points.
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
                out.parity = parity;
                out.mode_frequency = mode;
                out.amplitude = G_lpq * lmp_coeff;
                coeffs.push_back(out);
            }
        }
    }
    return coeffs;
}

// Complex-phasor engine for the SECULAR (cycle/orbit-averaged) 3D heating at a single point. the angular
// factor is a complex amplitude (the mode's time factor e^{i omega t} pulled out): U(t) = Re[U_c e^{i
// omega t}]. Thin wrapper: build the position-independent mode list once, then evaluate each mode's
// angular factor at (colatitude, longitude). NOTE the e^{i m phi} factor cancels in Im(sigma_c
// conj(eps_c)) for a single mode, so the per-mode secular heating is longitude-independent; longitude is
// still honored here for generality.
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
        planet_radius, semi_major_axis, orbital_frequency, spin_frequency, obliquity, eccentricity,
        host_mass, G_to_use, min_degree_l, max_degree_l, obliquity_truncation, eccentricity_truncation,
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

} // namespace tidalpy
