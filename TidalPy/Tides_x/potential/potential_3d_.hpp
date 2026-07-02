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

// Complex-phasor variant of the engine, for the SECULAR (cycle/orbit-averaged) 3D heating. the angular factor is a
// complex amplitude (the mode's time factor e^{i omega t} pulled out): U(t) = Re[U_c e^{i omega t}].
//   even parity ((l-m) even): U = A P cos(omega t + m phi)  -> U_c = A P e^{i m phi}
//   odd parity  ((l-m) odd):  U = A P sin(omega t + m phi)  -> U_c = -i A P e^{i m phi}
// The phi derivatives bring a factor i*m (d/dphi of e^{i m phi}); the theta derivatives act on P. With
// these complex amplitudes the cycle-average heating is exact: h_bar = (omega/2) Im(sigma_c : conj(eps_c)),
// with no time grid and the U-vs-phi-derivative phase quadrature handled automatically. NOTE the e^{i m
// phi} factor cancels in Im(sigma_c conj(eps_c)) for a single mode, so the per-mode secular heating is
// longitude-independent; longitude is still honored here for generality.
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
    error_code[0] = 0;
    std::vector<c_TidalPotential3DMode> modes;

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
        if (error_code[0] != 0) { return modes; }

        EccentricityFuncOutput eccentricity_funcs =
            c_eccentricity_func(error_code, eccentricity, degree_l, eccentricity_truncation);
        if (error_code[0] != 0) { return modes; }

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
                if (!found) { error_code[0] = -20; return modes; }
            }

            const c_LegendreValue legendre = c_legendre(degree_l, order_m, colatitude);
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

                const double amplitude = G_lpq * lmp_coeff;
                const double m_d = static_cast<double>(order_m);

                // Complex phasor: base (1 for cos / -i for sin) * e^{i m phi} * amplitude.
                const std::complex<double> base = (parity == 0)
                    ? std::complex<double>(1.0, 0.0)
                    : std::complex<double>(0.0, -1.0);
                const std::complex<double> e_imphi(std::cos(m_d * longitude), std::sin(m_d * longitude));
                const std::complex<double> phasor = amplitude * base * e_imphi;
                const std::complex<double> im(0.0, m_d);   // d/dphi -> factor i*m

                c_TidalPotential3DMode out;
                out.degree_l = degree_l;
                out.mode_frequency = mode;
                out.potential = c_PotentialPointC {
                    phasor * legendre.p,                     // U_c
                    phasor * legendre.dp_dtheta,             // dU/dtheta
                    phasor * (im * legendre.p),              // dU/dphi
                    phasor * legendre.d2p_dtheta2,           // d2U/dtheta2
                    phasor * (-m_d * m_d * legendre.p),      // d2U/dphi2
                    phasor * (im * legendre.dp_dtheta)       // d2U/dtheta_dphi
                };
                modes.push_back(out);
            }
        }
    }
    return modes;
}

} // namespace tidalpy
