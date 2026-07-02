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

// One active tidal mode's degree, signed forcing frequency, and evaluated potential angular factor.
struct c_TidalPotential3DMode {
    int degree_l = 0;
    double mode_frequency = 0.0;      // signed omega_lmpq [rad s-1]
    c_PotentialPoint potential;       // U and its theta/phi first and second derivatives at the point
};

// Evaluate every active tidal mode's potential angular factor at one point. The radius is the SURFACE
// radius (used only in the (R/a)^l coefficient); colatitude/longitude/time locate the point. Truncation
// levels select the eccentricity (G_lpq) and obliquity (F_lmp) function sets. On an unsupported degree /
// truncation *error_code is set nonzero and the partial result is returned.
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
        double time,
        int* error_code)
{
    error_code[0] = 0;
    std::vector<c_TidalPotential3DMode> modes;

    const double cos_lon = std::cos(longitude);   // reserved for future m-longitude fast paths
    (void)cos_lon;

    // (R/a) powers; multiply the outer (G M_host / a) into the l-power coefficient.
    const double R_a   = planet_radius / semi_major_axis;
    const double R_a_2 = R_a * R_a;
    double ra_l_coeff  = std::pow(R_a, static_cast<double>(min_degree_l)) * (G_to_use * host_mass / semi_major_axis);

    auto& lm_coeff_map = c_get_lm_coeff_map();
    c_Key2 lm_key;

    for (int degree_l = min_degree_l; degree_l <= max_degree_l; ++degree_l)
    {
        if (degree_l > min_degree_l)
        {
            ra_l_coeff *= R_a;   // grow (R/a)^l by one power per degree
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

        // Outer loop over the inclination-function (l, m, p) triples.
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

            // Angular (Legendre) part depends only on (l, m); evaluate once per (l, m, p) triple.
            const c_LegendreValue legendre = c_legendre(degree_l, order_m, colatitude);
            const int parity = (degree_l - order_m) & 1;   // 0 -> cos, 1 -> sin

            const double lmp_coeff = F_lmp * ra_l_coeff * lm_coeff;   // LINEAR in F (not squared)

            // Eccentricity vector G_lpq(q) for this (l, p).
            bool found = false;
            const c_IntMap<c_Key1, double>* ecc_by_q =
                eccentricity_funcs.second.get_ptr(found, c_Key2(lmp_key.a, lmp_key.c));   // (l, p)
            if (!found) { continue; }

            for (const auto& [q_key, G_lpq] : *ecc_by_q)
            {
                if (G_lpq == 0.0) { continue; }

                const int q = q_key.a;
                const double mode =
                    static_cast<double>(degree_l - 2 * lmp_key.c + q) * orbital_frequency
                    - static_cast<double>(order_m) * spin_frequency;

                const double amplitude = G_lpq * lmp_coeff;   // A_lmpq (linear in G)

                // Trig( omega_lmpq t + m phi ) with parity, and its phi derivatives (d/dphi -> factor m).
                const double arg = mode * time + static_cast<double>(order_m) * longitude;
                const double m_d = static_cast<double>(order_m);
                double trig, dtrig_dphi;
                if (parity == 0)
                {
                    trig = std::cos(arg);
                    dtrig_dphi = -m_d * std::sin(arg);
                }
                else
                {
                    trig = std::sin(arg);
                    dtrig_dphi =  m_d * std::cos(arg);
                }
                
                const double d2trig_dphi2 = -m_d * m_d * trig;

                c_TidalPotential3DMode out;
                out.degree_l = degree_l;
                out.mode_frequency = mode;
                out.potential = c_PotentialPoint {
                    amplitude * legendre.p * trig,               // U
                    amplitude * legendre.dp_dtheta * trig,       // dU/dtheta
                    amplitude * legendre.p * dtrig_dphi,         // dU/dphi
                    amplitude * legendre.d2p_dtheta2 * trig,     // d2U/dtheta2
                    amplitude * legendre.p * d2trig_dphi2,       // d2U/dphi2
                    amplitude * legendre.dp_dtheta * dtrig_dphi  // d2U/dtheta_dphi
                };
                modes.push_back(out);
            }
        }
    }
    return modes;
}

} // namespace tidalpy
