// kernel_.hpp - 3D tidal strain/stress/heating kernel.
//
// Combines the radial strain coefficients (strain_radial_.hpp) with the tidal potential's angular
// factor and its theta/phi derivatives, at a single point, to produce the 6 complex strain and
// stress tensor components and the volumetric heating.
//
// References: TB05 Eq. 10 / B13 Eqs. 8-9 with the Kervazo+2021 (A&A App. D) correction to the
// Tobie+2005 theta-phi / phi-phi typo, Takeuchi & Saito isotropic constitutive law, Europa-book
// Eq. 42 heating.
//
// Conventions: the tidal potential U and its derivatives are real (the potential
// is cos/sin in omega*t); the radial coefficients and moduli are complex (viscoelastic); strains and
// stresses are therefore complex. The 6 tensor components are ordered
//   [0]=rr, [1]=theta-theta, [2]=phi-phi, [3]=r-theta, [4]=r-phi, [5]=theta-phi.
#pragma once

#include <array>
#include <cmath>
#include <complex>

#include "strain_radial_.hpp"
#include "../potential/potential_point_.hpp"   // tidalpy::c_PotentialPoint (shared with the potential models)


namespace tidalpy {
namespace tides {

// The 6 complex strain (or stress) tensor components at a point.
struct c_Tensor6
{
    std::array<std::complex<double>, 6> c { };
};

// Compute the 6 strain and 6 stress components at one point from the radial coefficients, the
// potential point, and the colatitude (needed for sin/cot factors).
// If the radial coefficients are invalid (liquid / center), the components are NaN.
// Templated on the potential-point type so it serves both the real c_PotentialPoint (instantaneous
// path, U real) and the complex c_PotentialPointC (secular path, U is a complex phasor amplitude).
template <typename PotentialPointT>
inline void c_compute_strain_stress(
        const c_StrainRadialCoeffs& R,
        const PotentialPointT& P,
        double colatitude,
        c_Tensor6& strain_out,
        c_Tensor6& stress_out) noexcept
{
    if (!R.valid)
    {
        const std::complex<double> c_nan(
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
        );

        for (size_t k = 0; k < 6; ++k)
        {
            strain_out.c[k] = c_nan;
            stress_out.c[k] = c_nan;
        }
        return;
    }

    const double sin_theta = std::sin(colatitude);
    const double sin_inv   = (sin_theta == 0.0) ? std::numeric_limits<double>::quiet_NaN() : 1.0 / sin_theta;
    const double tan_theta = std::tan(colatitude);
    const double cot_theta = (tan_theta == 0.0) ? std::numeric_limits<double>::quiet_NaN() : 1.0 / tan_theta;

    // Angular helper combinations (auto: real for c_PotentialPoint, complex for c_PotentialPointC).
    const auto s2_t1 = (sin_inv * sin_inv) * P.d2U_dphi2 + cot_theta * P.dU_dtheta;
    const auto s4_t0 = P.dU_dphi * sin_inv;
    const auto s5_t0 = 2.0 * (P.d2U_dtheta_dphi - cot_theta * P.dU_dphi) * sin_inv;

    const std::complex<double> y1_r_U = R.y1_over_r * P.U;

    // Strain components.
    strain_out.c[0] = R.dy1_dr * P.U;                        // eps_rr
    strain_out.c[1] = R.y3_over_r * P.d2U_dtheta2 + y1_r_U;  // eps_thth
    strain_out.c[2] = y1_r_U + R.y3_over_r * s2_t1;          // eps_phph  (Kervazo-corrected)
    strain_out.c[3] = R.y4_over_2mu * P.dU_dtheta;           // eps_rth
    strain_out.c[4] = R.y4_over_2mu * s4_t0;                 // eps_rph
    strain_out.c[5] = R.y3_over_2r * s5_t0;                  // eps_thph  (Kervazo-corrected)

    // Stress (isotropic, Takeuchi & Saito): sigma = 2 mu eps + lame tr(eps) delta.
    const std::complex<double> trace = strain_out.c[0] + strain_out.c[1] + strain_out.c[2];
    const std::complex<double> trace_lame = R.lame * trace;
    const std::complex<double> two_mu = 2.0 * R.shear;
    for (size_t k = 0; k < 6; ++k)
    {
        stress_out.c[k] = two_mu * strain_out.c[k];
        if (k < 3) stress_out.c[k] += trace_lame;
    }
}

// Volumetric tidal heating [W m-3] at a point from the 6 stress and strain components.
// h = | sum_k [ Im(sigma_k) Re(eps_k) - Re(sigma_k) Im(eps_k) ] |, with factor 2 on the 3 off-diagonals.
inline double c_volumetric_heating(const c_Tensor6& stress, const c_Tensor6& strain) noexcept
{
    double h = 0.0;
    for (size_t k = 0; k < 6; ++k)
    {
        const double term = stress.c[k].imag() * strain.c[k].real()
                          - stress.c[k].real() * strain.c[k].imag();
        h += (k < 3) ? term : 2.0 * term;
    }
    // h can be complex. According to the Europa book eq. 42 we can take the abs to find the true volumetric heating.
    // TODO: this feels hacky though, would love to investigate this further.
    return std::abs(h);
}

// Signed cycle-average heating factor at a point: sum_k w_k [ Im(sigma_k) Re(eps_k) - Re(sigma_k)
// Im(eps_k) ] = sum_k w_k Im(sigma_k conj(eps_k)), with factor 2 on the three off-diagonals and NO
// abs(). For the secular (cycle/orbit-averaged) heating this is called with the complex-phasor
// stress/strain of a single mode; the per-mode volumetric heating is (omega_mode/2) times this, and
// modes sum with sign (distinct-frequency cross terms average to zero, so they are simply omitted).
inline double c_volumetric_heating_signed(const c_Tensor6& stress, const c_Tensor6& strain) noexcept
{
    double h = 0.0;
    for (size_t k = 0; k < 6; ++k)
    {
        const double term = stress.c[k].imag() * strain.c[k].real()
                          - stress.c[k].real() * strain.c[k].imag();
        h += (k < 3) ? term : 2.0 * term;
    }
    return h;
}

}  // namespace tides
}  // namespace tidalpy
