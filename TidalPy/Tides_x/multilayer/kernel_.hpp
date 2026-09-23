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
#include <limits>

#include "constants_.hpp"                        // TidalPyConstants::d_EPS
#include "strain_radial_.hpp"
#include "../potential/potential_point_.hpp"   // tidalpy::c_PotentialPointC (shared with the potential engine)


namespace tidalpy {
namespace tides {

// The 6 complex strain (or stress) tensor components at a point.
struct c_Tensor6
{
    std::array<std::complex<double>, 6> c { };
};

// 1 / sin(theta) and 1 / tan(theta) as the strain uses them, NaN at a pole. A colatitude whose sine is within machine
// epsilon of zero is a pole, so 0 and pi (whose sine is 1.2e-16, not 0) are treated alike.
struct c_ColatitudeTrig
{
    double sin_inv   = 0.0;
    double cot_theta = 0.0;
};

inline c_ColatitudeTrig c_colatitude_trig(double colatitude) noexcept
{
    const double sin_theta = std::sin(colatitude);
    c_ColatitudeTrig trig;
    if (std::abs(sin_theta) <= TidalPyConstants::d_EPS) {
        trig.sin_inv   = std::numeric_limits<double>::quiet_NaN();
        trig.cot_theta = std::numeric_limits<double>::quiet_NaN();
        return trig;
    }
    trig.sin_inv   = 1.0 / sin_theta;
    trig.cot_theta = 1.0 / std::tan(colatitude);
    return trig;
}

// The part of the strain at a point that depends on the potential point and the colatitude but not on radius. A grid
// forms these once per (colatitude, longitude) and reuses them at every radius.
struct c_AngularStrainFactors
{
    std::complex<double> U;             // U
    std::complex<double> dU_dtheta;     // dU/dtheta
    std::complex<double> d2U_dtheta2;   // d2U/dtheta2
    std::complex<double> s2_t1;         // d2U/dphi2 / sin^2(theta) + cot(theta) dU/dtheta
    std::complex<double> s4_t0;         // dU/dphi / sin(theta)
    std::complex<double> s5_t0;         // 2 (d2U/dtheta dphi - cot(theta) dU/dphi) / sin(theta)
};

// The angular strain factors of a potential point at a colatitude. Templated on the potential-point type like
// c_compute_strain_stress.
template <typename PotentialPointT>
inline c_AngularStrainFactors c_angular_strain_factors(const PotentialPointT& P, const c_ColatitudeTrig& trig) noexcept
{
    c_AngularStrainFactors A;
    A.U           = P.U;
    A.dU_dtheta   = P.dU_dtheta;
    A.d2U_dtheta2 = P.d2U_dtheta2;
    A.s2_t1       = (trig.sin_inv * trig.sin_inv) * P.d2U_dphi2 + trig.cot_theta * P.dU_dtheta;
    A.s4_t0       = P.dU_dphi * trig.sin_inv;
    A.s5_t0       = 2.0 * (P.d2U_dtheta_dphi - trig.cot_theta * P.dU_dphi) * trig.sin_inv;
    return A;
}

// Compute the 6 strain and 6 stress components at one point from the radial coefficients and the angular strain
// factors. If the radial coefficients are invalid (liquid / center), the components are NaN.
inline void c_compute_strain_stress_from_factors(
        const c_StrainRadialCoeffs& R,
        const c_AngularStrainFactors& A,
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

    const std::complex<double> y1_r_U = R.y1_over_r * A.U;

    // Strain components.
    strain_out.c[0] = R.dy1_dr * A.U;                        // eps_rr
    strain_out.c[1] = R.y3_over_r * A.d2U_dtheta2 + y1_r_U;  // eps_thth
    strain_out.c[2] = y1_r_U + R.y3_over_r * A.s2_t1;        // eps_phph  (Kervazo-corrected)
    strain_out.c[3] = R.y4_over_2mu * A.dU_dtheta;           // eps_rth
    strain_out.c[4] = R.y4_over_2mu * A.s4_t0;               // eps_rph
    strain_out.c[5] = R.y3_over_2r * A.s5_t0;                // eps_thph  (Kervazo-corrected)

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

// As above, but taking the potential point and the colatitude (needed for the sin and cot factors).
// Templated on the potential-point type; every caller passes the complex phasor c_PotentialPointC, so the
// strains and stresses are complex amplitudes at the mode's frequency.
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
        c_compute_strain_stress_from_factors(R, c_AngularStrainFactors{}, strain_out, stress_out);
        return;
    }
    c_compute_strain_stress_from_factors(
        R,
        c_angular_strain_factors(P, c_colatitude_trig(colatitude)),
        strain_out,
        stress_out);
}

// Magnitude of the weighted bilinear form [Pa] at a point from the 6 stress and strain amplitudes:
// h = | sum_k [ Im(sigma_k) Re(eps_k) - Re(sigma_k) Im(eps_k) ] |, with factor 2 on the 3 off-diagonals.
// For the summed amplitudes of one frequency, (|omega|/2) h is the cycle-averaged heating [W m-3].
inline double c_volumetric_heating(const c_Tensor6& stress, const c_Tensor6& strain) noexcept
{
    double h = 0.0;
    for (size_t k = 0; k < 6; ++k)
    {
        const double term = stress.c[k].imag() * strain.c[k].real()
                          - stress.c[k].real() * strain.c[k].imag();
        h += (k < 3) ? term : 2.0 * term;
    }
    // The weighted sum is real. For the summed amplitudes of one frequency it is non-negative for dissipative
    // moduli, so there the magnitude equals the signed form below (Europa book Eq. 42).
    return std::abs(h);
}

// Signed form of c_volumetric_heating: sum_k w_k Im(sigma_k conj(eps_k)) with no abs(). The secular
// (cycle or orbit-averaged) heating calls it with the total complex-phasor stress and strain at one
// frequency (every wave at that |omega| summed first, since same-frequency cross terms survive the
// average); that frequency's volumetric heating is (|omega|/2) times this, and the frequencies add
// (distinct-frequency cross terms average to zero, so they are omitted).
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

// The 3 complex displacement components at a point: [0]=radial u_r, [1]=polar u_theta, [2]=azimuthal u_phi.
struct c_Vector3
{
    std::array<std::complex<double>, 3> c { };
};

// Tidal displacements at one point from the radial functions y1 (radial) and y3 (tangential) and the
// potential point (TB05 Eq. 9; SVC16):
//   u_r     = y1 * U
//   u_theta = y3 * dU/dtheta
//   u_phi   = y3 * dU/dphi / sin(theta)
// y1, y3 [s^2 m^-1] times U [m^2 s^-2] give metres. Templated like c_compute_strain_stress so the real
// (instantaneous) and complex-phasor potential points both work. At a pole (sin theta = 0) u_phi is 0 when
// dU/dphi is 0 (every m = 0 mode) and NaN otherwise.
template <typename PotentialPointT>
inline void c_compute_displacements(
        const std::complex<double>& y1,
        const std::complex<double>& y3,
        const PotentialPointT& P,
        double colatitude,
        c_Vector3& out) noexcept
{
    out.c[0] = y1 * P.U;
    out.c[1] = y3 * P.dU_dtheta;
    // The same pole test as c_colatitude_trig.
    const double sin_theta = std::sin(colatitude);
    if (std::abs(sin_theta) > TidalPyConstants::d_EPS)
    {
        out.c[2] = y3 * P.dU_dphi / sin_theta;
    }
    else if (std::abs(P.dU_dphi) == 0.0)
    {
        out.c[2] = std::complex<double>(0.0, 0.0);
    }
    else
    {
        const double nan_val = std::numeric_limits<double>::quiet_NaN();
        out.c[2] = std::complex<double>(nan_val, nan_val);
    }
}

}  // namespace tides
}  // namespace tidalpy
