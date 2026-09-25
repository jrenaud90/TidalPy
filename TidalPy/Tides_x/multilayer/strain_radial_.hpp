// strain_radial_.hpp - Radial coefficients of the 3D tidal strain tensor.
//
// The 3D tidal strain tensor factorizes, per mode, as
//     eps_c(r, theta, phi) = R_c(r) * A_c(theta, phi)
// where R_c(r) is a radial coefficient built from the radial-solver y-functions + the complex
// (viscoelastic) moduli, and A_c(theta, phi) is the tidal potential's angular factor (and its
// theta/phi derivatives). This header computes the radial coefficients R_c(r); the angular factors
// A_c live in the potentials submodule, and the kernel that combines them (and forms stress +
// heating) is built on top of both.
//
// References: The theta-phi and phi-phi forms use the Kervazo+2021 (A&A Appendix D) correction to
// the Tobie+2005 typo; that correction is purely angular, so it does not affect these radial coefficients.
#pragma once

#include <complex>
#include <limits>


namespace tidalpy {
namespace tides {

// The distinct radial multipliers the 6 strain components need (all complex), plus the complex
// moduli (carried through so the kernel can build the isotropic stress sigma = 2 mu eps + lame tr(eps)).
// Mapping to the strain components (U = tidal potential angular factor; subscripts = its derivatives):
//   eps_rr     = dy1_dr      * U
//   eps_thth   = y3_over_r   * d2U/dtheta2           + y1_over_r * U
//   eps_phph   = y1_over_r   * U                     + y3_over_r * [ sin^-2(theta) d2U/dphi2 + cot(theta) dU/dtheta ]
//   eps_rtheta = y4_over_2mu * dU/dtheta
//   eps_rphi   = y4_over_2mu * sin^-1(theta) dU/dphi
//   eps_thphi  = y3_over_2r  * [ 2 ( d2U/dtheta_dphi - cot(theta) dU/dphi ) sin^-1(theta) ]
struct c_StrainRadialCoeffs
{
    c_StrainRadialCoeffs() {}
    c_StrainRadialCoeffs(
        std::complex<double> dy1_dr_,
        std::complex<double> y1_over_r_,
        std::complex<double> y3_over_r_,
        std::complex<double> y4_over_2mu_,
        std::complex<double> y3_over_2r_,
        std::complex<double> shear_,
        std::complex<double> lame_,
        bool valid_):
            dy1_dr(dy1_dr_),
            y1_over_r(y1_over_r_),
            y3_over_r(y3_over_r_),
            y4_over_2mu(y4_over_2mu_),
            y3_over_2r(y3_over_2r_),
            shear(shear_),
            lame(lame_),
            valid(valid_) {}

    std::complex<double> dy1_dr      {0.0, 0.0};   // eps_rr
    std::complex<double> y1_over_r   {0.0, 0.0};   // eps_thth, eps_phph (the U term)
    std::complex<double> y3_over_r   {0.0, 0.0};   // eps_thth (d2U/dtheta2), eps_phph (bracket)
    std::complex<double> y4_over_2mu {0.0, 0.0};   // eps_rtheta, eps_rphi
    std::complex<double> y3_over_2r  {0.0, 0.0};   // eps_thphi
    std::complex<double> shear       {0.0, 0.0};   // complex mu  (for the stress constitutive law)
    std::complex<double> lame        {0.0, 0.0};   // complex lambda = kappa - 2/3 mu
    bool valid {false};                            // false -> liquid / degenerate (no shear kernel here)
};

// Compute the radial strain coefficients at one radius for one (l, omega) radial solution.
//
// Inputs:
//   y1, y2, y3, y4 : the collapsed radial-solver y-functions at this radius (complex; from the dense
//                    calling system, evaluated at the mode's frequency). y2 is only needed for the
//                    solid-compressible dy1/dr.
//   shear, bulk    : complex viscoelastic moduli at this radius and the mode's frequency.
//   radius         : radius [m] (must be > 0).
//   degree_l       : harmonic degree l (the radial problem depends on l, not on m).
//   layer_is_solid, layer_is_incompressible : layer assumptions for the layer containing this radius.
inline c_StrainRadialCoeffs c_compute_strain_radial_coeffs(
        const std::complex<double>& y1,
        const std::complex<double>& y2,
        const std::complex<double>& y3,
        const std::complex<double>& y4,
        const std::complex<double>& shear,
        const std::complex<double>& bulk,
        double radius,
        double degree_l,
        bool layer_is_solid,
        bool layer_is_incompressible) noexcept
{
    c_StrainRadialCoeffs out;
    const std::complex<double> c_zero(0.0, 0.0);

    out.shear = shear;
    out.lame  = bulk - (2.0 / 3.0) * shear;

    // Solid-only shear kernel: a liquid (mu = 0) or the singular center (r = 0) has no well-defined
    // shear strain here. Mark invalid and NaN-fill the multipliers.
    if (!layer_is_solid || shear == c_zero || !(radius > 0.0))
    {
        const double nan_val = std::numeric_limits<double>::quiet_NaN();
        const std::complex<double> c_nan(nan_val, nan_val);
        out.dy1_dr = out.y1_over_r = out.y3_over_r = out.y4_over_2mu = out.y3_over_2r = c_nan;
        out.valid = false;
        return out;
    }

    const double r_inv = 1.0 / radius;
    const double llp1  = degree_l * (degree_l + 1.0);
    const std::complex<double> y1_y3_term = 2.0 * y1 - llp1 * y3;   // (2 y1 - l(l+1) y3), as in the ODEs

    // dy1/dr - layer-type specific (mirrors RadialSolver_x/derivatives/odes_.hpp exactly).
    if (layer_is_incompressible)
    {
        // Solid incompressible (static & dynamic): dy1/dr = -(2 y1 - l(l+1) y3)/r.
        out.dy1_dr = -y1_y3_term * r_inv;
    }
    else
    {
        // Solid compressible (static & dynamic): dy1/dr = (1/(lame+2mu)) [ y2 - (lame/r)(2 y1 - l(l+1) y3) ].
        const std::complex<double> lame_2mu = out.lame + 2.0 * shear;
        out.dy1_dr = (1.0 / lame_2mu) * (y2 - out.lame * r_inv * y1_y3_term);
    }

    out.y1_over_r   = y1 * r_inv;
    out.y3_over_r   = y3 * r_inv;
    out.y4_over_2mu = y4 / (2.0 * shear);
    out.y3_over_2r  = y3 * (0.5 * r_inv);
    out.valid = true;
    return out;
}

}  // namespace tides
}  // namespace tidalpy

