// common_.hpp - Spherical Bessel z function and Takeuchi phi/psi functions
//
// References
// ----------
// TS72: Takeuchi & Saito (1972)
// KMN15: Kamata, Matsuyama, & Nimmo (2015)
#pragma once

#include <cmath>
#include <complex>
#include "xsf/bessel.h"
#include "xsf/sph_bessel.h"

// Calculates the z function from the spherical Bessel functions.
//
// References
// ----------
// TS72 Eqs. 96, 97; KMN15 Eq. B14
//
// Parameters
// ----------
// x_squared : complex
//     Expression passed to the Bessel function.
// degree_l : int
//     Tidal harmonic order.
//
// Returns
// -------
// complex : z
inline std::complex<double> c_z_calc(
        const std::complex<double>& x_squared,
        const int degree_l) noexcept
{
    // The recursion of TS72 Eq. 97, z_l = x^2 / (2l + 3 - z_{l+1}), is the spherical Bessel recurrence for this same
    // ratio and agrees with it to about 1e-13 (checked for l = 2..10, real and complex x^2 up to 30; Issue #42).
    // The direct ratio is used because it needs no truncation depth.

    // Taylor expansion works well when x_squared is small and is faster.
    const double l_dbl = static_cast<double>(degree_l);
    const double l2 = l_dbl * 2.0;

    std::complex<double> z;

    if (std::abs(x_squared) > 0.1)
    {
        const std::complex<double> x = std::sqrt(x_squared);
        z = x * xsf::sph_bessel_j(degree_l + 1, x) / xsf::sph_bessel_j(degree_l, x);
    } else
    {
        // Use Taylor series; JPR derived this on 2024-02-05
        const double l2_3   = l2 + 3.0;
        const double l2_3sq = l2_3 * l2_3;
        const double l2_3cb = l2_3sq * l2_3;
        const double l2_5   = l2 + 5.0;
        const double l2_5sq = l2_5 * l2_5;
        const double l2_7   = l2 + 7.0;
        const double l2_9   = l2 + 9.0;
        const double l2_11  = l2 + 11.0;
        // Powers of u = x^2: the series runs u, u^2, u^3, u^4, u^5.
        const std::complex<double> u_2 = x_squared * x_squared;
        const std::complex<double> u_3 = u_2 * x_squared;
        const std::complex<double> u_4 = u_2 * u_2;
        const std::complex<double> u_5 = u_4 * x_squared;

        z = (x_squared                   / l2_3) +
            (u_2                         / (l2_3sq * l2_5)) +
            (u_3 * 2.0                   / (l2_3cb * l2_5 * l2_7)) +
            (u_4 * (27.0 + 10.0 * l_dbl) / (l2_3cb * l2_3 * l2_5sq * l2_7 * l2_9)) +
            (u_5 * (90.0 + 28.0 * l_dbl) / (l2_3cb * l2_3sq * l2_5sq * l2_7 * l2_9 * l2_11));
    }

    return z;
}


// Calculate phi, phi_{l+1}, and psi functions used to find initial conditions for shooting method.
//
// TS72 Eq. 103 defines phi_l(x) = (2l + 1)!! j_l(x) / x^l and psi_l(x) = 2 (2l + 3) (1 - phi_l(x)) / x^2, with
// z2 = x^2. psi's 1 - phi_l cancels as z2 -> 0, so small |z2| uses a series (Issue #41). The series is exact to
// rounding below |z2| ~ 0.3, but its truncation error grows quickly beyond (1e-11 at 1, 2e-5 at 10, useless past
// ~30), and |z2| reaches those values at high degree, at short periods, and in dynamic liquids at tidal periods.
// Above |z2| = 0.5 the definitions are used instead, which hold about 1e-13 from there up.
//
// References
// ----------
// TS72 Eq. 103
//
// Parameters
// ----------
// z2 : complex
//     z^2 argument
// degree_l : int
//     Tidal harmonic order
// phi_ptr : complex*, output
// phi_lplus1_ptr : complex*, output
// psi_ptr : complex*, output
inline void c_takeuchi_phi_psi(
        const std::complex<double>& z2,
        const int degree_l,
        std::complex<double>* phi_ptr,
        std::complex<double>* phi_lplus1_ptr,
        std::complex<double>* psi_ptr) noexcept
{
    if (std::abs(z2) > 0.5)
    {
        const std::complex<double> x = std::sqrt(z2);
        // (2l + 1)!! / x^l, built one factor at a time so neither part overflows at high degree.
        std::complex<double> dfact_over_xl(1.0, 0.0);
        for (int k = 1; k <= degree_l; ++k)
        {
            dfact_over_xl *= (2.0 * static_cast<double>(k) + 1.0) / x;
        }
        const double dlp3 = 2.0 * static_cast<double>(degree_l) + 3.0;
        phi_ptr[0]        = dfact_over_xl * xsf::sph_bessel_j(degree_l, x);
        phi_lplus1_ptr[0] = dfact_over_xl * (dlp3 / x) * xsf::sph_bessel_j(degree_l + 1, x);
        psi_ptr[0]        = 2.0 * dlp3 * (1.0 - phi_ptr[0]) / z2;
        return;
    }

    const std::complex<double> z4   = z2 * z2;
    const std::complex<double> z6   = z4 * z2;
    const std::complex<double> z8   = z6 * z2;
    const std::complex<double> z10  = z8 * z2;
    const double l    = static_cast<double>(degree_l);
    const double l2   = 2.0 * l;
    const double l_3  = 3.0  + l2;
    const double l_5  = 5.0  + l2;
    const double l_7  = 7.0  + l2;
    const double l_9  = 9.0  + l2;
    const double l_11 = 11.0 + l2;
    const double l_13 = 13.0 + l2;

    // Series expansion was done in `Takeuchi Starting Conditions.nb` on 2024-01-31 by JPR.
    phi_ptr[0] = (
         1.0 +
        -z2  / (2.0    * l_3) +
         z4  / (8.0    * l_3 * l_5) +
        -z6  / (48.0   * l_3 * l_5 * l_7) +
         z8  / (384.0  * l_3 * l_5 * l_7 * l_9) +
        -z10 / (3840.0 * l_3 * l_5 * l_7 * l_9 * l_11)
    );

    // Note that the l+1 function skips a denominator term compared to the function above.
    phi_lplus1_ptr[0] = (
         1.0 +
        -z2  / (2.0    * l_5) +
         z4  / (8.0    * l_5 * l_7) +
        -z6  / (48.0   * l_5 * l_7 * l_9) +
         z8  / (384.0  * l_5 * l_7 * l_9 * l_11) +
        -z10 / (3840.0 * l_5 * l_7 * l_9 * l_11 * l_13)
    );

    psi_ptr[0] = (
         1.0 +
        -z2  / (4.0     * l_5) +
        // TS72 prints this term with 12 in place of 24; the definition of psi above gives 24.
         z4  / (24.0    * l_5 * l_7) +
        -z6  / (192.0   * l_5 * l_7 * l_9) +
         z8  / (1920.0  * l_5 * l_7 * l_9 * l_11) +
        -z10 / (23040.0 * l_5 * l_7 * l_9 * l_11 * l_13)
     );
}
