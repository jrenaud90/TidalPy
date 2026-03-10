// common_.hpp - Spherical Bessel z function and Takeuchi phi/psi functions
// Ported from TidalPy/RadialSolver/starting/common.pyx
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

// Calculates the z function using spherical Bessel function, see Eq. B14 of KMN15.
//
// References
// ----------
// TS72 Eqs. 96, 97
// KMN15 Eq. B14
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
// z : complex
//     Result
inline std::complex<double> c_z_calc(
        const std::complex<double>& x_squared,
        const int degree_l) noexcept
{
    // QUESTION: (Issue #42) The recursion formula shown in TS72 and KMN15 (See TS72 Eq. 97) does not reproduce
    // the full version of this function (TS72 Eq. 96). Off by ~30% at l=2, improves to ~10% at l=10.
    // We do not use the recursive formula here.

    // Taylor expansion works well when x_squared is small and is faster.
    const double l_dbl = static_cast<double>(degree_l);
    const double l2 = l_dbl * 2.0;

    std::complex<double> z;

    if (std::abs(x_squared) > 0.1)
    {
        // Use full spherical Bessel function
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
        const std::complex<double> x_sq2 = x_squared * x_squared;
        const std::complex<double> x_sq4 = x_sq2 * x_sq2;
        const std::complex<double> x_sq6 = x_sq4 * x_sq2;
        const std::complex<double> x_sq8 = x_sq6 * x_sq2;

        z = (x_squared                         / l2_3) +
            (x_sq2                             / (l2_3sq * l2_5)) +
            (x_sq4 * 2.0                       / (l2_3cb * l2_5 * l2_7)) +
            (x_sq6 * (27.0 + 10.0 * l_dbl)     / (l2_3cb * l2_3 * l2_5sq * l2_7 * l2_9)) +
            (x_sq8 * (90.0 + 28.0 * l_dbl)     / (l2_3cb * l2_3sq * l2_5sq * l2_7 * l2_9 * l2_11));
    }

    return z;
}


// Calculate phi, phi_{l+1}, and psi functions used to find initial conditions for shooting method.
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
    // Floating point errors prevent us from using the exact definition of these functions (See Issue #41)
    // Instead we use the limiting version of the functions.
    // However, we leave the full definition at the bottom of this function for reference.

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
        // NOTE: Error? TS72 quotes the next term as z4 / (12. * (5. + 2. * l) * (7 + 2. * l)); factor of 2 off in denom.
         z4  / (24.0    * l_5 * l_7) +
        -z6  / (192.0   * l_5 * l_7 * l_9) +
         z8  / (1920.0  * l_5 * l_7 * l_9 * l_11) +
        -z10 / (23040.0 * l_5 * l_7 * l_9 * l_11 * l_13)
     );

    /* Full version
    # cdef char lp1 = degree_l + 1
    # cdef double l_dbl_factorial   = cf_double_factorial(2 * degree_l + 1)
    # cdef double lp1_dbl_factorial = cf_double_factorial(2 * lp1 + 1)
    
    # cdef double complex z    = cf_csqrt(z2)
    # cdef double complex zl   = cf_cipow(z, degree_l)
    # cdef double complex zlp1 = cf_cipow(z, lp1)

    # phi_ptr[0]        = l_dbl_factorial * (spherical_jn(degree_l, z) / zl)
    # phi_lplus1_ptr[0] = lp1_dbl_factorial * (spherical_jn(lp1, z) / zlp1)
    # psi_ptr[0]        = (2. * (2. * degree_l + 3.) / (z * z)) * (1. - phi_ptr[0])

    # DEBUG: Uncomment the below and comment the above to use the method that TidalPy v0.4.0 used for this function.
    # Left here for comparison/debug purpouses.
    # cdef double degree_lp1 = degree_l + 1.
    # cdef double complex z4 = z2 * z2

    # phi_ptr[0] = 1. - \
    #       z2 / (2. * (2. * degree_l + 3.)) + \
    #       z4 / (8. * (2. * degree_l + 3.) * (2. * degree_l + 5.))
    # phi_lplus1_ptr[0] = 1. - \
    #              z2 / (2. * (2. * degree_lp1 + 3.)) + \
    #              z4 / (8. * (2. * degree_lp1 + 3.) * (2. * degree_lp1 + 5.))
    # psi_ptr[0] = 1. - \
    #       z2 / (4. * (2. * degree_l + 5.)) + \
    #       z4 / (12. * (2. * degree_l + 5.) * (2. * degree_l + 7.))
    */
}
