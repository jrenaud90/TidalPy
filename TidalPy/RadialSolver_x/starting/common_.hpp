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

// Compute spherical Bessel function of the first kind j_n(z) for complex argument
// using forward recurrence:
//   j_0(z) = sin(z)/z
//   j_1(z) = sin(z)/z^2 - cos(z)/z
//   j_{n+1}(z) = (2n+1)/z * j_n(z) - j_{n-1}(z)
inline std::complex<double> c_spherical_jn(
        int n,
        std::complex<double> z) noexcept
{
    // Handle z ~ 0 edge case
    if (std::abs(z) < 1.0e-30) {
        if (n == 0) return std::complex<double>(1.0, 0.0);
        return std::complex<double>(0.0, 0.0);
    }

    std::complex<double> j0 = std::sin(z) / z;
    if (n == 0) return j0;

    std::complex<double> j1 = std::sin(z) / (z * z) - std::cos(z) / z;
    if (n == 1) return j1;

    // Forward recurrence
    std::complex<double> j_prev = j0;
    std::complex<double> j_curr = j1;
    for (int k = 1; k < n; ++k) {
        std::complex<double> j_next = (static_cast<double>(2 * k + 1) / z) * j_curr - j_prev;
        j_prev = j_curr;
        j_curr = j_next;
    }
    return j_curr;
}


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
        std::complex<double> x_squared,
        int degree_l) noexcept
{
    // QUESTION: (Issue #42) The recursion formula shown in TS72 and KMN15 (See TS72 Eq. 97) does not reproduce
    // the full version of this function (TS72 Eq. 96). Off by ~30% at l=2, improves to ~10% at l=10.
    // We do not use the recursive formula here.

    // Taylor expansion works well when x_squared is small and is faster.
    double l2 = static_cast<double>(degree_l) * 2.0;

    std::complex<double> z;

    if (std::abs(x_squared) > 0.1) {
        // Use full spherical Bessel function
        std::complex<double> x = std::sqrt(x_squared);
        z = x * c_spherical_jn(degree_l + 1, x) / c_spherical_jn(degree_l, x);
    } else {
        // Use Taylor series; JPR derived this on 2024-02-05
        double l2_3  = l2 + 3.0;
        double l2_5  = l2 + 5.0;
        double l2_7  = l2 + 7.0;
        double l2_9  = l2 + 9.0;
        double l2_11 = l2 + 11.0;

        std::complex<double> x2 = x_squared;
        std::complex<double> x4 = x2 * x2;
        std::complex<double> x6 = x4 * x2;  // x_squared**4 in original = x^8? No, look carefully at original:
        // Original uses x_squared**2 = (x^2)^2 = x^4, x_squared**4 = (x^2)^4 = x^8, etc.
        // But that doesn't match. Let me re-read...
        // Line 63: x_squared**2 / (l2_3^2 * l2_5) -- this is x_squared^2
        // Line 64: x_squared**4 * 2. / (...) -- but this looks wrong, should be x_squared^3
        // Actually looking at the pattern more carefully, the original code has a specific Taylor expansion
        // Let me just port it exactly as written:
        //   x_squared / l2_3
        //   + x_squared^2 / (l2_3^2 * l2_5)
        //   + x_squared^4 * 2 / (l2_3^3 * l2_5 * l2_7)       -- NOTE: x_squared**4 in original Cython = (x_squared)^4
        //   + x_squared^6 * (27 + 10*l) / (l2_3^4 * l2_5^2 * l2_7 * l2_9)
        //   + x_squared^8 * (90 + 28*l) / (l2_3^5 * l2_5^2 * l2_7 * l2_9 * l2_11)
        // OPT: The exponents in the original code seem unusual (jumping from ^2 to ^4 to ^6 to ^8).
        // Porting exactly as written.
        double l_dbl = static_cast<double>(degree_l);
        std::complex<double> x_sq2 = x_squared * x_squared;
        std::complex<double> x_sq4 = x_sq2 * x_sq2;
        std::complex<double> x_sq6 = x_sq4 * x_sq2;
        std::complex<double> x_sq8 = x_sq6 * x_sq2;

        z = (x_squared                         / l2_3) +
            (x_sq2                              / (l2_3 * l2_3 * l2_5)) +
            (x_sq4 * 2.0                        / (l2_3 * l2_3 * l2_3 * l2_5 * l2_7)) +
            (x_sq6 * (27.0 + 10.0 * l_dbl)     / (l2_3 * l2_3 * l2_3 * l2_3 * l2_5 * l2_5 * l2_7 * l2_9)) +
            (x_sq8 * (90.0 + 28.0 * l_dbl)     / (l2_3 * l2_3 * l2_3 * l2_3 * l2_3 * l2_5 * l2_5 * l2_7 * l2_9 * l2_11));
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
        std::complex<double> z2,
        int degree_l,
        std::complex<double>* phi_ptr,
        std::complex<double>* phi_lplus1_ptr,
        std::complex<double>* psi_ptr) noexcept
{
    // Floating point errors prevent us from using the exact definition of these functions (See Issue #41)
    // Instead we use the limiting version of the functions.

    std::complex<double> z4   = z2 * z2;
    std::complex<double> z6   = z4 * z2;
    std::complex<double> z8   = z6 * z2;
    std::complex<double> z10  = z8 * z2;
    double l    = static_cast<double>(degree_l);
    double l2   = 2.0 * l;
    double l_3  = 3.0  + l2;
    double l_5  = 5.0  + l2;
    double l_7  = 7.0  + l2;
    double l_9  = 9.0  + l2;
    double l_11 = 11.0 + l2;
    double l_13 = 13.0 + l2;

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
        // RECORD: Error? TS72 quotes the next term as z4 / (12. * (5. + 2. * l) * (7 + 2. * l)); factor of 2 off in denom.
         z4  / (24.0    * l_5 * l_7) +
        -z6  / (192.0   * l_5 * l_7 * l_9) +
         z8  / (1920.0  * l_5 * l_7 * l_9 * l_11) +
        -z10 / (23040.0 * l_5 * l_7 * l_9 * l_11 * l_13)
     );
}
