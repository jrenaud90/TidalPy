// solid_matrix_.hpp - Fundamental matrix and its inverse
// Ported from TidalPy/RadialSolver/matrix_types/solid_matrix.pyx
//
// Fundamental Matrix and its inverse as defined in Sabadini, Vermeerson, & Cambiotti (2016)
//
// Assumptions
// -----------
// These matrices assume an incompressible body.
//
// References
// ----------
// SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
// HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
// ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
#pragma once

#include <complex>
#include <cmath>

#include "../../constants_.hpp"


/// Construct the fundamental matrix and its inverse for a generic order-l.
///
/// See Eq. 2.42 of SVC16
///
/// Assumptions:
///   These matrices assume an incompressible body.
///
/// Parameters
/// ----------
/// first_slice_index : size_t
///     Initial radial index to start populating matrices at.
/// num_radial_slices : size_t
///     Number of radial slices.
/// radius_array_ptr : double*
///     Pointer to array of Radius values [m].
/// density_array_ptr : double*
///     Pointer to array of Density at each radius [kg m-3].
/// gravity_array_ptr : double*
///     Pointer to array of acceleration due to gravity at each radius [m s-2].
/// complex_shear_array_ptr : std::complex<double>*
///     Pointer to array of Complex shear modulus at each radius [Pa].
/// fundamental_mtx_ptr : std::complex<double>*
///     (Output) 6x6 fundamental matrix per slice.
/// inverse_fundamental_mtx_ptr : std::complex<double>*
///     (Output) Inverse of the fundamental matrix per slice.
/// derivative_mtx_ptr : std::complex<double>*
///     (Output) The matrix A such that dy/dr = A * y.
/// degree_l : int
///     Harmonic degree.
/// G_to_use : double
///     Gravitational constant.
inline void c_fundamental_matrix(
    size_t first_slice_index,
    size_t num_radial_slices,
    double* radius_array_ptr,
    double* density_array_ptr,
    double* gravity_array_ptr,
    std::complex<double>* complex_shear_array_ptr,
    std::complex<double>* fundamental_mtx_ptr,
    std::complex<double>* inverse_fundamental_mtx_ptr,
    std::complex<double>* derivative_mtx_ptr,
    int degree_l,
    double G_to_use) noexcept
{
    // Degree-l optimizations
    const double degree_l_dbl = static_cast<double>(degree_l);
    const double dlm1         = 2.0 * degree_l_dbl - 1.0;
    const double l2p3lm1      = degree_l_dbl * degree_l_dbl + 3.0 * degree_l_dbl - 1.0;
    const double l2mlm3       = degree_l_dbl * degree_l_dbl - degree_l_dbl - 3.0;
    const double lp1          = degree_l_dbl + 1.0;
    const double lp2          = degree_l_dbl + 2.0;
    const double lp3          = degree_l_dbl + 3.0;
    const double l2m1         = degree_l_dbl * degree_l_dbl - 1.0;
    const double dlp1         = 2.0 * degree_l_dbl + 1.0;
    const double dlp3         = 2.0 * degree_l_dbl + 3.0;
    const double dlp1_inverse = 1.0 / dlp1;

    for (size_t slice_i = first_slice_index; slice_i < num_radial_slices; ++slice_i)
    {
        // Shift index by 36 (for the inner 6x6 matrix)
        const size_t index_shift = slice_i * 36;

        // Unpack radially dependent variables
        const double radius  = radius_array_ptr[slice_i];
        const double gravity = gravity_array_ptr[slice_i];
        const double density = density_array_ptr[slice_i];
        const std::complex<double>complex_shear = complex_shear_array_ptr[slice_i];

        // Radius-based optimizations
        const double r_inv  = 1.0 / radius;
        const double rl     = std::pow(radius, degree_l_dbl);
        const double rlp1   = std::pow(radius, degree_l_dbl + 1.0);
        const double rlp2   = std::pow(radius, degree_l_dbl + 2.0);
        const double rlp3   = std::pow(radius, degree_l_dbl + 3.0);
        const double rnl    = std::pow(radius, -degree_l_dbl);
        const double rnlm2  = std::pow(radius, -degree_l_dbl - 2.0);
        const double rlm1   = std::pow(radius, degree_l_dbl - 1.0);
        const double rgp    = radius * gravity * density;
        const double piGp   = TidalPyConstants::d_PI * G_to_use * density;
        const std::complex<double> mu_inv = 1.0 / complex_shear;
        const std::complex<double> rgp_s  = rgp * mu_inv;
        const std::complex<double> r_s    = radius * mu_inv;
        const std::complex<double> pr_s   = density * r_s;

        // Build Fundamental Matrix (Eq. 2.42 in SVC16)
        // Row 1
        fundamental_mtx_ptr[index_shift + 0]  = degree_l_dbl * rlp1 / (2.0 * dlp3);
        fundamental_mtx_ptr[index_shift + 1]  = rlm1;
        fundamental_mtx_ptr[index_shift + 2]  = 0.0;
        fundamental_mtx_ptr[index_shift + 3]  = lp1 * rnl / (2.0 * dlm1);
        fundamental_mtx_ptr[index_shift + 4]  = rnlm2;
        fundamental_mtx_ptr[index_shift + 5]  = 0.0;

        // Row 2
        fundamental_mtx_ptr[index_shift + 6]  = lp3 * rlp1 / (2.0 * dlp3 * lp1);
        fundamental_mtx_ptr[index_shift + 7]  = rlm1 / degree_l_dbl;
        fundamental_mtx_ptr[index_shift + 8]  = 0.0;
        fundamental_mtx_ptr[index_shift + 9]  = (2.0 - degree_l_dbl) * rnl / (2.0 * degree_l_dbl * dlm1);
        fundamental_mtx_ptr[index_shift + 10] = -rnlm2 / lp1;
        fundamental_mtx_ptr[index_shift + 11] = 0.0;

        // Row 3
        // RECORD: Believe there is a typo in HH14, they have the radius^l only on one term instead of both.
        fundamental_mtx_ptr[index_shift + 12] = (degree_l_dbl * rgp + 2.0 * l2mlm3 * complex_shear) * rl / (2.0 * dlp3);
        fundamental_mtx_ptr[index_shift + 13] = (rgp + 2.0 * (degree_l_dbl - 1.0) * complex_shear) * std::pow(radius, degree_l_dbl - 2.0);
        fundamental_mtx_ptr[index_shift + 14] = -density * rl;
        fundamental_mtx_ptr[index_shift + 15] = (lp1 * rgp - 2.0 * l2p3lm1 * complex_shear) / (2.0 * dlm1 * rlp1);
        fundamental_mtx_ptr[index_shift + 16] = (rgp - 2.0 * lp2 * complex_shear) / rlp3;
        fundamental_mtx_ptr[index_shift + 17] = -density / rlp1;

        // Row 4
        fundamental_mtx_ptr[index_shift + 18] = degree_l_dbl * lp2 * complex_shear * rl / (dlp3 * lp1);
        fundamental_mtx_ptr[index_shift + 19] = 2.0 * (degree_l_dbl - 1.0) * complex_shear * std::pow(radius, degree_l_dbl - 2.0) / degree_l_dbl;
        fundamental_mtx_ptr[index_shift + 20] = 0.0;
        fundamental_mtx_ptr[index_shift + 21] = l2m1 * complex_shear / (degree_l_dbl * dlm1 * rlp1);
        fundamental_mtx_ptr[index_shift + 22] = 2.0 * lp2 * complex_shear / (lp1 * rlp3);
        fundamental_mtx_ptr[index_shift + 23] = 0.0;

        // Row 5
        fundamental_mtx_ptr[index_shift + 24] = 0.0;
        fundamental_mtx_ptr[index_shift + 25] = 0.0;
        fundamental_mtx_ptr[index_shift + 26] = -rl;
        fundamental_mtx_ptr[index_shift + 27] = 0.0;
        fundamental_mtx_ptr[index_shift + 28] = 0.0;
        fundamental_mtx_ptr[index_shift + 29] = -1.0 / rlp1;

        // Row 6
        fundamental_mtx_ptr[index_shift + 30] = 2.0 * piGp * degree_l_dbl * rlp1 / dlp3;
        fundamental_mtx_ptr[index_shift + 31] = 4.0 * piGp * rlm1;
        fundamental_mtx_ptr[index_shift + 32] = -dlp1 * rlm1;
        fundamental_mtx_ptr[index_shift + 33] = 2.0 * piGp * lp1 / (dlm1 * rl);
        fundamental_mtx_ptr[index_shift + 34] = 4.0 * piGp / rlp2;
        fundamental_mtx_ptr[index_shift + 35] = 0.0;

        // Inverse of the Fundamental Matrix
        // From SVC16 Eq. 2.45: Fundamental Inverse = D_Mtx * Y^Bar_Mtx
        // D_Mtx is a diagonal matrix. We multiply first and write down the product directly.

        // D Coefficients
        const double d_coeff_1 = dlp1_inverse * lp1 / rlp1;
        const double d_coeff_2 = dlp1_inverse * degree_l_dbl * lp1 / (2.0 * dlm1 * rlm1);
        const double d_coeff_3 = dlp1_inverse * 1.0 / rlm1;
        const double d_coeff_4 = dlp1_inverse * degree_l_dbl * rl;
        const double d_coeff_5 = dlp1_inverse * rlp2 * degree_l_dbl * lp1 / (2.0 * dlp3);
        const double d_coeff_6 = dlp1_inverse * -rlp1;

        // Row 1
        inverse_fundamental_mtx_ptr[index_shift + 0]  = d_coeff_1 * (rgp_s - 2.0 * lp2);
        inverse_fundamental_mtx_ptr[index_shift + 1]  = d_coeff_1 * (2.0 * degree_l_dbl * lp2);
        inverse_fundamental_mtx_ptr[index_shift + 2]  = d_coeff_1 * (-r_s);
        inverse_fundamental_mtx_ptr[index_shift + 3]  = d_coeff_1 * (degree_l_dbl * r_s);
        inverse_fundamental_mtx_ptr[index_shift + 4]  = d_coeff_1 * (pr_s);
        inverse_fundamental_mtx_ptr[index_shift + 5]  = 0.0;

        // Row 2
        inverse_fundamental_mtx_ptr[index_shift + 6]  = d_coeff_2 * (-rgp_s + 2.0 * l2p3lm1 / lp1);
        inverse_fundamental_mtx_ptr[index_shift + 7]  = d_coeff_2 * (-2.0 * l2m1);
        inverse_fundamental_mtx_ptr[index_shift + 8]  = d_coeff_2 * (r_s);
        inverse_fundamental_mtx_ptr[index_shift + 9]  = d_coeff_2 * ((2.0 - degree_l_dbl) * r_s);
        inverse_fundamental_mtx_ptr[index_shift + 10] = d_coeff_2 * (-pr_s);
        inverse_fundamental_mtx_ptr[index_shift + 11] = 0.0;

        // Row 3
        inverse_fundamental_mtx_ptr[index_shift + 12] = d_coeff_3 * (4.0 * piGp);
        inverse_fundamental_mtx_ptr[index_shift + 13] = 0.0;
        inverse_fundamental_mtx_ptr[index_shift + 14] = 0.0;
        inverse_fundamental_mtx_ptr[index_shift + 15] = 0.0;
        inverse_fundamental_mtx_ptr[index_shift + 16] = 0.0;
        inverse_fundamental_mtx_ptr[index_shift + 17] = -d_coeff_3;

        // Row 4
        inverse_fundamental_mtx_ptr[index_shift + 18] = d_coeff_4 * (rgp_s + 2.0 * (degree_l_dbl - 1.0));
        inverse_fundamental_mtx_ptr[index_shift + 19] = d_coeff_4 * (2.0 * l2m1);
        inverse_fundamental_mtx_ptr[index_shift + 20] = d_coeff_4 * (-r_s);
        inverse_fundamental_mtx_ptr[index_shift + 21] = d_coeff_4 * (-lp1 * r_s);
        inverse_fundamental_mtx_ptr[index_shift + 22] = d_coeff_4 * (pr_s);
        inverse_fundamental_mtx_ptr[index_shift + 23] = 0.0;

        // Row 5
        inverse_fundamental_mtx_ptr[index_shift + 24] = d_coeff_5 * (-rgp_s - 2.0 * l2mlm3 / degree_l_dbl);
        inverse_fundamental_mtx_ptr[index_shift + 25] = d_coeff_5 * (-2.0 * degree_l_dbl * lp2);
        inverse_fundamental_mtx_ptr[index_shift + 26] = d_coeff_5 * (r_s);
        inverse_fundamental_mtx_ptr[index_shift + 27] = d_coeff_5 * (lp3 * r_s);
        inverse_fundamental_mtx_ptr[index_shift + 28] = d_coeff_5 * (-pr_s);
        inverse_fundamental_mtx_ptr[index_shift + 29] = 0.0;

        // Row 6
        inverse_fundamental_mtx_ptr[index_shift + 30] = d_coeff_6 * (4.0 * piGp * radius);
        inverse_fundamental_mtx_ptr[index_shift + 31] = 0.0;
        inverse_fundamental_mtx_ptr[index_shift + 32] = 0.0;
        inverse_fundamental_mtx_ptr[index_shift + 33] = 0.0;
        inverse_fundamental_mtx_ptr[index_shift + 34] = d_coeff_6 * dlp1;
        inverse_fundamental_mtx_ptr[index_shift + 35] = d_coeff_6 * (-radius);

        // Build derivative matrix
        // Defined in SV04 -- Only valid for the incompressible case.
        // See SVC16 Eq. 1.95
        //    Note: the lambda in SVC16 is defined as bulk_mod - (2/3)*shear (Eq. 1.77; 2nd Lame parameter),
        //    for the incompressible assumption we assume the ratio (lambda/beta) -> 1 as K -> inf
        //    See SVC16 Eq. 1.95 for a compressible version. Take limit as K->inf to find below.

        // Row 1
        derivative_mtx_ptr[index_shift + 0]  = -2.0 * r_inv;
        derivative_mtx_ptr[index_shift + 1]  = degree_l_dbl * lp1 * r_inv;
        derivative_mtx_ptr[index_shift + 2]  = 0.0;
        derivative_mtx_ptr[index_shift + 3]  = 0.0;
        derivative_mtx_ptr[index_shift + 4]  = 0.0;
        derivative_mtx_ptr[index_shift + 5]  = 0.0;

        // Row 2
        derivative_mtx_ptr[index_shift + 6]  = -1.0 * r_inv;
        derivative_mtx_ptr[index_shift + 7]  = r_inv;
        derivative_mtx_ptr[index_shift + 8]  = 0.0;
        derivative_mtx_ptr[index_shift + 9]  = mu_inv;
        derivative_mtx_ptr[index_shift + 10] = 0.0;
        derivative_mtx_ptr[index_shift + 11] = 0.0;

        // Row 3
        derivative_mtx_ptr[index_shift + 12] = (4.0 * r_inv) * (3.0 * complex_shear * r_inv - density * gravity);
        derivative_mtx_ptr[index_shift + 13] = -degree_l_dbl * lp1 * r_inv * (6.0 * complex_shear * r_inv - density * gravity);
        derivative_mtx_ptr[index_shift + 14] = 0.0;
        derivative_mtx_ptr[index_shift + 15] = degree_l_dbl * lp1 * r_inv;
        derivative_mtx_ptr[index_shift + 16] = -density * lp1 * r_inv;
        derivative_mtx_ptr[index_shift + 17] = density;

        // Row 4
        derivative_mtx_ptr[index_shift + 18] = (-1.0 * r_inv) * (6.0 * complex_shear * r_inv - density * gravity);
        derivative_mtx_ptr[index_shift + 19] = 2.0 * (2.0 * degree_l_dbl * degree_l_dbl + 2.0 * degree_l_dbl - 1.0) * complex_shear * (r_inv * r_inv);
        derivative_mtx_ptr[index_shift + 20] = -r_inv;
        derivative_mtx_ptr[index_shift + 21] = -3.0 * r_inv;
        derivative_mtx_ptr[index_shift + 22] = density * r_inv;
        derivative_mtx_ptr[index_shift + 23] = 0.0;

        // Row 5
        derivative_mtx_ptr[index_shift + 24] = -4.0 * piGp;
        derivative_mtx_ptr[index_shift + 25] = 0.0;
        derivative_mtx_ptr[index_shift + 26] = 0.0;
        derivative_mtx_ptr[index_shift + 27] = 0.0;
        derivative_mtx_ptr[index_shift + 28] = -lp1 * r_inv;
        derivative_mtx_ptr[index_shift + 29] = 1.0;

        // Row 6
        derivative_mtx_ptr[index_shift + 30] = -4.0 * piGp * (degree_l_dbl + 1.0) * r_inv;
        derivative_mtx_ptr[index_shift + 31] = 4.0 * piGp * degree_l_dbl * lp1 * r_inv;
        derivative_mtx_ptr[index_shift + 32] = 0.0;
        derivative_mtx_ptr[index_shift + 33] = 0.0;
        derivative_mtx_ptr[index_shift + 34] = 0.0;
        derivative_mtx_ptr[index_shift + 35] = (degree_l_dbl - 1.0) * r_inv;
    }
}
