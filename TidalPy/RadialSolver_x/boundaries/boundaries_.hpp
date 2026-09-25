// boundaries_.hpp: surface boundary condition solve.
//
// References
// ----------
// KMN15: Kamata et al. (2015; JGR-P)
// S74: Saito (1974; JPE)
// KTC21: Kervazo et al. (2021; A&A)
#pragma once

#include <cmath>
#include <cstddef>
#include <complex>
#include <limits>
#include <Eigen/Dense>

#include "../../constants_.hpp"


// Solve the surface linear system for the collapse constants of the top layer (written to constant_vector_ptr,
// unused entries NaN). bc_solution_info_ptr is 0 on success, 1 when the system is singular. layer_type: 0 =
// solid (3x3 system), 1 = liquid (1x1 static, 2x2 dynamic).
inline void c_apply_surface_bc(
        std::complex<double>* constant_vector_ptr,
        int* bc_solution_info_ptr,
        double* bc_pointer,
        std::complex<double>* uppermost_y_per_solution_ptr,
        double surface_gravity,
        double G_to_use,
        size_t num_sols,
        size_t max_num_y,
        size_t ytype_i,
        int layer_type,
        bool layer_is_static,
        bool layer_is_incomp) noexcept
{
    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    *bc_solution_info_ptr = 0;

    if (layer_type == 0)
    {
        Eigen::Matrix3cd A;
        Eigen::Vector3cd B;

        // At the surface y_2 = S_1, y_4 = S_4, y_6 = S_6 (KTC21 Eq. B.37; KMN15 Eq. 16).
        B(0) = bc_pointer[ytype_i * 3 + 0];
        B(1) = bc_pointer[ytype_i * 3 + 1];
        B(2) = bc_pointer[ytype_i * 3 + 2];

        A(0, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 1];
        A(1, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 3];
        A(2, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 5];

        A(0, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 1];
        A(1, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 3];
        A(2, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 5];

        A(0, 2) = uppermost_y_per_solution_ptr[2 * max_num_y + 1];
        A(1, 2) = uppermost_y_per_solution_ptr[2 * max_num_y + 3];
        A(2, 2) = uppermost_y_per_solution_ptr[2 * max_num_y + 5];

        // Fails only when the solution is non-finite (truly singular).
        Eigen::PartialPivLU<Eigen::Matrix3cd> lu(A);
        Eigen::Vector3cd X = lu.solve(B);
        if (!X.allFinite()) {
            *bc_solution_info_ptr = 1;
            return;
        }

        constant_vector_ptr[0] = X(0);
        constant_vector_ptr[1] = X(1);
        constant_vector_ptr[2] = X(2);

    } else
    {
        if (layer_is_static)
        {
            // y_7 = y_6 + (4 pi G / g) y_2
            std::complex<double> B_val =
                bc_pointer[ytype_i * 3 + 2] +
                bc_pointer[ytype_i * 3 + 0] * (4.0 * TidalPyConstants::d_PI * G_to_use / surface_gravity);

            // y_7 is at index 1 (index 0 is y_5).
            std::complex<double> A_val = uppermost_y_per_solution_ptr[0 * max_num_y + 1];

            if (std::abs(A_val) == 0.0) {
                *bc_solution_info_ptr = 1;
                return;
            }

            constant_vector_ptr[0] = B_val / A_val;
            constant_vector_ptr[1] = std::complex<double>(nan_val, nan_val);
            constant_vector_ptr[2] = std::complex<double>(nan_val, nan_val);

        } else
        {
            Eigen::Matrix2cd A;
            Eigen::Vector2cd B;

            B(0) = bc_pointer[ytype_i * 3 + 0];
            B(1) = bc_pointer[ytype_i * 3 + 2];

            // y_2 and y_6 are at indices 1 and 3.
            A(0, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 1];
            A(1, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 3];

            A(0, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 1];
            A(1, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 3];

            Eigen::PartialPivLU<Eigen::Matrix2cd> lu(A);
            Eigen::Vector2cd X = lu.solve(B);
            if (!X.allFinite()) {
                *bc_solution_info_ptr = 1;
                return;
            }

            constant_vector_ptr[0] = X(0);
            constant_vector_ptr[1] = X(1);
            constant_vector_ptr[2] = std::complex<double>(nan_val, nan_val);
        }
    }
}


// Error amplification of the surface solve: the collapsed y_k = sum_s c_s y_{k,s} cancels when the constants
// are large, and roundoff in y_{k,s} is amplified by about max_k sum_s |c_s| |y_{k,s}| / max_k |y_k|, which
// floors the relative accuracy of the Love numbers at that factor times machine epsilon. max_num_y is the
// stride between solutions.
inline double c_estimate_surface_amplification(
        const std::complex<double>* constant_vector_ptr,
        const std::complex<double>* uppermost_y_per_solution_ptr,
        size_t num_sols,
        size_t num_ys,
        size_t max_num_y) noexcept
{
    double max_cancellation_scale = 0.0;
    double max_collapsed_mag      = 0.0;
    for (size_t y_i = 0; y_i < num_ys; ++y_i)
    {
        double cancellation_scale = 0.0;
        std::complex<double> collapsed_value(0.0, 0.0);
        for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
        {
            const std::complex<double> constant = constant_vector_ptr[solution_i];
            const std::complex<double> y_value  = uppermost_y_per_solution_ptr[solution_i * max_num_y + y_i];
            cancellation_scale += std::abs(constant) * std::abs(y_value);
            collapsed_value    += constant * y_value;
        }
        max_cancellation_scale = std::fmax(max_cancellation_scale, cancellation_scale);
        max_collapsed_mag      = std::fmax(max_collapsed_mag, std::abs(collapsed_value));
    }

    if (max_cancellation_scale <= 0.0)
    {
        // All-zero surface.
        return 1.0;
    }
    if (max_collapsed_mag <= TidalPyConstants::d_EPS * max_cancellation_scale)
    {
        // Complete cancellation; cap at the largest meaningful amplification.
        return 1.0 / TidalPyConstants::d_EPS;
    }
    return max_cancellation_scale / max_collapsed_mag;
}


// Reciprocal 1-norm condition number of an N x N complex matrix, 1 / (||A||_1 ||A^-1||_1), with the inverse
// formed explicitly (N <= 3 here). 0 when A is zero or its inverse is not finite.
template <int N>
inline double c_reciprocal_condition_1norm(const Eigen::Matrix<std::complex<double>, N, N>& matrix) noexcept
{
    const double matrix_norm = matrix.cwiseAbs().colwise().sum().maxCoeff();
    if (!(matrix_norm > 0.0))
    {
        return 0.0;
    }
    const Eigen::PartialPivLU<Eigen::Matrix<std::complex<double>, N, N>> lu(matrix);
    const Eigen::Matrix<std::complex<double>, N, N> inverse = lu.inverse();
    if (!inverse.allFinite())
    {
        return 0.0;
    }
    const double inverse_norm = inverse.cwiseAbs().colwise().sum().maxCoeff();
    if (!(inverse_norm > 0.0))
    {
        return 0.0;
    }
    return 1.0 / (matrix_norm * inverse_norm);
}


// Rank measure of the surface solve: the reciprocal 1-norm condition number of the matrix c_apply_surface_bc
// solves (its rows are y2, y4, y6 for a solid, y2, y6 for a dynamic liquid, y7 for a static liquid; one column
// per independent solution), after equilibration so the result does not depend on units or on how each starting
// solution was normalized:
//     1. each radial function is divided by its largest magnitude across the solutions, and
//     2. each solution is divided by its largest scaled radial function (all of its ys, not only the three the
//        boundary conditions constrain).
// Step 2 is what exposes a combination of solutions that is O(1) in size but carries no surface traction or
// potential gradient, such as the rigid translation of a static body at degree 1: no surface condition can fix
// its amplitude. Near 1 is well conditioned (at most 1); a value near machine epsilon means the solution
// constants are undetermined; roughly, the constants lose log10(1 / rcond) digits relative to the error in the
// integrated solutions. c_estimate_surface_amplification measures cancellation in the collapse instead and can
// read 1 for a singular system. max_num_y is the stride between solutions; num_ys is 2 * num_sols.
inline double c_estimate_surface_rcond(
        const std::complex<double>* uppermost_y_per_solution_ptr,
        size_t num_sols,
        size_t num_ys,
        size_t max_num_y,
        int layer_type,
        bool layer_is_static) noexcept
{
    // The rows the boundary conditions constrain, in the layer's own y storage order. Liquid storage is
    // (y1, y2, y5, y6) when dynamic, constraining y2 and y6 (indices 1 and 3), and (y5, y7) when static,
    // constraining y7 (index 1), so the leading entries of the solid list serve both.
    const size_t bc_rows[3] = {1, 3, 5};
    size_t num_bc_rows = 3;
    if (layer_type != 0)
    {
        num_bc_rows = layer_is_static ? 1 : 2;
    }
    if ((num_bc_rows != num_sols) || (num_ys > max_num_y) || (num_ys != 2 * num_sols))
    {
        return TidalPyConstants::d_NAN;
    }

    // Step 1: the largest magnitude of each radial function across the solutions.
    double y_scale[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    for (size_t y_i = 0; y_i < num_ys; ++y_i)
    {
        double largest = 0.0;
        for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
        {
            const double magnitude = std::abs(uppermost_y_per_solution_ptr[solution_i * max_num_y + y_i]);
            if (!std::isfinite(magnitude))
            {
                return 0.0;
            }
            largest = std::fmax(largest, magnitude);
        }
        // A radial function that is zero in every solution stays zero; it constrains nothing.
        y_scale[y_i] = (largest > 0.0) ? largest : 1.0;
    }

    // Step 2: the size of each solution in the scaled radial functions.
    double solution_scale[3] = {1.0, 1.0, 1.0};
    for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
    {
        double largest = 0.0;
        for (size_t y_i = 0; y_i < num_ys; ++y_i)
        {
            largest = std::fmax(
                largest, std::abs(uppermost_y_per_solution_ptr[solution_i * max_num_y + y_i]) / y_scale[y_i]);
        }
        if (!(largest > 0.0))
        {
            // An identically zero solution: the basis itself is rank deficient.
            return 0.0;
        }
        solution_scale[solution_i] = largest;
    }

    const auto scaled_entry = [&](size_t row_i, size_t solution_i) {
        const size_t y_i = bc_rows[row_i];
        return uppermost_y_per_solution_ptr[solution_i * max_num_y + y_i] /
            (y_scale[y_i] * solution_scale[solution_i]);
    };

    if (num_bc_rows == 3)
    {
        Eigen::Matrix3cd equilibrated;
        for (size_t row_i = 0; row_i < 3; ++row_i)
        {
            for (size_t solution_i = 0; solution_i < 3; ++solution_i)
            {
                equilibrated(row_i, solution_i) = scaled_entry(row_i, solution_i);
            }
        }
        return c_reciprocal_condition_1norm<3>(equilibrated);
    }
    else if (num_bc_rows == 2)
    {
        Eigen::Matrix2cd equilibrated;
        for (size_t row_i = 0; row_i < 2; ++row_i)
        {
            for (size_t solution_i = 0; solution_i < 2; ++solution_i)
            {
                equilibrated(row_i, solution_i) = scaled_entry(row_i, solution_i);
            }
        }
        return c_reciprocal_condition_1norm<2>(equilibrated);
    }
    // Static liquid: one entry, so the only question is whether it vanishes.
    const std::complex<double> entry = scaled_entry(0, 0);
    return (std::abs(entry) > 0.0) ? 1.0 : 0.0;
}
