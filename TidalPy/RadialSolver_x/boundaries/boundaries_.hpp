// boundaries_.hpp: surface boundary condition solve.
//
// References
// ----------
// KMN15: Kamata et al. (2015; JGR-P)
// S74: Saito (1974; JPE)
// KTC21: Kervazo et al. (2021; A&A)
#pragma once

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
