// boundaries_.hpp - Apply surface boundary conditions using Eigen
// Ported from TidalPy/RadialSolver/boundaries/boundaries.pyx
//
// References
// ----------
// KMN15: Kamata et al. (2015; JGR-P)
// S74: Saito (1974; JPE)
// KTC21: Kervazo et al. (2021; A&A)
#pragma once

#include <complex>
#include <limits>
#include <Eigen/Dense> // Replaced lapack_.hpp with Eigen

#include "../../constants_.hpp"


// Apply boundary conditions at the planet's surface by solving a linear system.
//
// Parameters
// ----------
// constant_vector_ptr : complex*, output
//     Constant vector (overwritten with solution).
// bc_solution_info_ptr : int*, output
//     Solver status (0 = success, >0 = error/singular).
// bc_pointer : double*
//     Boundary condition values.
// uppermost_y_per_solution_ptr : complex*
//     Y values at surface for each solution.
// surface_gravity : double
// G_to_use : double
// num_sols : size_t
// max_num_y : size_t
// ytype_i : size_t
// layer_type : int
//     0=solid, 1=liquid.
// layer_is_static : bool
// layer_is_incomp : bool
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
    *bc_solution_info_ptr = 0; // Default to success

    if (layer_type == 0)
    {
        // Solid layer (3x3 system)
        Eigen::Matrix3cd A;
        Eigen::Vector3cd B;

        // At the surface: y_2 = S_1; y_4 = S_4; y_6 = S_6 [See: B.37 in KTC21; 16 in KMN15]
        B(0) = bc_pointer[ytype_i * 3 + 0];
        B(1) = bc_pointer[ytype_i * 3 + 1];
        B(2) = bc_pointer[ytype_i * 3 + 2];

        // Fill A matrix (Eigen is column-major by default, just like Fortran/LAPACK)
        // We assign directly to (row, col) indices for clarity.
        A(0, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 1];
        A(1, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 3];
        A(2, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 5];

        A(0, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 1];
        A(1, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 3];
        A(2, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 5];

        A(0, 2) = uppermost_y_per_solution_ptr[2 * max_num_y + 1];
        A(1, 2) = uppermost_y_per_solution_ptr[2 * max_num_y + 3];
        A(2, 2) = uppermost_y_per_solution_ptr[2 * max_num_y + 5];

        // Solve A * X = B using partial-pivot LU (mirrors LAPACK zgesv).
        // Only fail if the solution is non-finite (truly singular).
        Eigen::PartialPivLU<Eigen::Matrix3cd> lu(A);
        Eigen::Vector3cd X = lu.solve(B);
        if (!X.allFinite()) {
            *bc_solution_info_ptr = 1; // Singular
            return;
        }

        constant_vector_ptr[0] = X(0);
        constant_vector_ptr[1] = X(1);
        constant_vector_ptr[2] = X(2);

    } else
    {
        if (layer_is_static)
        {
            // Static liquid layer: 1 solution, 1 BC (1x1 system)
            
            // y_7 = y_6 + (4 pi G / g) y_2
            std::complex<double> B_val = 
                bc_pointer[ytype_i * 3 + 2] +
                bc_pointer[ytype_i * 3 + 0] * (4.0 * TidalPyConstants::d_PI * G_to_use / surface_gravity);

            // y_7 held in index 1 (index 0 is y_5)
            std::complex<double> A_val = uppermost_y_per_solution_ptr[0 * max_num_y + 1];

            // 1D solve is just division
            if (std::abs(A_val) == 0.0) {
                *bc_solution_info_ptr = 1; // Singular
                return;
            }

            constant_vector_ptr[0] = B_val / A_val;
            constant_vector_ptr[1] = std::complex<double>(nan_val, nan_val);
            constant_vector_ptr[2] = std::complex<double>(nan_val, nan_val);

        } else
        {
            // Dynamic liquid layer: 2 solutions, 2 BCs (2x2 system)
            Eigen::Matrix2cd A;
            Eigen::Vector2cd B;

            B(0) = bc_pointer[ytype_i * 3 + 0];
            B(1) = bc_pointer[ytype_i * 3 + 2];

            // y_2 and y_6 at indices 1 and 3
            A(0, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 1];
            A(1, 0) = uppermost_y_per_solution_ptr[0 * max_num_y + 3];
            
            A(0, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 1];
            A(1, 1) = uppermost_y_per_solution_ptr[1 * max_num_y + 3];

            // Solve A * X = B using partial-pivot LU (mirrors LAPACK zgesv).
            Eigen::PartialPivLU<Eigen::Matrix2cd> lu(A);
            Eigen::Vector2cd X = lu.solve(B);
            if (!X.allFinite()) {
                *bc_solution_info_ptr = 1; // Singular
                return;
            }

            constant_vector_ptr[0] = X(0);
            constant_vector_ptr[1] = X(1);
            constant_vector_ptr[2] = std::complex<double>(nan_val, nan_val);
        }
    }
}


// Estimate how strongly the surface boundary condition solve amplifies error.
//
// The collapsed surface solution is y_k = sum_s c_s y_{k,s} over the independent solutions s. When the
// solution constants c_s are large and cancel (deep starting radii, high harmonic degrees), roundoff and
// integration error in the y_{k,s} are amplified into the collapsed values by roughly the returned factor:
// the largest cancellation scale sum_s |c_s| |y_{k,s}| across the y rows divided by the largest collapsed
// magnitude |y_k|. The relative accuracy of the surface solution, and of the Love numbers derived from it,
// is then floor limited to about (returned factor) * machine epsilon regardless of integration tolerance.
//
// Parameters
// ----------
// constant_vector_ptr : complex*
//     Solution constants from the surface boundary condition solve.
// uppermost_y_per_solution_ptr : complex*
//     Y values at the surface for each independent solution.
// num_sols : size_t
//     Number of independent solutions (1, 2, or 3).
// num_ys : size_t
//     Number of y values per solution.
// max_num_y : size_t
//     Stride between solutions in uppermost_y_per_solution_ptr.
//
// Returns
// -------
// double
//     Amplification factor (1 when a single solution leaves no room for cancellation).
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
        // Degenerate all-zero surface; no cancellation information.
        return 1.0;
    }
    if (max_collapsed_mag <= TidalPyConstants::d_EPS * max_cancellation_scale)
    {
        // Every y row cancels completely; cap at the largest meaningful amplification.
        return 1.0 / TidalPyConstants::d_EPS;
    }
    return max_cancellation_scale / max_collapsed_mag;
}
