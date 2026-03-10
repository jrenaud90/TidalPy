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

        // Solve A * X = B
        Eigen::FullPivLU<Eigen::Matrix3cd> lu(A);
        if (!lu.isInvertible()) {
            *bc_solution_info_ptr = 1; // Singular
            return;
        }
        
        Eigen::Vector3cd X = lu.solve(B);
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

            // Solve A * X = B
            Eigen::FullPivLU<Eigen::Matrix2cd> lu(A);
            if (!lu.isInvertible()) {
                *bc_solution_info_ptr = 1; // Singular
                return;
            }

            Eigen::Vector2cd X = lu.solve(B);
            constant_vector_ptr[0] = X(0);
            constant_vector_ptr[1] = X(1);
            constant_vector_ptr[2] = std::complex<double>(nan_val, nan_val);
        }
    }
}
