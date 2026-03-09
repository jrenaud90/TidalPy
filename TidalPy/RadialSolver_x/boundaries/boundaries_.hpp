// boundaries_.hpp - Apply surface boundary conditions using LAPACK zgesv
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

#include "../../constants_.hpp"
#include "../lapack_.hpp"


// Apply boundary conditions at the planet's surface by solving a linear system.
//
// Parameters
// ----------
// constant_vector_ptr : complex*, output
//     Constant vector (overwritten with solution).
// bc_solution_info_ptr : int*, output
//     LAPACK solver status.
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

    // LAPACK parameters
    int lapack_nrhs = 1;
    int lapack_ipiv[10];
    int num_sols_int = static_cast<int>(num_sols);

    // Allocate surface matrices on stack
    std::complex<double> surface_matrix_solid[3][3];
    std::complex<double> surface_matrix_liquid_dynamic[2][2];
    std::complex<double> surface_matrix_liquid_static[1][1];
    std::complex<double>* surface_matrix_ptr;

    if (layer_type == 0) {
        // Solid layer
        surface_matrix_ptr = &surface_matrix_solid[0][0];

        // At the surface: y_2 = S_1; y_4 = S_4; y_6 = S_6 [See: B.37 in KTC21; 16 in KMN15]
        constant_vector_ptr[0] = bc_pointer[ytype_i * 3 + 0];
        constant_vector_ptr[1] = bc_pointer[ytype_i * 3 + 1];
        constant_vector_ptr[2] = bc_pointer[ytype_i * 3 + 2];

        // Transposed for FORTRAN-ordered zgesv
        surface_matrix_ptr[0] = uppermost_y_per_solution_ptr[0 * max_num_y + 1];
        surface_matrix_ptr[1] = uppermost_y_per_solution_ptr[0 * max_num_y + 3];
        surface_matrix_ptr[2] = uppermost_y_per_solution_ptr[0 * max_num_y + 5];
        surface_matrix_ptr[3] = uppermost_y_per_solution_ptr[1 * max_num_y + 1];
        surface_matrix_ptr[4] = uppermost_y_per_solution_ptr[1 * max_num_y + 3];
        surface_matrix_ptr[5] = uppermost_y_per_solution_ptr[1 * max_num_y + 5];
        surface_matrix_ptr[6] = uppermost_y_per_solution_ptr[2 * max_num_y + 1];
        surface_matrix_ptr[7] = uppermost_y_per_solution_ptr[2 * max_num_y + 3];
        surface_matrix_ptr[8] = uppermost_y_per_solution_ptr[2 * max_num_y + 5];
    } else {
        if (layer_is_static) {
            // Static liquid layer: 1 solution, 1 BC
            surface_matrix_ptr = &surface_matrix_liquid_dynamic[0][0];

            // y_7 = y_6 + (4 pi G / g) y_2
            constant_vector_ptr[0] =
                bc_pointer[ytype_i * 3 + 2] +
                bc_pointer[ytype_i * 3 + 0] * (4.0 * TidalPyConstants::d_PI * G_to_use / surface_gravity);

            constant_vector_ptr[1] = nan_val;
            constant_vector_ptr[2] = nan_val;

            // y_7 held in index 1 (index 0 is y_5)
            surface_matrix_ptr[0] = uppermost_y_per_solution_ptr[0 * max_num_y + 1];
        } else {
            // Dynamic liquid layer: 2 solutions, 2 BCs
            surface_matrix_ptr = &surface_matrix_liquid_static[0][0];

            constant_vector_ptr[0] = bc_pointer[ytype_i * 3 + 0];
            constant_vector_ptr[1] = bc_pointer[ytype_i * 3 + 2];

            constant_vector_ptr[2] = nan_val;

            // y_2 and y_6 at indices 1 and 3
            surface_matrix_ptr[0] = uppermost_y_per_solution_ptr[0 * max_num_y + 1];
            surface_matrix_ptr[1] = uppermost_y_per_solution_ptr[0 * max_num_y + 3];
            surface_matrix_ptr[2] = uppermost_y_per_solution_ptr[1 * max_num_y + 1];
            surface_matrix_ptr[3] = uppermost_y_per_solution_ptr[1 * max_num_y + 3];
        }
    }

    // Solve A * X = B using LAPACK zgesv
    zgesv_(
        &num_sols_int,         // N
        &lapack_nrhs,          // NRHS
        surface_matrix_ptr,    // A (overwritten)
        &num_sols_int,         // LDA
        lapack_ipiv,           // IPIV (output)
        constant_vector_ptr,   // B -> X (overwritten)
        &num_sols_int,         // LDB
        bc_solution_info_ptr   // INFO (output)
        );
}
