// matrix_.hpp: propagation of the tidal solution with the fundamental matrix.
//
// References
// ----------
// SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
// HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
// ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
// B13   : Beuthe (2013, DOI: 10.1016/j.icarus.2012.11.020)
#pragma once

#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <limits>
#include <Eigen/Dense>

#include "../constants_.hpp"
#include "../Material_x/eos/eos_solution_.hpp"
#include "../constants_.hpp"
#include "rs_solution_.hpp"
#include "boundaries/surface_bc_.hpp"
#include "matrix_types/solid_matrix_.hpp"


/// Propagation matrix solver for radial tidal solutions; solid, static, incompressible layers only.
///
/// Parameters
/// ----------
/// solution_storage_ptr : c_RadialSolutionStorage*
///     Holds the EOS data and receives the gridded solution.
/// frequency : double
///     Forcing frequency [rad s-1].
/// planet_bulk_density : double
///     [kg m-3].
/// slices_per_layer : size_t
///     Slices this call lays down in each layer. The method propagates from one slice to the next, so it is
///     discretized by construction; this is the only grid left in the radial solver.
/// num_bc_models, bc_models_ptr
///     Boundary condition models (free = 0, tidal = 1, loading = 2).
/// G_to_use : double
///     Gravitational constant in solve units.
/// degree_l : int
///     Harmonic degree.
/// starting_radius : double
///     0 selects the automatic choice governed by start_radius_tolerance.
/// core_model : int
///     Core starting condition (0 to 4).
/// verbose : bool
///     Print status messages.
///
/// Returns
/// -------
/// int
///     Error code, 0 on success.
inline int c_matrix_propagate(
    c_RadialSolutionStorage* solution_storage_ptr,
    double frequency,
    double planet_bulk_density,
    size_t slices_per_layer,
    size_t num_bc_models,
    int* bc_models_ptr,
    double G_to_use,
    int degree_l,
    double starting_radius,
    double start_radius_tolerance,
    int core_model,
    bool verbose) noexcept
{
    const std::complex<double> cmplx_zero(0.0, 0.0);
    const std::complex<double> cmplx_one(1.0, 0.0);
    const std::complex<double> cmplx_NAN(
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN());

    c_EOSSolution* eos_solution_storage_ptr = solution_storage_ptr->get_eos_solution_ptr();

    solution_storage_ptr->message = "RadialSolver.PropMatrixMethod:: Propagator Matrix Method Called.\n";

    // The matrix method fills the grid, so a prior shooting solve's interpolant state must not survive on
    // reused storage or get_radial_solution would dispatch to the wrong branch.
    solution_storage_ptr->reset_interpolant_storage();

    const size_t num_layers     = eos_solution_storage_ptr->num_layers;
    const size_t total_slices   = num_layers * slices_per_layer;
    const size_t top_slice_i    = total_slices - 1;

    if (num_layers == 0 || slices_per_layer < 5)
    {
        solution_storage_ptr->message = "RadialSolver.PropMatrixMethod:: at least 5 slices per layer are required.";
        solution_storage_ptr->error_code = -5;
        solution_storage_ptr->success    = false;
        return solution_storage_ptr->error_code;
    }

    // The method propagates from one slice to the next, so it needs a grid; built here.
    std::vector<double>& radius_grid = solution_storage_ptr->p_matrix_radius_solve;
    radius_grid.assign(total_slices, 0.0);
    std::vector<double> gravity_grid(total_slices, 0.0);
    std::vector<double> density_grid(total_slices, 0.0);
    std::vector<std::complex<double>> shear_grid(total_slices, cmplx_zero);
    std::vector<size_t> first_slice_index_by_layer(num_layers, 0);
    std::vector<size_t> num_slices_by_layer(num_layers, slices_per_layer);
    {
        c_EOSMaterialState eos_material_state;
        const double last_span = static_cast<double>(slices_per_layer - 1);
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
        {
            const double radius_top = eos_solution_storage_ptr->upper_radius_bylayer_vec[layer_i];
            const double radius_bot = (layer_i == 0)
                ? 0.0
                : eos_solution_storage_ptr->upper_radius_bylayer_vec[layer_i - 1];
            first_slice_index_by_layer[layer_i] = layer_i * slices_per_layer;
            for (size_t slice_j = 0; slice_j < slices_per_layer; ++slice_j)
            {
                const size_t slice_i = first_slice_index_by_layer[layer_i] + slice_j;
                const double radius_here =
                    radius_bot + (static_cast<double>(slice_j) / last_span) * (radius_top - radius_bot);
                radius_grid[slice_i] = radius_here;
                eos_solution_storage_ptr->call_material(layer_i, radius_here, eos_material_state);
                gravity_grid[slice_i] = eos_material_state.gravity;
                density_grid[slice_i] = eos_material_state.density;
                shear_grid[slice_i]   = eos_material_state.shear_modulus;
            }
        }
    }
    // The y-grid follows the grid this call chose.
    solution_storage_ptr->num_slices = total_slices;
    solution_storage_ptr->total_size =
        static_cast<size_t>(C_MAX_NUM_Y_REAL) * total_slices * solution_storage_ptr->num_ytypes;
    solution_storage_ptr->full_solution_vec.resize(solution_storage_ptr->total_size);

    const size_t* first_slice_index_by_layer_ptr = first_slice_index_by_layer.data();
    const size_t* num_slices_by_layer_ptr        = num_slices_by_layer.data();
    double* radius_array_ptr  = radius_grid.data();
    double* gravity_array_ptr = gravity_grid.data();
    double* density_array_ptr = density_grid.data();
    std::complex<double>* complex_shear_array_ptr = shear_grid.data();

    const double planet_radius = radius_array_ptr[top_slice_i];

    const double degree_l_dbl      = static_cast<double>(degree_l);
    const size_t num_ytypes        = num_bc_models;
    const size_t num_output_ys     = C_MAX_NUM_Y * num_ytypes;

    // 15 = 5 (max solve_for entries) * 3 (surface conditions)
    double boundary_conditions[15];
    double* bc_pointer = &boundary_conditions[0];
    int bc_error = c_get_surface_bc(
        bc_pointer,
        bc_models_ptr,
        num_ytypes,
        planet_radius,
        planet_bulk_density,
        degree_l_dbl);

    if (bc_error != 0)
    {
        solution_storage_ptr->message = "RadialSolver.PropMatrixMethod:: Error computing surface boundary conditions.";
        solution_storage_ptr->error_code = bc_error;
        solution_storage_ptr->success = false;
        return bc_error;
    }

    // TS72 to SVC16 sign convention: the last component flips for the tidal and loading cases (SVC16
    // Eq. 1.127); a free surface needs no change.
    for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
    {
        const size_t full_shift = 3 * ytype_i;
        if (bc_models_ptr[ytype_i] == 0)
        {
        }
        else if (bc_models_ptr[ytype_i] == 1)
        {
            bc_pointer[full_shift + 2] *= -1.0;
        }
        else if (bc_models_ptr[ytype_i] == 2)
        {
            bc_pointer[full_shift + 2] *= -1.0;
        }
    }

    // Automatic starting radius after Martens' thesis and the LoadDef manual, capped by the config.
    if (starting_radius == 0.0)
    {
        starting_radius = planet_radius * std::pow(start_radius_tolerance, 1.0 / degree_l_dbl);
        starting_radius = std::fmin(
            starting_radius, tidalpy_config_ptr->d_MAX_START_RADIUS_FRAC * planet_radius);
    }

    // Find the layer holding the starting radius.
    double layer_upper_radius, last_layer_upper_radius, radius_check;
    size_t last_index_before_start = 0;
    size_t first_slice_index       = 0;
    double last_radius_check       = 0.0;

    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        layer_upper_radius = eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i];
        if (current_layer_i == 0)
            last_layer_upper_radius = 0.0;
        else
            last_layer_upper_radius = eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i - 1];

        if (last_layer_upper_radius < starting_radius && starting_radius <= layer_upper_radius)
        {
            first_slice_index = first_slice_index_by_layer_ptr[current_layer_i];

            for (size_t slice_i = first_slice_index;
                 slice_i < first_slice_index + num_slices_by_layer_ptr[current_layer_i];
                 ++slice_i)
            {
                radius_check = radius_array_ptr[slice_i];
                if (last_radius_check < starting_radius && starting_radius <= radius_check)
                {
                    if (slice_i == 0)
                        last_index_before_start = 0;
                    else
                        last_index_before_start = slice_i - 1;
                    break;
                }
                else
                {
                    last_radius_check = radius_check;
                }
            }
            break;
        }
        else
        {
            last_radius_check = last_layer_upper_radius;
        }
    }

    // Start at index 2 at the earliest: index 0 is r = 0, where fundamental matrix elements are NaN or inf.
    first_slice_index = last_index_before_start + 1;
    if (first_slice_index == 0 || first_slice_index == 1)
        first_slice_index = 2;

    const size_t matrix_size   = 6 * 6 * total_slices;
    const size_t prop_mat_size = 6 * 3 * total_slices;

    std::vector<std::complex<double>> fundamental_mtx_vec(matrix_size);
    std::vector<std::complex<double>> inverse_fundamental_mtx_vec(matrix_size);
    std::vector<std::complex<double>> derivative_mtx_vec(matrix_size);
    std::vector<std::complex<double>> propagation_mtx_vec(prop_mat_size);

    std::complex<double>* fundamental_mtx_ptr         = fundamental_mtx_vec.data();
    std::complex<double>* inverse_fundamental_mtx_ptr = inverse_fundamental_mtx_vec.data();
    std::complex<double>* derivative_mtx_ptr          = derivative_mtx_vec.data();
    std::complex<double>* propagation_mtx_ptr         = propagation_mtx_vec.data();

    // Matrices are filled from first_slice_index - 1 upward. TODO: only solid, static, incompressible layers.
    c_fundamental_matrix(
        first_slice_index - 1,
        total_slices,
        radius_array_ptr,
        density_array_ptr,
        gravity_array_ptr,
        complex_shear_array_ptr,
        fundamental_mtx_ptr,
        inverse_fundamental_mtx_ptr,
        derivative_mtx_ptr,
        degree_l,
        G_to_use);

    // Each slice i is a shell (r_{i-1}, r_i] of its own material, and inside it y = Y_i(r) c_i. Continuity at
    // r_{i-1} gives c_i = Y_i(r_{i-1})^-1 y(r_{i-1}), so the propagator needs the inverse of slice i's matrix at the
    // shell's lower radius: slice i's density and shear modulus with the radius and gravity of slice i - 1. Using
    // slice i - 1's own inverse instead would cancel every interior factor and leave a uniform body.
    std::vector<double> radius_lower_grid(total_slices, 0.0);
    std::vector<double> gravity_lower_grid(total_slices, 0.0);
    for (size_t slice_i = 1; slice_i < total_slices; ++slice_i)
    {
        radius_lower_grid[slice_i]  = radius_array_ptr[slice_i - 1];
        gravity_lower_grid[slice_i] = gravity_array_ptr[slice_i - 1];
    }
    std::vector<std::complex<double>> fundamental_lower_mtx_vec(matrix_size);
    std::vector<std::complex<double>> inverse_lower_mtx_vec(matrix_size);
    c_fundamental_matrix(
        first_slice_index,
        total_slices,
        radius_lower_grid.data(),
        density_array_ptr,
        gravity_lower_grid.data(),
        complex_shear_array_ptr,
        fundamental_lower_mtx_vec.data(),
        inverse_lower_mtx_vec.data(),
        derivative_mtx_ptr,
        degree_l,
        G_to_use);
    const std::complex<double>* inverse_lower_mtx_ptr = inverse_lower_mtx_vec.data();

    // Seed with the core starting conditions. From IcyDwarf: "They are inconsequential on the rest of the
    // solution, so false assumptions are OK."
    size_t index_shift_18 = (first_slice_index - 1) * 18;
    size_t index_shift_36 = (first_slice_index - 1) * 36;

    if (core_model == 0)
    {
        // Henning & Hurford (2014): seed matrix from first three columns of Y at base layer
        for (size_t j = 0; j < 6; ++j)
        {
            const size_t row_shift_index = index_shift_18 + (j * 3);
            for (size_t k = 0; k < 3; ++k)
                propagation_mtx_ptr[row_shift_index + k] = fundamental_mtx_ptr[index_shift_36 + j * 6 + k];
        }
    }
    else if (core_model == 1)
    {
        // Roberts & Nimmo (2008): liquid innermost zone
        for (size_t j = 0; j < 6; ++j)
        {
            const size_t row_shift_index = index_shift_18 + (j * 3);
            for (size_t k = 0; k < 3; ++k)
            {
                if ((j == 2) && (k == 0))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else if ((j == 3) && (k == 1))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else if ((j == 5) && (k == 2))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_zero;
            }
        }
    }
    else if (core_model == 2)
    {
        // Solid Inner Core (Based on Henning & Hurford 2014)
        for (size_t j = 0; j < 6; ++j)
        {
            const size_t row_shift_index = index_shift_18 + (j * 3);
            for (size_t k = 0; k < 3; ++k)
            {
                if ((j == 0) && (k == 0))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else if ((j == 1) && (k == 1))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else if ((j == 2) && (k == 2))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_zero;
            }
        }
    }
    else if (core_model == 3)
    {
        // Liquid Inner Core (based on Tobie+2005; as determined by Marc Neveu for IcyDwarf)
        for (size_t j = 0; j < 6; ++j)
        {
            const size_t row_shift_index = index_shift_18 + (j * 3);
            for (size_t k = 0; k < 3; ++k)
            {
                if ((j == 0) && (k == 0))
                    propagation_mtx_ptr[row_shift_index + k] = std::complex<double>(0.05, 0.0);
                else if ((j == 1) && (k == 1))
                    propagation_mtx_ptr[row_shift_index + k] = std::complex<double>(0.01, 0.0);
                else if ((j == 5) && (k == 2))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_zero;
            }
        }
    }
    else if (core_model == 4)
    {
        // Interface matrix from SVC Eq. 1.150
        const double grav_constant = (4.0 / 3.0) * TidalPyConstants::d_PI * G_to_use * density_array_ptr[first_slice_index - 1];
        for (size_t j = 0; j < 6; ++j)
        {
            const size_t row_shift_index = index_shift_18 + (j * 3);
            for (size_t k = 0; k < 3; ++k)
            {
                if ((j == 0) && (k == 0))
                    propagation_mtx_ptr[row_shift_index + k] = -std::pow(radius_array_ptr[first_slice_index - 1], degree_l - 1) / grav_constant;
                else if ((j == 0) && (k == 2))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else if ((j == 1) && (k == 1))
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one;
                else if ((j == 2) && (k == 2))
                    propagation_mtx_ptr[row_shift_index + k] = density_array_ptr[first_slice_index - 1] * grav_constant * radius_array_ptr[first_slice_index - 1];
                else if ((j == 4) && (k == 0))
                    propagation_mtx_ptr[row_shift_index + k] = std::pow(radius_array_ptr[first_slice_index - 1], degree_l);
                else if ((j == 5) && (k == 0))
                    propagation_mtx_ptr[row_shift_index + k] = 2.0 * (degree_l - 1) * std::pow(radius_array_ptr[first_slice_index - 1], degree_l - 1);
                else if ((j == 5) && (k == 2))
                    propagation_mtx_ptr[row_shift_index + k] = 3.0 * grav_constant;
                else
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_zero;
            }
        }
    }
    else
    {
        solution_storage_ptr->message =
            "RadialSolver.PropMatrixMethod:: Unknown starting core conditions encountered in `c_matrix_propagate`: " +
            std::to_string(core_model) + " (acceptable values: 0, 1, 2, 3, 4)\n";
        solution_storage_ptr->error_code = -20;
        solution_storage_ptr->success = false;
        if (verbose)
            std::printf("%s", solution_storage_ptr->message.c_str());
        return solution_storage_ptr->error_code;
    }

    Eigen::Matrix3cd surface_matrix;
    surface_matrix.setConstant(cmplx_NAN);

    std::complex<double> temp_cmplx;
    std::complex<double> temp_matrix[18];

    std::complex<double> surface_solution[3];
    std::complex<double> bc_copy[3];

    // NaN so uninitialized reads are visible.
    for (size_t j = 0; j < 18; ++j)
    {
        if (j < 3)
        {
            surface_solution[j] = cmplx_NAN;
            bc_copy[j] = cmplx_NAN;
        }
        temp_matrix[j] = cmplx_NAN;
    }

    // Build the propagation matrix shell by shell.
    for (size_t slice_i = first_slice_index; slice_i < total_slices; ++slice_i)
    {
        index_shift_36 = slice_i * 36;
        index_shift_18 = slice_i * 18;
        const size_t last_index_shift_18 = (slice_i - 1) * 18;

        // P_{i} = Y_{i}(r_i) @ ( Y_{i}(r_{i-1})^{-1} @ P_{i-1} )

        // First matrix multiplication: A = Y_{i}(r_{i-1})^{-1} @ P_{i-1}
        for (size_t j = 0; j < 6; ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                temp_cmplx = std::complex<double>(0.0, 0.0);
                for (size_t jj = 0; jj < 6; ++jj)
                {
                    temp_cmplx += (
                        inverse_lower_mtx_ptr[index_shift_36 + j * 6 + jj] *
                        propagation_mtx_ptr[last_index_shift_18 + jj * 3 + k]);
                }
                temp_matrix[j * 3 + k] = temp_cmplx;
            }
        }

        // Outer matrix multiplication: P_{i} = Y_{i} @ A
        for (size_t j = 0; j < 6; ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                temp_cmplx = std::complex<double>(0.0, 0.0);
                for (size_t jj = 0; jj < 6; ++jj)
                {
                    temp_cmplx += (
                        fundamental_mtx_ptr[index_shift_36 + j * 6 + jj] *
                        temp_matrix[jj * 3 + k]);
                }
                propagation_mtx_ptr[index_shift_18 + j * 3 + k] = temp_cmplx;
            }
        }

        // Surface condition matrix: rows 3, 4, 6 of the propagation matrix.
        if (slice_i == (total_slices - 1))
        {
            for (size_t i = 0; i < 3; ++i)
            {
                surface_matrix(0, i) = propagation_mtx_ptr[index_shift_18 + (2 * 3) + i];
                surface_matrix(1, i) = propagation_mtx_ptr[index_shift_18 + (3 * 3) + i];
                surface_matrix(2, i) = propagation_mtx_ptr[index_shift_18 + (5 * 3) + i];
            }
        }
    }

    std::complex<double> ts_conversion[6];
    for (size_t i = 0; i < 6; ++i)
        ts_conversion[i] = cmplx_NAN;

    double* solution_dbl_ptr = solution_storage_ptr->full_solution_vec.data();
    std::complex<double>* solution_ptr = reinterpret_cast<std::complex<double>*>(solution_dbl_ptr);

    size_t ytype_i = 0;
    while (solution_storage_ptr->error_code == 0)
    {
        if (ytype_i == num_bc_models)
            break;

        Eigen::Vector3cd B_vec;
        for (size_t i = 0; i < 3; ++i)
        {
            B_vec(i) = std::complex<double>(bc_pointer[ytype_i * 3 + i], 0.0);
        }

        // Solve U = S^-1 B. FullPivLU::isInvertible's relative threshold called a badly scaled SI matrix singular,
        // so a solve is judged by whether its answer is finite, as the shooting method's surface solve does.
        Eigen::PartialPivLU<Eigen::Matrix3cd> lu(surface_matrix);
        Eigen::Vector3cd X = lu.solve(B_vec);

        if (!X.allFinite())
        {
            solution_storage_ptr->message =
                "RadialSolver.PropMatrixMethod:: Error encountered while applying surface boundary condition.\n"
                "Eigen FullPivLU: Surface matrix is singular or poorly conditioned.\n"
                "The solutions may not be valid at the surface.\n";
            solution_storage_ptr->error_code = -21;
            solution_storage_ptr->success = false;
            if (verbose)
                std::printf("%s", solution_storage_ptr->message.c_str());
            return solution_storage_ptr->error_code;
        }

        for (size_t i = 0; i < 3; ++i)
        {
            bc_copy[i] = X(i);
        }

        // Apply the propagation matrix to the surface solution at every slice.
        for (size_t slice_i = 0; slice_i < total_slices; ++slice_i)
        {
            index_shift_18              = slice_i * 18;
            const size_t ytype_shift    = ytype_i * C_MAX_NUM_Y;
            const size_t full_shift     = num_output_ys * slice_i + ytype_shift;

            if (slice_i < first_slice_index)
            {
                for (size_t i = 0; i < 6; ++i)
                    solution_ptr[full_shift + i] = cmplx_NAN;
            }
            else
            {
                for (size_t j = 0; j < 6; ++j)
                {
                    temp_cmplx = std::complex<double>(0.0, 0.0);
                    for (size_t jj = 0; jj < 3; ++jj)
                    {
                        temp_cmplx += (
                            propagation_mtx_ptr[index_shift_18 + j * 3 + jj] *
                            bc_copy[jj]);
                    }
                    solution_ptr[full_shift + j] = temp_cmplx;
                }

                // SVC16 to TS72 convention (B13 Eq. 7): swap y2 and y3, negate y5 and y6.
                ts_conversion[0] = solution_ptr[full_shift + 0];
                ts_conversion[1] = solution_ptr[full_shift + 2];
                ts_conversion[2] = solution_ptr[full_shift + 1];
                ts_conversion[3] = solution_ptr[full_shift + 3];
                ts_conversion[4] = -1.0 * solution_ptr[full_shift + 4];
                ts_conversion[5] = -1.0 * solution_ptr[full_shift + 5];

                for (size_t i = 0; i < 6; ++i)
                    solution_ptr[full_shift + i] = ts_conversion[i];
            }
        }

        ++ytype_i;
    }

    if (solution_storage_ptr->error_code != 0)
    {
        solution_storage_ptr->success = false;
    }
    else
    {
        solution_storage_ptr->success = true;
        solution_storage_ptr->message = "RadialSolver.MatrixPropagation: Completed without any noted issues.";
    }

    return solution_storage_ptr->error_code;
}
