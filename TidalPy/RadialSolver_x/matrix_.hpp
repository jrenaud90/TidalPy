// matrix_.hpp - Propagation of tidal solution using the fundamental matrix
// Ported from TidalPy/RadialSolver/matrix.pyx
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
#include <Eigen/Dense> // Replaced lapack_.hpp with Eigen

#include "../constants_.hpp"
#include "../Material_x/eos/eos_solution_.hpp"
#include "../constants_.hpp"
#include "rs_solution_.hpp"
#include "boundaries/surface_bc_.hpp"
#include "matrix_types/solid_matrix_.hpp"


/// Propagation matrix solver for radial tidal solutions.
///
/// Currently only supports solid, static, incompressible layers.
///
/// Parameters
/// ----------
/// solution_storage_ptr : c_RadialSolutionStorage*
///     Pointer to the solution storage (contains EOS data and output arrays).
/// frequency : double
///     Forcing frequency [rad s-1].
/// planet_bulk_density : double
///     Bulk density of the planet [kg m-3].
/// first_slice_index_by_layer_ptr : size_t*
///     Array of first radial slice index per layer.
/// num_slices_by_layer_ptr : size_t*
///     Array of number of slices per layer.
/// num_layers_for_slices : size_t
///     Number of elements in the above arrays.
/// num_bc_models : size_t
///     Number of boundary condition models (tidal, loading, free).
/// bc_models_ptr : int*
///     Array of boundary condition model IDs.
/// G_to_use : double
///     Gravitational constant.
/// degree_l : int
///     Harmonic degree.
/// starting_radius : double
///     Starting radius for solution (0 = auto-determine).
/// start_radius_tolerance : double
///     Tolerance for starting radius formula.
/// core_model : int
///     Core model type (0-4).
/// verbose : bool
///     Print status messages.
///
/// Returns
/// -------
/// int : error code (0 = success)
inline int c_matrix_propagate(
    c_RadialSolutionStorage* solution_storage_ptr,
    double frequency,
    double planet_bulk_density,
    size_t* first_slice_index_by_layer_ptr,
    size_t* num_slices_by_layer_ptr,
    size_t num_layers_for_slices,
    size_t num_bc_models,
    int* bc_models_ptr,
    double G_to_use,
    int degree_l,
    double starting_radius,
    double start_radius_tolerance,
    int core_model,
    bool verbose) noexcept
{
    // Setup
    const std::complex<double> cmplx_zero(0.0, 0.0);
    const std::complex<double> cmplx_one(1.0, 0.0);
    const std::complex<double> cmplx_NAN(
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN());

    // Get raw pointer of radial solver storage and eos storage
    c_EOSSolution* eos_solution_storage_ptr = solution_storage_ptr->get_eos_solution_ptr();

    solution_storage_ptr->message = "RadialSolver.PropMatrixMethod:: Propagator Matrix Method Called.\n";

    // Pull out key information
    const size_t num_layers     = eos_solution_storage_ptr->num_layers;
    const size_t total_slices   = eos_solution_storage_ptr->radius_array_size;
    const size_t top_slice_i    = total_slices - 1;

    // Alias pointers to EOS properties
    double* radius_array_ptr  = eos_solution_storage_ptr->radius_array_vec.data();
    double* gravity_array_ptr = eos_solution_storage_ptr->gravity_array_vec.data();
    double* density_array_ptr = eos_solution_storage_ptr->density_array_vec.data();

    // Need to recast the storage's shear/bulk double arrays to double complex for local use
    std::complex<double>* complex_shear_array_ptr =
        reinterpret_cast<std::complex<double>*>(eos_solution_storage_ptr->complex_shear_array_vec.data());

    // Pull out constants
    const double planet_radius = radius_array_ptr[top_slice_i];

    // Find boundary condition at the top of the planet
    const double degree_l_dbl      = static_cast<double>(degree_l);
    const size_t num_ytypes        = num_bc_models;
    const size_t num_output_ys     = C_MAX_NUM_Y * num_ytypes;

    // Boundary condition array: 15 = 5 (max models) * 3 (conditions per model)
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

    // Make conversions from TS72 to SV04
    // The propagation matrix uses different sign conventions so we need to convert
    for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
    {
        const size_t full_shift = 3 * ytype_i;
        if (bc_models_ptr[ytype_i] == 0)
        {
            // Free surface, no conversion needed.
        }
        else if (bc_models_ptr[ytype_i] == 1)
        {
            // Tidal boundary condition - last component is negative of TS74 (see Eq. 1.127 of S&V2004)
            bc_pointer[full_shift + 2] *= -1.0;
        }
        else if (bc_models_ptr[ytype_i] == 2)
        {
            // Loading boundary condition.
            bc_pointer[full_shift + 2] *= -1.0;
        }
    }

    // Determine starting radius slice
    if (starting_radius == 0.0)
    {
        // Use a model involving planet radius and degree l (based on H. Martens thesis and LoadDef manual).
        starting_radius = planet_radius * std::pow(start_radius_tolerance, 1.0 / degree_l_dbl);
        // Ensure not too close to the surface
        starting_radius = std::fmin(starting_radius, 0.95 * planet_radius);
    }

    // Determine which layer this starting radius resides in
    double layer_upper_radius, last_layer_upper_radius, radius_check;
    size_t start_layer_i           = 0;
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
            start_layer_i = current_layer_i;
            first_slice_index = first_slice_index_by_layer_ptr[current_layer_i];

            // Find the last radial slice before the starting radius
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

    // For the propagation matrix we have to start at index 2.
    // Index 0 corresponds to r=0 where many fundamental matrix elements are nan/inf.
    first_slice_index = last_index_before_start + 1;
    if (first_slice_index == 0 || first_slice_index == 1)
        first_slice_index = 2;

    // Define memory for fundamental matrices (heap allocated since size is runtime-dependent)
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

    // Only populate matrix values starting at first_slice_index - 1
    // TODO: Currently only solid, static, incompressible layers are supported for matrix propagation.
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

    // Initialize the base of the propagation matrix to the initial conditions
    // From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
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

    // Step through the planet's shells and build the propagation matrix
    Eigen::Matrix3cd surface_matrix;
    surface_matrix.setConstant(cmplx_NAN); // Initialize to NaN for debugging safety

    std::complex<double> temp_cmplx;
    std::complex<double> temp_matrix[18];

    // Create the downstream solution arrays
    std::complex<double> surface_solution[3];
    std::complex<double> bc_copy[3];

    // Initialize to NaN for debugging
    for (size_t j = 0; j < 18; ++j)
    {
        if (j < 3)
        {
            surface_solution[j] = cmplx_NAN;
            bc_copy[j] = cmplx_NAN;
        }
        temp_matrix[j] = cmplx_NAN;
    }

    for (size_t slice_i = first_slice_index; slice_i < total_slices; ++slice_i)
    {
        // Need to start the index for this radial slice
        index_shift_36 = slice_i * 36;
        const size_t last_index_shift_36 = (slice_i - 1) * 36;
        index_shift_18 = slice_i * 18;
        const size_t last_index_shift_18 = (slice_i - 1) * 18;

        // P_{i} = Y_{i} @ ( Y_{i-1}^{-1} @ P_{i-1} )

        // First matrix multiplication: A = Y_{i-1}^{-1} @ P_{i-1}
        for (size_t j = 0; j < 6; ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                temp_cmplx = std::complex<double>(0.0, 0.0);
                for (size_t jj = 0; jj < 6; ++jj)
                {
                    temp_cmplx += (
                        inverse_fundamental_mtx_ptr[last_index_shift_36 + j * 6 + jj] *
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

        // At the surface, extract the surface condition matrix into the Eigen Matrix (3x3 from rows [3, 4, 6])
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

    // Used to convert from SVC radial solutions to T&S format
    std::complex<double> ts_conversion[6];
    for (size_t i = 0; i < 6; ++i)
        ts_conversion[i] = cmplx_NAN;

    // Cast solution pointer from double to complex
    double* solution_dbl_ptr = solution_storage_ptr->full_solution_vec.data();
    std::complex<double>* solution_ptr = reinterpret_cast<std::complex<double>*>(solution_dbl_ptr);

    size_t ytype_i = 0;
    while (solution_storage_ptr->error_code == 0)
    {
        if (ytype_i == num_bc_models)
            break;

        // Set up the RHS vector B using Eigen
        Eigen::Vector3cd B_vec;
        for (size_t i = 0; i < 3; ++i)
        {
            B_vec(i) = std::complex<double>(bc_pointer[ytype_i * 3 + i], 0.0);
        }

        // Solve the linear equation U = S^-1 @ B using Eigen's FullPivLU
        Eigen::FullPivLU<Eigen::Matrix3cd> lu(surface_matrix);

        // Check for singularity/errors
        if (!lu.isInvertible())
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

        // Extract the solution and push back to bc_copy for downstream loop
        Eigen::Vector3cd X = lu.solve(B_vec);
        for (size_t i = 0; i < 3; ++i)
        {
            bc_copy[i] = X(i);
        }

        // Step through each radial step and apply the propagation matrix to the surface solution
        for (size_t slice_i = 0; slice_i < total_slices; ++slice_i)
        {
            index_shift_18              = slice_i * 18;
            const size_t ytype_shift    = ytype_i * C_MAX_NUM_Y;
            const size_t full_shift     = num_output_ys * slice_i + ytype_shift;

            if (slice_i < first_slice_index)
            {
                // Not in the solution region yet. NaN the results.
                for (size_t i = 0; i < 6; ++i)
                    solution_ptr[full_shift + i] = cmplx_NAN;
            }
            else
            {
                // Matrix multiplication: prop_matrix @ surface_solution
                for (size_t j = 0; j < 6; ++j)
                {
                    temp_cmplx = std::complex<double>(0.0, 0.0);
                    for (size_t jj = 0; jj < 3; ++jj)
                    {
                        temp_cmplx += (
                            propagation_mtx_ptr[index_shift_18 + j * 3 + jj] *
                            bc_copy[jj]);  // bc_copy now contains the solution X
                    }
                    solution_ptr[full_shift + j] = temp_cmplx;
                }

                // Convert from SVC16 to TS72 sign convention (B13 Eq. 7)
                ts_conversion[0] = solution_ptr[full_shift + 0];         // No Change
                ts_conversion[1] = solution_ptr[full_shift + 2];         // Flip y3 for y2
                ts_conversion[2] = solution_ptr[full_shift + 1];         // Flip y2 for y3
                ts_conversion[3] = solution_ptr[full_shift + 3];         // No Change
                ts_conversion[4] = -1.0 * solution_ptr[full_shift + 4];  // Change sign
                ts_conversion[5] = -1.0 * solution_ptr[full_shift + 5];  // Change sign

                // Store converted values back
                for (size_t i = 0; i < 6; ++i)
                    solution_ptr[full_shift + i] = ts_conversion[i];
            }
        }

        // Next y-type
        ++ytype_i;
    }

    // Update solution status and return
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
