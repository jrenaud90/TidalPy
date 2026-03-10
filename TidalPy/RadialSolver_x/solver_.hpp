// solver_.hpp - Top-level radial solver declarations
// Ported from TidalPy/RadialSolver/solver.pyx
//
// NOTE: The actual solver (cf_radial_solver) and Python wrapper (radial_solver)
// live in solver.pyx because they depend on:
//   - CyRK's Cython-only API (ODEMethod, PreEvalFunc)
//   - Material_x EOS solver (c_solve_eos via Cython)
//   - The shooting method (cf_shooting_solver in shooting.pyx)
//
// This header exists as a placeholder for any future C++ helper functions.
#pragma once

#include "constants_.hpp"
#include "rs_solution_.hpp"

// CyRK imports
#include "cysolution.hpp"
#include "c_events.hpp"
#include "cysolve.hpp"

// TidalPy imports
#include "../constants_.hpp"
#include "../utilities/math/numerics_.hpp"
#include "../utilities/dimensions/nondimensional_.hpp"

// RadialSolver imports
#include "rs_constants_.hpp"
#include "rs_solution_.hpp"
#include "love_.hpp"
#include "shooting_.hpp"
#include "matrix_.hpp"
#include "derivatives/odes_.hpp"

// Material imports
#include "../Material_x/eos/eos_solution_.hpp"
#include "../Material_x/eos/solver_.hpp"
#include "../Material_x/eos/methods/interpolate_.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <complex>
#include <stdexcept>
#include <cctype>


constexpr int C_EOS_INTERPOLATE_METHOD_INT = 0;

int c_radial_solver(
    c_RadialSolutionStorage* solution_storage_ptr,
    size_t total_slices,
    double* radius_array_in_ptr,
    double* density_array_in_ptr,
    std::complex<double>* complex_bulk_modulus_in_ptr,
    std::complex<double>* complex_shear_modulus_in_ptr,
    double frequency,
    double planet_bulk_density,
    size_t num_layers,
    int* layer_types_ptr,
    bool* is_static_bylayer_ptr,
    bool* is_incompressible_bylayer_ptr,
    double surface_pressure,
    int degree_l,
    size_t num_bc_models,
    int* bc_models_ptr,
    int core_model,
    bool use_kamata,
    double starting_radius,
    double start_radius_tolerance,
    ODEMethod integration_method_int,
    double integration_rtol,
    double integration_atol,
    bool scale_rtols_bylayer_type,
    size_t max_num_steps,
    size_t expected_size,
    size_t max_ram_MB,
    double max_step,
    bool nondimensionalize,
    bool use_prop_matrix,
    int* eos_integration_method_int_bylayer_ptr,
    ODEMethod eos_integration_method,
    double eos_rtol,
    double eos_atol,
    double eos_pressure_tol,
    int eos_max_iters,
    bool verbose,
    bool warnings
) noexcept
{
    // Check if we have configs.
    if (tidalpy_config_ptr == nullptr)
    {
        solution_storage_ptr->error_code = -999;
        solution_storage_ptr->message = std::string(
            "RadialSolver_x:: Fatal Error. TidalPyConfig pointer is uninitialized. "
            "initialize_tidalpy_config() must be called before the solver runs.\n"
        );
        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }
        return solution_storage_ptr->error_code;
    }

    // Figure out how many slices are in each layer
    std::vector<size_t> first_slice_index_by_layer_vec(num_layers);
    std::vector<size_t> num_slices_by_layer_vec(num_layers);

    // Pull out raw pointer
    c_EOSSolution* eos_solution_storage_ptr = solution_storage_ptr->get_eos_solution_ptr();

    // EOS variables
    size_t bottom_slice_index;
    std::vector<PreEvalFunc> eos_function_bylayer_vec(num_layers);
    c_EOS_ODEInput eos_input;
    std::vector<c_EOS_ODEInput> eos_inputs_bylayer_vec;
    eos_inputs_bylayer_vec.reserve(num_layers);
    std::vector<c_InterpolateEOSInput> specific_eos_input_bylayer_vec;
    specific_eos_input_bylayer_vec.reserve(num_layers);
    char* specific_eos_char_ptr = nullptr;

    // Ensure there is at least one layer
    if (num_layers == 0)
    {
        solution_storage_ptr->error_code = -5;
        solution_storage_ptr->message    =
            std::string("RadialSolver_x:: requires at least one layer, zero provided.\n");
        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }
        return solution_storage_ptr->error_code;
    }

    if (solution_storage_ptr->error_code == 0)
    {
        bool top_layer = false;
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
        {
            if (layer_i == num_layers - 1)
            {
                top_layer = true;
            }

            // Determine starting slice index
            if (layer_i == 0)
            {
                first_slice_index_by_layer_vec[layer_i] = 0;
            }
            else
            {
                first_slice_index_by_layer_vec[layer_i] = first_slice_index_by_layer_vec[layer_i - 1] + num_slices_by_layer_vec[layer_i - 1];
            }

            const double layer_upper_radius = eos_solution_storage_ptr->upper_radius_bylayer_vec[layer_i];

            size_t layer_slices    = 0;
            size_t interface_check = 0;
            for (size_t slice_i = first_slice_index_by_layer_vec[layer_i]; slice_i < total_slices; ++slice_i)
            {
                const double radius_check = radius_array_in_ptr[slice_i];

                if (c_isclose(radius_check, layer_upper_radius, 1.0e-9, 0.0))
                {
                    interface_check += 1;
                    if (interface_check > 1)
                    {
                        break;
                    }
                }
                else if (radius_check > layer_upper_radius)
                {
                    break;
                }
                layer_slices += 1;
            }

            if (layer_slices < 5)
            {
                solution_storage_ptr->error_code = -5;
                solution_storage_ptr->message    = std::string("RadialSolver_x:: At least five layer slices per layer are required.\n");
                if (verbose)
                {
                    printf("%s", solution_storage_ptr->message.c_str());
                }
                return solution_storage_ptr->error_code;
            }

            num_slices_by_layer_vec[layer_i] = layer_slices;
        }
    }

    // Get other needed inputs
    const double radius_planet = radius_array_in_ptr[total_slices - 1];

    double G_to_use                = tidalpy_config_ptr->d_G;
    double radius_planet_to_use    = radius_planet;
    double bulk_density_to_use     = planet_bulk_density;
    double frequency_to_use        = frequency;
    double surface_pressure_to_use = surface_pressure;
    double starting_radius_to_use  = starting_radius;

    c_NonDimensionalScales non_dim_scales(
        frequency,
        radius_planet,
        planet_bulk_density
    );

    if (nondimensionalize && solution_storage_ptr->error_code == 0)
    {
        for (size_t slice_i = 0; slice_i < total_slices; ++slice_i)
        {
            radius_array_in_ptr[slice_i]          /= non_dim_scales.length_conversion;
            density_array_in_ptr[slice_i]         /= non_dim_scales.density_conversion;
            complex_bulk_modulus_in_ptr[slice_i]  /= non_dim_scales.pascal_conversion;
            complex_shear_modulus_in_ptr[slice_i] /= non_dim_scales.pascal_conversion;
        }

        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
        {
            eos_solution_storage_ptr->upper_radius_bylayer_vec[layer_i] /= non_dim_scales.length_conversion;
        }

        G_to_use                = tidalpy_config_ptr->d_G / (non_dim_scales.length3_conversion / 
            (non_dim_scales.mass_conversion * non_dim_scales.second2_conversion));
        radius_planet_to_use    = radius_planet / non_dim_scales.length_conversion;
        bulk_density_to_use     = planet_bulk_density / non_dim_scales.density_conversion;
        frequency_to_use        = frequency / (1.0 / non_dim_scales.second_conversion);
        surface_pressure_to_use = surface_pressure / non_dim_scales.pascal_conversion;
        starting_radius_to_use  = starting_radius / non_dim_scales.length_conversion;

        solution_storage_ptr->change_radius_array(radius_array_in_ptr, total_slices, true);
    }

    // Solve the equation of state
    if (solution_storage_ptr->error_code == 0)
    {
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
        {
            if (eos_integration_method_int_bylayer_ptr[layer_i] == C_EOS_INTERPOLATE_METHOD_INT)
            {
                eos_function_bylayer_vec[layer_i] = c_preeval_interpolate;
                bottom_slice_index                = first_slice_index_by_layer_vec[layer_i];

                specific_eos_input_bylayer_vec.emplace_back(
                    num_slices_by_layer_vec[layer_i],
                    &radius_array_in_ptr[bottom_slice_index],
                    &density_array_in_ptr[bottom_slice_index],
                    &complex_bulk_modulus_in_ptr[bottom_slice_index],
                    &complex_shear_modulus_in_ptr[bottom_slice_index]
                );
                specific_eos_char_ptr = reinterpret_cast<char*>(&specific_eos_input_bylayer_vec.back());

                eos_inputs_bylayer_vec.emplace_back(
                    G_to_use,
                    radius_planet_to_use,
                    specific_eos_char_ptr,
                    false,
                    false,
                    false
                );
            }
            else
            {
                solution_storage_ptr->error_code = -250;
                break;
            }
        }
    }

    if (solution_storage_ptr->error_code == 0)
    {
        c_solve_eos(
            eos_solution_storage_ptr,
            eos_function_bylayer_vec,
            eos_inputs_bylayer_vec,
            bulk_density_to_use,
            surface_pressure_to_use,
            G_to_use,
            eos_integration_method,
            eos_rtol,
            eos_atol,
            eos_pressure_tol,
            eos_max_iters,
            verbose
        );
    }

    // Run requested radial solver method
    int sub_process_error_code = 0;
    if (eos_solution_storage_ptr->success && solution_storage_ptr->error_code == 0)
    {
        if (use_prop_matrix)
        {
            sub_process_error_code = c_matrix_propagate(
                solution_storage_ptr,
                frequency_to_use,
                bulk_density_to_use,
                first_slice_index_by_layer_vec.data(),
                num_slices_by_layer_vec.data(),
                num_layers,
                num_bc_models,
                bc_models_ptr,
                G_to_use,
                degree_l,
                starting_radius_to_use,
                start_radius_tolerance,
                core_model,
                verbose
            );
        }
        else
        {
            sub_process_error_code = c_shooting_solver(
                solution_storage_ptr,
                frequency_to_use,
                bulk_density_to_use,
                layer_types_ptr,
                is_static_bylayer_ptr,
                is_incompressible_bylayer_ptr,
                first_slice_index_by_layer_vec,
                num_slices_by_layer_vec,
                num_bc_models,
                bc_models_ptr,
                G_to_use,
                degree_l,
                use_kamata,
                starting_radius_to_use,
                start_radius_tolerance,
                integration_method_int,
                integration_rtol,
                integration_atol,
                scale_rtols_bylayer_type,
                max_num_steps,
                expected_size,
                max_ram_MB,
                max_step,
                verbose
            );
        }
    }

    // Finalize
    if (nondimensionalize)
    {
        solution_storage_ptr->dimensionalize_data(&non_dim_scales, true);

        for (size_t slice_i = 0; slice_i < total_slices; ++slice_i)
        {
            radius_array_in_ptr[slice_i]          *= non_dim_scales.length_conversion;
            density_array_in_ptr[slice_i]         *= non_dim_scales.density_conversion;
            complex_bulk_modulus_in_ptr[slice_i]  *= non_dim_scales.pascal_conversion;
            complex_shear_modulus_in_ptr[slice_i] *= non_dim_scales.pascal_conversion;
        }
    }

    if (solution_storage_ptr->success)
    {
        solution_storage_ptr->find_love();
    }

    return solution_storage_ptr->error_code;
}


// =================================================================================================
// MOVED FROM `DEF` FUNC
// =================================================================================================
// This helper function captures the input validation, string mapping, and sanity checks originally 
// performed in the `radial_solver` Python wrapper.

std::string to_lower(const std::string& input)
{
    std::string result = input;
    for (char& c : result)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return result;
}

void c_validate_and_prep_radial_inputs(
    size_t total_slices,
    const double* radius_array,
    const double* density_array,
    double frequency,
    size_t num_layers,
    const std::vector<std::string>& layer_types,
    const bool* is_static_bylayer,
    const bool* is_incompressible_bylayer,
    const double* upper_radius_bylayer_array,
    bool use_prop_matrix,
    double starting_radius,
    const std::vector<std::string>& solve_for,
    const std::string& integration_method,
    const std::vector<std::string>& eos_method_bylayer,
    const std::string& eos_integration_method,
    bool warnings,
    int* layer_types_out_ptr,
    int* bc_models_out_ptr,
    size_t& num_bc_models_out,
    ODEMethod& integration_method_out,
    std::vector<int>& eos_integration_method_int_bylayer_out,
    ODEMethod& eos_integration_method_out
)
{
    // Layer/Array dimension checks
    if (layer_types.size() != num_layers)
        throw std::invalid_argument("Number of `layer_types` must match `num_layers`.");

    // Ascending layer order
    double last_layer_r = 0.0;
    for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        if (upper_radius_bylayer_array[layer_i] <= last_layer_r)
        {
            throw std::invalid_argument("`upper_radius_bylayer_array` must be in ascending order.");
        }
        last_layer_r = upper_radius_bylayer_array[layer_i];
    }

    // Frequency checks
    if (std::abs(frequency) < tidalpy_config_ptr->d_MIN_FREQUENCY)
        throw std::invalid_argument("Forcing frequency is too small (are you sure you are in rad s-1?).");
    else if (std::abs(frequency) > tidalpy_config_ptr->d_MAX_FREQUENCY)
        throw std::invalid_argument("Forcing frequency is too large (are you sure you are in rad s-1?).");

    // Propagation matrix limitation checks
    if (use_prop_matrix)
    {
        if (num_layers > 1)
            throw std::logic_error("Currently, TidalPy's propagation matrix technique only works for 1-layer worlds.");
        if (to_lower(layer_types[0]) != "solid")
            throw std::invalid_argument("The Propagation matrix technique only works for solid layers.");
        if (not is_static_bylayer[0])
            throw std::invalid_argument("The Propagation matrix technique does not allow for dynamic layers.");
        if (not is_incompressible_bylayer[0])
            throw std::invalid_argument("The Propagation matrix technique does not allow for compressible layers.");
    }

    if ((starting_radius != 0.0) && (starting_radius > 0.90 * radius_array[total_slices - 1]))
    {
        throw std::invalid_argument("Starting radius is above 90% of the planet radius. Try a lower radius.");
    }

    if (radius_array[0] != 0.0)
    {
        throw std::invalid_argument("Radius array must start at zero.");
    }

    // Slice count and interface verification
    double last_layer_radius = 0.0;
    for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        bool top_layer      = (layer_i == num_layers - 1);
        double layer_radius = upper_radius_bylayer_array[layer_i];

        size_t slice_check       = 0;
        size_t layer_check       = 0;
        double last_radius_check = 0.0;

        for (size_t slice_i = 0; slice_i < total_slices; ++slice_i)
        {
            double radius_check = radius_array[slice_i];
            if (radius_check < 0.0)
                throw std::invalid_argument("A negative radius value was found in `radius_array`.");
            if (radius_check < last_radius_check)
                throw std::invalid_argument("Radius array must be in ascending order.");
            if (c_isclose(radius_check, layer_radius, 1.0e-9, 0.0))
                layer_check += 1;
            if (last_layer_radius <= radius_check && radius_check <= layer_radius)
                slice_check += 1;
            
            last_radius_check = radius_check;
        }

        last_layer_radius = layer_radius;

        if (slice_check < 5)
            throw std::invalid_argument("A minimum of 5 sub-slices (including top and bottom) are required for each layer.");

        if (top_layer)
        {
            if (layer_check != 1)
                throw std::invalid_argument("Radius of layer " + std::to_string(layer_i) + " found " + std::to_string(layer_check) + " times. Expected 1 time (non-interface layer).");
        }
        else
        {
            if (layer_check != 2)
                throw std::invalid_argument("Radius of layer " + std::to_string(layer_i) + " found " + std::to_string(layer_check) + " times. Expected 2 times (interface layer).");
        }
    }

    // Layer type string to int processing
    bool dynamic_liquid = false;
    for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        const std::string l_type = to_lower(layer_types[layer_i]);

        if (l_type == "solid")
        {
            layer_types_out_ptr[layer_i] = 0;
        }
        else if (l_type == "liquid")
        {
            layer_types_out_ptr[layer_i] = 1;
            if (not is_static_bylayer[layer_i])
            {
                dynamic_liquid = true;
            }
        }
        else
        {
            layer_types_out_ptr[layer_i] = -1;
            throw std::invalid_argument("Layer type " + layer_types[layer_i] + " is not supported.");
        }
    }

    if (dynamic_liquid && std::abs(frequency) < 2.5e-5 && warnings)
    {
        printf("WARNING: Dynamic liquid layer detected in RadialSolver for a small frequency. Results may be unstable. Extra care is advised!\n");
    }

    // Integration methods parsing
    std::string int_method_lower = to_lower(integration_method);
    if (int_method_lower == "rk45")        integration_method_out = ODEMethod::RK45;
    else if (int_method_lower == "rk23")   integration_method_out = ODEMethod::RK23;
    else if (int_method_lower == "dop853") integration_method_out = ODEMethod::DOP853;
    else throw std::invalid_argument("Unsupported integration method provided: " + int_method_lower);

    std::string eos_int_method_lower = to_lower(eos_integration_method);
    if (eos_int_method_lower == "rk45")        eos_integration_method_out = ODEMethod::RK45;
    else if (eos_int_method_lower == "rk23")   eos_integration_method_out = ODEMethod::RK23;
    else if (eos_int_method_lower == "dop853") eos_integration_method_out = ODEMethod::DOP853;
    else throw std::invalid_argument("Unsupported EOS integration method provided: " + eos_int_method_lower);

    // EOS integration methods by layer
    eos_integration_method_int_bylayer_out.resize(num_layers);
    if (eos_method_bylayer.empty())
    {
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
        {
            eos_integration_method_int_bylayer_out[layer_i] = C_EOS_INTERPOLATE_METHOD_INT;
        }
    }
    else
    {
        for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
        {
            if (to_lower(eos_method_bylayer[layer_i]) == "interpolate")
            {
                eos_integration_method_int_bylayer_out[layer_i] = C_EOS_INTERPOLATE_METHOD_INT;
            }
            else
            {
                throw std::logic_error("Unknown EOS method provided: " + eos_method_bylayer[layer_i]);
            }
        }
    }

    // Boundary condition models map
    num_bc_models_out = solve_for.size();
    if (num_bc_models_out == 0)
    {
        num_bc_models_out = 1;
        bc_models_out_ptr[0] = 1;
    }
    else
    {
        for (size_t i = 0; i < num_bc_models_out; ++i)
        {
            std::string solve_tmp = to_lower(solve_for[i]);
            if (solve_tmp == "free")         bc_models_out_ptr[i] = 0;
            else if (solve_tmp == "tidal")   bc_models_out_ptr[i] = 1;
            else if (solve_tmp == "loading") bc_models_out_ptr[i] = 2;
            else throw std::invalid_argument("Unsupported value provided for `solve_for`: " + solve_tmp);
        }
    }
}