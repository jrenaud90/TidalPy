#pragma once

#include <cstdio>
#include <cstring>
#include <vector>
#include <memory>
#include <string>
#include <cmath>

#include "c_common.hpp"    // CyRK: DiffeqFuncType, PreEvalFunc, ODEMethod
#include "cysolution.hpp"  // CyRK: CySolverResult
#include "cysolve.hpp"     // CyRK: baseline_cysolve_ivp_noreturn
#include "c_events.hpp"    // CyRK: Event

#include "constants_.hpp"  // TidalPy: TidalPyConstants, d_PI, d_INF, d_EPS_100

#include "ode_.hpp"            // c_eos_diffeq, c_EOS_ODEInput, C_EOS_Y_VALUES, C_EOS_EXTRA_VALUES
#include "eos_solution_.hpp"   // c_EOSSolution


/// Solve the equation of state for a layered planet.
///
/// Integrates gravity, pressure, mass, and moment of inertia radially from center to surface.
/// Uses a convergence loop on surface pressure to determine the correct central pressure.
///
/// Parameters
/// ----------
/// eos_solution_ptr : c_EOSSolution*
///     Output solution object (must be constructed with radius array and layer info).
/// eos_function_bylayer_ptr_vec : vector of PreEvalFunc
///     EOS evaluation function for each layer (called during ODE integration).
/// eos_input_bylayer_vec : vector of c_EOS_ODEInput
///     Input parameters for each layer's EOS function.
/// planet_bulk_density : double
///     Bulk density of the planet [kg m-3], used for initial pressure guess.
/// surface_pressure : double
///     Expected surface pressure [Pa] (default: 0.0).
/// G_to_use : double
///     Gravitational constant [m3 kg-1 s-2].
/// integration_method : ODEMethod
///     CyRK integration method (default: DOP853).
/// rtol : double
///     Relative tolerance for integration.
/// atol : double
///     Absolute tolerance for integration.
/// pressure_tol : double
///     Convergence tolerance for surface pressure iteration.
/// max_iters : size_t
///     Maximum number of convergence iterations.
/// verbose : bool
///     Print status messages if true.
inline void c_solve_eos(
        c_EOSSolution* eos_solution_ptr,
        std::vector<PreEvalFunc>& eos_function_bylayer_ptr_vec,
        std::vector<c_EOS_ODEInput>& eos_input_bylayer_vec,
        double planet_bulk_density,
        double surface_pressure = 0.0,
        double G_to_use = 6.674015e-11,
        ODEMethod integration_method = ODEMethod::DOP853,
        double rtol = 1.0e-6,
        double atol = 1.0e-10,
        double pressure_tol = 1.0e-3,
        size_t max_iters = 100,
        bool verbose = true
        ) noexcept
{
    // Set the message assuming success, it will be updated if we run into failure
    eos_solution_ptr->message = std::string("Equation of state solver finished without issue.");

    // We will just use one rtol and one atol for all y's but still need to provide it as a vector.
    std::vector<double> rtols_vec = {rtol};
    std::vector<double> atols_vec = {atol};

    // Determine planetary properties
    const size_t len_radius_array  = eos_solution_ptr->radius_array_vec.size();
    const double planet_radius     = eos_solution_ptr->radius_array_vec.back();
    const double r0_gravity        = 0.0;
    const double r0_pressure_guess = (
        (2.0 / 3.0) * TidalPyConstants::d_PI * G_to_use * planet_radius * 
        planet_radius * planet_bulk_density * planet_bulk_density)
        + surface_pressure;

    const double r0_mass = 0.0;
    const double r0_moi  = 0.0;

    // Setup bound variables
    double radius_start   = 0.0;
    double radius_stop    = 0.0;
    DiffeqFuncType diffeq = c_eos_diffeq;

    // We need the central pressure of the planet. Use the global bulk density to calculate this.
    double y0[4] = {r0_gravity, r0_pressure_guess, r0_mass, r0_moi};

    // y information
    const size_t num_y = C_EOS_Y_VALUES;
    size_t num_extra   = 0;

    // Layer information
    size_t top_of_last_layer_index = 0;
    std::vector<double> y0_bylayer_vec(4);

    // EOS functions and inputs
    std::vector<char> args_vec(sizeof(c_EOS_ODEInput));
    c_EOS_ODEInput* args_vec_eos_ptr = reinterpret_cast<c_EOS_ODEInput*>(args_vec.data());
    c_EOS_ODEInput* eos_input_layer_ptr = nullptr;

    // Other integration information
    double max_step              = 0.0;
    size_t max_num_steps         = 10000;
    size_t max_ram_MB            = 512;
    bool use_dense_output        = false;
    double first_step            = 0.0;
    size_t expected_size         = 64;
    std::vector<double> t_eval_vec(0);
    PreEvalFunc layer_eos_func   = nullptr;

    // Events (empty — not used by EOS solver)
    std::vector<Event> events_vec;

    // Pressure convergence variables
    double calculated_surf_pressure = TidalPyConstants::d_INF;
    double pressure_diff            = TidalPyConstants::d_INF;
    double pressure_diff_abs        = TidalPyConstants::d_INF;
    int iterations                  = 0;
    bool failed                     = false;
    bool max_iters_hit              = false;
    bool final_run                  = false;

    // Integration solution variables
    size_t last_solution_size = 0;
    std::unique_ptr<CySolverResult> integration_result_uptr = std::make_unique<CySolverResult>(integration_method);
    CySolverResult* integration_result_ptr = nullptr;

    // Loop variables
    const size_t num_layers = eos_solution_ptr->num_layers;

    // Solve the equation of state in a convergence loop based on the surface pressure.
    while (true)
    {
        // Reset calculated pressure
        calculated_surf_pressure = TidalPyConstants::d_INF;

        if (!final_run)
        {
            iterations++;
        }

        // Step through each macro layer of the planet and solve the equation of state starting from bottom to top
        for (size_t layer_i = 0; layer_i < num_layers; layer_i++)
        {
            // Setup bounds and initial conditions for next layer's integration
            radius_stop = eos_solution_ptr->upper_radius_bylayer_vec[layer_i];
            if (layer_i == 0)
            {
                radius_start = 0.0;
                // Set y0 for bottom-most layer equal to the global y0
                for (size_t y_i = 0; y_i < num_y; y_i++)
                {
                    y0_bylayer_vec[y_i] = y0[y_i];
                }
            }

            // Set the maximum step size equal to 1/3 the layer's thickness
            max_step = 0.33 * (radius_stop - radius_start);

            // Get eos function and inputs for this layer
            eos_input_layer_ptr = &eos_input_bylayer_vec[layer_i];

            if (final_run)
            {
                // We now want to make sure that all final calculations are performed.
                eos_input_layer_ptr->update_bulk  = true;
                eos_input_layer_ptr->update_shear = true;
                eos_input_layer_ptr->final_solve  = true;
                // Capture extra outputs and store interpolators
                num_extra        = C_EOS_EXTRA_VALUES;
                use_dense_output = true;
            }
            else
            {
                // During the iterations we do not need to update the complex bulk or shear
                eos_input_layer_ptr->update_bulk  = false;
                eos_input_layer_ptr->update_shear = false;
                // We also are not at the final call step.
                eos_input_layer_ptr->final_solve = false;
                num_extra        = 0;
                use_dense_output = false;
            }

            // Store additional arguments to our char vector
            std::memcpy(args_vec_eos_ptr, eos_input_layer_ptr, sizeof(c_EOS_ODEInput));

            // Get layer-specific eos function, called during the EOS diffeq.
            layer_eos_func = eos_function_bylayer_ptr_vec[layer_i];

            ///// Radial Integrate the EOS Through the Planet /////
            if (!integration_result_uptr)
            {
                integration_result_uptr = std::make_unique<CySolverResult>(integration_method);
            }
            integration_result_ptr = integration_result_uptr.get();

            baseline_cysolve_ivp_noreturn(
                integration_result_ptr,
                diffeq,            // Differential equation [DiffeqFuncType]
                radius_start,      // Start radius for this layer
                radius_stop,       // Stop radius for this layer
                y0_bylayer_vec,    // y0 array vector<double>
                expected_size,     // Expected final integration size [size_t]
                num_extra,         // Number of extra outputs tracked [size_t]
                args_vec,          // Extra input args to diffeq vector[char]
                max_num_steps,     // Max number of steps [size_t]
                max_ram_MB,        // Max amount of RAM allowed [size_t]
                use_dense_output,  // Use dense output [bool]
                t_eval_vec,        // Interpolate at radius array vector[double]
                layer_eos_func,    // Pre-eval function used in diffeq [PreEvalFunc]
                events_vec,        // Events vector [vector<Event>]
                rtols_vec,         // Relative Tolerance vector[double]
                atols_vec,         // Absolute Tolerance vector[double]
                max_step,          // Maximum step size [double]
                first_step,        // Initial step size [double]
                true               // Force retain solver [bool]
            );
            /////////////////////////////////////////////////////
            last_solution_size = integration_result_ptr->size;
            eos_solution_ptr->save_steps_taken(integration_result_ptr->steps_taken);

            if (!integration_result_ptr->success)
            {
                failed = true;
            }

            if (final_run && !failed)
            {
                // Save the current cysolver result
                // (we need to save the whole object so we can make interpolator calls to it later)
                eos_solution_ptr->save_cyresult(std::move(integration_result_uptr));

                // Change where the integrator result pointer is pointing to since we moved the unique pointer
                integration_result_ptr = eos_solution_ptr->cysolver_results_uptr_bylayer_vec.back().get();
            }
            else if ((layer_i == num_layers - 1) && !failed)
            {
                // Find planet surface pressure for this iteration
                // (Total number of slices) - (num_y - location of pressure) - 1
                size_t surface_pressure_index = (last_solution_size * num_y) - (num_y - 2) - 1;
                calculated_surf_pressure      = integration_result_ptr->solution[surface_pressure_index];
            }

            // Prepare for next layer
            if ((num_layers > 1) && !failed)
            {
                // Bottom radius value equals top of lower layer's radius
                radius_start = eos_solution_ptr->upper_radius_bylayer_vec[layer_i];
                top_of_last_layer_index = (num_extra + num_y) * (last_solution_size - 1);

                // y0 for this layer equals the top most result of the lower layer
                if (integration_result_ptr)
                {
                    std::memcpy(
                        y0_bylayer_vec.data(),
                        &integration_result_ptr->solution[top_of_last_layer_index],
                        sizeof(double) * num_y);
                }
                else
                {
                    // Not sure why that would be null but in any case we are in a fail state.
                    failed = true;
                }
            }

            // Clear pointers
            integration_result_ptr = nullptr;

            if (failed)
            {
                break;
            }
        }

        if (failed)
        {
            break;
        }

        if (final_run)
        {
            // We are done!
            break;
        }
        else
        {
            // Update the central pressure using the error at the surface as the correction factor
            pressure_diff     = surface_pressure - calculated_surf_pressure;
            pressure_diff_abs = pressure_diff;
            if (pressure_diff < 0.0)
            {
                pressure_diff_abs = -pressure_diff;
            }

            // Calculate percent difference to use in convergence check.
            if (surface_pressure > TidalPyConstants::d_EPS_100)
            {
                pressure_diff_abs /= surface_pressure;
            }

            // Check if we are done next iteration
            if (pressure_diff_abs <= pressure_tol)
            {
                final_run = true;
            }
            else
            {
                y0[1] += pressure_diff;
            }
        }

        if (iterations >= static_cast<int>(max_iters))
        {
            max_iters_hit = true;
            eos_solution_ptr->max_iters_hit = true;
            // To ensure that there is some output we will go ahead and do a final run.
            final_run = true;
        }
    }

    // Done with convergence loop.
    eos_solution_ptr->iterations = iterations;

    // Display any warnings
    if (max_iters_hit)
    {
        eos_solution_ptr->message = std::string("Warning in `c_solve_eos`: Maximum number of iterations hit without convergence.");
        if (verbose)
        {
            std::printf("%s", eos_solution_ptr->message.c_str());
        }
    }

    if (failed)
    {
        eos_solution_ptr->success = false;
        eos_solution_ptr->message = std::string("Warning in `c_solve_eos`: Integrator failed at iteration ") + std::to_string(iterations);
        if (integration_result_ptr)
        {
            eos_solution_ptr->message += std::string(". Message: ") + integration_result_ptr->message;
        }
        if (verbose)
        {
            std::printf("%s", eos_solution_ptr->message.c_str());
        }
    }
    else
    {
        // Set feedback attributes
        eos_solution_ptr->success = true;

        // Set other final parameters
        eos_solution_ptr->pressure_error = pressure_diff_abs;

        // Tell the eos solution to perform a full planet interpolation and store the results. Including surface results.
        eos_solution_ptr->interpolate_full_planet();
    }

    if (integration_result_uptr)
    {
        integration_result_uptr.reset();
    }

    integration_result_ptr = nullptr;
}
