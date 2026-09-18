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
/// Integrates gravity, pressure, mass, and moment of inertia from center to surface in the caller's units (SI or
/// non-dimensional). The central pressure is found by a secant iteration on the surface-pressure mismatch: the
/// first update assumes a unit slope (exact for an incompressible planet), later ones use the measured slope.
///
/// Parameters
/// ----------
/// eos_solution_ptr : c_EOSSolution*
///     Output; must be constructed with the radius array and layer info.
/// eos_function_bylayer_ptr_vec : vector of PreEvalFunc
///     EOS evaluation function per layer (used during integration and kept for later evaluation).
/// eos_input_bylayer_vec : vector of c_EOS_ODEInput
///     Input parameters per layer.
/// planet_bulk_density : double
///     Seeds the central-pressure guess.
/// surface_pressure : double
///     Target surface pressure.
/// G_to_use : double
///     Gravitational constant in solve units; negative selects the shared config's SI value.
/// integration_method, rtol, atol
///     CyRK method and tolerances.
/// pressure_tol : double
///     Surface-pressure convergence tolerance relative to the central-pressure scale
///     (2/3) pi G rho_bulk^2 R^2 + surface_pressure. Must exceed rtol, the integrator's own noise, to converge.
/// max_iters : size_t
///     Iteration cap; reported through max_iters_hit and the message, the solution is still returned.
/// verbose : bool
///     Print status messages.
///
/// Assumptions
/// -----------
/// - Spherical symmetry and hydrostatic equilibrium.
/// - The surface pressure rises monotonically with the central pressure (any EOS with positive density and
///   compressibility), which the secant iteration relies on.
inline void c_solve_eos(
        c_EOSSolution* eos_solution_ptr,
        std::vector<PreEvalFunc>& eos_function_bylayer_ptr_vec,
        std::vector<c_EOS_ODEInput>& eos_input_bylayer_vec,
        double planet_bulk_density,
        double surface_pressure,
        double G_to_use,
        ODEMethod integration_method,
        double rtol,
        double atol,
        double pressure_tol,
        size_t max_iters,
        bool verbose
        ) noexcept
{
    eos_solution_ptr->message = std::string("Equation of state solver finished without issue.");

    if (G_to_use < 0.0) { G_to_use = c_get_G(); }

    // One tolerance for every y; CyRK still takes vectors.
    std::vector<double> rtols_vec = {rtol};
    std::vector<double> atols_vec = {atol};

    const size_t len_radius_array  = eos_solution_ptr->radius_array_vec.size();
    const double planet_radius     = eos_solution_ptr->radius_array_vec.back();
    const double r0_gravity        = 0.0;
    const double r0_pressure_guess = (
        (2.0 / 3.0) * TidalPyConstants::d_PI * G_to_use * planet_radius * 
        planet_radius * planet_bulk_density * planet_bulk_density)
        + surface_pressure;

    const double r0_mass = 0.0;
    const double r0_moi  = 0.0;

    double radius_start   = 0.0;
    double radius_stop    = 0.0;
    DiffeqFuncType diffeq = c_eos_diffeq;

    // The initial central pressure is that of a uniform sphere at the bulk density.
    double y0[4] = {r0_gravity, r0_pressure_guess, r0_mass, r0_moi};

    // Secant iteration state on f(P_c) = P_surface(P_c) - surface_pressure. The convergence test is relative to
    // the central-pressure scale, which is the size of the integrator's own noise on the surface pressure.
    const double pressure_scale   = (r0_pressure_guess > 0.0) ? r0_pressure_guess : 1.0;
    double previous_central       = TidalPyConstants::d_NAN;
    double previous_diff          = TidalPyConstants::d_NAN;

    // Only the four structure variables are integrated; density, moduli, and viscosities are evaluated afterwards
    // from the retained dense output, so the dense interpolant stays a plain polynomial evaluation.
    const size_t num_y     = C_EOS_Y_VALUES;
    const size_t num_extra = 0;

    size_t top_of_last_layer_index = 0;
    std::vector<double> y0_bylayer_vec(4);

    std::vector<char> args_vec(sizeof(c_EOS_ODEInput));
    c_EOS_ODEInput* args_vec_eos_ptr = reinterpret_cast<c_EOS_ODEInput*>(args_vec.data());
    c_EOS_ODEInput* eos_input_layer_ptr = nullptr;

    double max_step              = 0.0;
    size_t max_num_steps         = 10000;
    size_t max_ram_MB            = 512;
    bool use_dense_output        = false;
    double first_step            = 0.0;
    size_t expected_size         = 64;
    std::vector<double> t_eval_vec(0);
    PreEvalFunc layer_eos_func   = nullptr;

    std::vector<Event> events_vec;  // unused by the EOS solve

    double calculated_surf_pressure = TidalPyConstants::d_INF;
    double pressure_diff            = TidalPyConstants::d_INF;
    double pressure_diff_abs        = TidalPyConstants::d_INF;
    int iterations                  = 0;
    bool failed                     = false;
    bool max_iters_hit              = false;
    bool final_run                  = false;
    std::string integrator_failure_message;

    size_t last_solution_size = 0;
    std::unique_ptr<CySolverResult> integration_result_uptr = std::make_unique<CySolverResult>(integration_method);
    CySolverResult* integration_result_ptr = nullptr;

    const size_t num_layers = eos_solution_ptr->num_layers;

    // Surface-pressure convergence loop.
    while (true)
    {
        calculated_surf_pressure = TidalPyConstants::d_INF;

        if (!final_run)
        {
            iterations++;
        }

        // Integrate layer by layer from the center outward.
        for (size_t layer_i = 0; layer_i < num_layers; layer_i++)
        {
            radius_stop = eos_solution_ptr->upper_radius_bylayer_vec[layer_i];
            if (layer_i == 0)
            {
                radius_start = 0.0;
                for (size_t y_i = 0; y_i < num_y; y_i++)
                {
                    y0_bylayer_vec[y_i] = y0[y_i];
                }
            }

            // Maximum step of one third of the layer thickness.
            max_step = 0.33 * (radius_stop - radius_start);

            eos_input_layer_ptr = &eos_input_bylayer_vec[layer_i];

            // The integration needs only the density and captures no extra outputs, so the diffeq must not write
            // past the structure variables. Only the final run keeps its dense output.
            eos_input_layer_ptr->update_bulk  = false;
            eos_input_layer_ptr->update_shear = false;
            eos_input_layer_ptr->final_solve  = false;
            use_dense_output = final_run;

            std::memcpy(args_vec_eos_ptr, eos_input_layer_ptr, sizeof(c_EOS_ODEInput));

            layer_eos_func = eos_function_bylayer_ptr_vec[layer_i];

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
                true,              // Force retain solver [bool]
                nullptr            // Analytic jacobian (null = numerical; used by implicit methods) [JacobianFuncType]
            );
            last_solution_size = integration_result_ptr->size;
            eos_solution_ptr->save_steps_taken(integration_result_ptr->steps_taken);

            if (!integration_result_ptr->success)
            {
                failed = true;
                // The result pointer is cleared every layer, so capture the message now.
                integrator_failure_message = integration_result_ptr->message;
            }

            if (final_run && !failed)
            {
                // The whole result is kept for later dense calls; repoint after the move.
                eos_solution_ptr->save_cyresult(std::move(integration_result_uptr));
                integration_result_ptr = eos_solution_ptr->cysolver_results_uptr_bylayer_vec.back().get();
            }
            else if ((layer_i == num_layers - 1) && !failed)
            {
                // Pressure of the last step: (total slices) - (num_y - pressure index) - 1.
                size_t surface_pressure_index = (last_solution_size * num_y) - (num_y - 2) - 1;
                calculated_surf_pressure      = integration_result_ptr->solution[surface_pressure_index];
            }

            if ((num_layers > 1) && !failed)
            {
                radius_start = eos_solution_ptr->upper_radius_bylayer_vec[layer_i];
                top_of_last_layer_index = (num_extra + num_y) * (last_solution_size - 1);

                // The next layer starts from the top of this one.
                if (integration_result_ptr)
                {
                    std::memcpy(
                        y0_bylayer_vec.data(),
                        &integration_result_ptr->solution[top_of_last_layer_index],
                        sizeof(double) * num_y);
                }
                else
                {
                    failed = true;
                }
            }

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
            break;
        }
        else
        {
            pressure_diff     = calculated_surf_pressure - surface_pressure;
            pressure_diff_abs = std::fabs(pressure_diff);

            if (pressure_diff_abs <= pressure_tol * pressure_scale)
            {
                // Converged: the next pass keeps the dense output.
                final_run = true;
            }
            else
            {
                // Secant update; the unit slope is the fallback when the measured slope is unusable.
                double step = -pressure_diff;
                if (std::isfinite(previous_central))
                {
                    const double slope = (pressure_diff - previous_diff) / (y0[1] - previous_central);
                    if (std::isfinite(slope) && (slope > 0.0))
                    {
                        step = -pressure_diff / slope;
                    }
                }
                previous_central = y0[1];
                previous_diff    = pressure_diff;

                // Keep the central pressure positive: halve an overshooting step.
                double next_central = y0[1] + step;
                while ((next_central <= 0.0) && (std::fabs(step) > TidalPyConstants::d_EPS * pressure_scale))
                {
                    step        *= 0.5;
                    next_central = y0[1] + step;
                }
                y0[1] = next_central;
            }
        }

        if (iterations >= static_cast<int>(max_iters))
        {
            max_iters_hit = true;
            eos_solution_ptr->max_iters_hit = true;
            // Still produce output.
            final_run = true;
        }
    }

    eos_solution_ptr->iterations = iterations;

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
        if (!integrator_failure_message.empty())
        {
            eos_solution_ptr->message += std::string(". Message: ") + integrator_failure_message;
        }
        if (verbose)
        {
            std::printf("%s", eos_solution_ptr->message.c_str());
        }
    }
    else
    {
        eos_solution_ptr->success = true;
        eos_solution_ptr->pressure_error = pressure_diff_abs;

        // Keep the layer EOS functions for on-demand evaluation, then sample the planet onto the radius array.
        eos_solution_ptr->save_eos_functions(eos_function_bylayer_ptr_vec, eos_input_bylayer_vec);
        eos_solution_ptr->interpolate_full_planet();
    }

    if (integration_result_uptr)
    {
        integration_result_uptr.reset();
    }

    integration_result_ptr = nullptr;
}
