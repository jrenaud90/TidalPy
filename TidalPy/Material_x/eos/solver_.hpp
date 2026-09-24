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
/// Integrates gravity, pressure, mass, and moment of inertia from center to surface in the caller's units.
/// The central pressure comes from a safeguarded secant iteration on the surface-pressure mismatch. The first
/// update assumes a unit slope, exact for an incompressible planet, and later ones use the measured slope. Where
/// that slope is not positive (the surface pressure of a compressible planet can fall as the central pressure
/// first rises), the step grows geometrically in the unit-slope direction until the mismatch changes sign, and
/// until then no secant step may grow more than fourfold over the one before. Once
/// it has, the root is bracketed: a secant step that would leave the bracket is replaced by an Illinois
/// false-position step inside it, and any step not under half the one two passes back by a bisection (Brent's
/// safeguard), so a mismatch with a near-jump cannot stall the iteration.
///
/// Only the pass that converges needs its dense output, and capturing it roughly doubles the cost of a
/// pass, so it is switched on for the first pass of a warm start (which usually converges there) and for
/// any pass the secant's own error model expects to converge. A pass that converges without it is repeated.
///
/// Parameters
/// ----------
/// eos_solution_ptr : c_EOSSolution*
///     Output; must be constructed with the radius array and layer info.
/// eos_function_bylayer_ptr_vec, eos_input_bylayer_vec
///     EOS evaluation function and arguments per layer; used during the integration and kept for later.
/// planet_bulk_density : double
///     Seeds the central-pressure guess.
/// surface_pressure : double
///     Target surface pressure.
/// G_to_use : double
///     Gravitational constant in solve units; negative selects the shared config's SI value.
/// integration_method, rtol, atol
///     CyRK method and tolerances.
/// pressure_tol : double
///     Surface-pressure tolerance relative to the central-pressure scale (2/3) pi G rho_bulk^2 R^2 +
///     surface_pressure. Must exceed rtol, the integrator's own noise, to converge.
/// max_iters : size_t
///     Iteration cap; reported through max_iters_hit, and the solution is still returned.
/// verbose : bool
///     Print status messages.
/// segment_vec_ptr : const vector of c_EOSSegment, optional
///     The radial segments to integrate, ascending, covering the planet. Null integrates one segment per
///     layer, the layout of a solve that does not carry temperature.
/// integrate_temperature : bool
///     Carry temperature and heat flow as two extra states, each segment's gradient form taken from the
///     layout. Requires a segment layout.
/// central_pressure_guess : double
///     First central pressure to try: the last converged value when the caller has one. NaN or a
///     non-positive value starts from a uniform sphere at the bulk density.
///
/// Assumptions
/// -----------
/// - Spherical symmetry and hydrostatic equilibrium.
/// - The surface pressure rises with the central pressure near the root, and the unit-slope direction leads
///   toward it: too high a surface pressure lowers the central pressure. A structure still off its target surface
///   pressure at `max_iters` is reported as a failure.
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
        bool verbose,
        const std::vector<c_EOSSegment>* segment_vec_ptr = nullptr,
        bool integrate_temperature = false,
        double central_pressure_guess = TidalPyConstants::d_NAN
        )
{
    eos_solution_ptr->message = std::string("Equation of state solver finished without issue.");

    if (G_to_use < 0.0) { G_to_use = c_get_G(); }

    // One tolerance for every y; CyRK still takes vectors.
    std::vector<double> rtols_vec = {rtol};
    std::vector<double> atols_vec = {atol};

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

    // One segment per layer unless the caller supplied a layout.
    std::vector<c_EOSSegment> segment_vec;
    if (segment_vec_ptr == nullptr)
    {
        segment_vec.resize(eos_solution_ptr->num_layers);
        for (size_t layer_i = 0; layer_i < eos_solution_ptr->num_layers; ++layer_i)
        {
            segment_vec[layer_i].upper_radius = eos_solution_ptr->upper_radius_bylayer_vec[layer_i];
            segment_vec[layer_i].layer_index  = layer_i;
        }
        integrate_temperature = false;
    }
    else
    {
        segment_vec = *segment_vec_ptr;
    }
    const size_t num_segments = segment_vec.size();
    eos_solution_ptr->set_segments(segment_vec);

    DiffeqFuncType diffeq = integrate_temperature ? c_eos_diffeq_thermal : c_eos_diffeq;

    // The caller's central pressure, or that of a uniform sphere at the bulk density. A thermal solve
    // starts at the innermost segment's temperature with whatever heat flow it carries, zero at a regular
    // center.
    const bool warm_start    = std::isfinite(central_pressure_guess) && (central_pressure_guess > 0.0);
    const double r0_pressure = warm_start ? central_pressure_guess : r0_pressure_guess;
    double y0[C_EOS_THERMAL_Y_VALUES] = {r0_gravity, r0_pressure, r0_mass, r0_moi, 0.0, 0.0};
    if (integrate_temperature && (num_segments > 0))
    {
        y0[4] = segment_vec[0].start_temperature;
        y0[5] = segment_vec[0].start_heat_flow;
    }

    // Secant state on f(P_c) = P_surface(P_c) - surface_pressure. The convergence test is relative to the
    // central-pressure scale, the size of the integrator's own noise on the surface pressure.
    const double pressure_scale   = (r0_pressure_guess > 0.0) ? r0_pressure_guess : 1.0;
    const double pressure_tol_abs = pressure_tol * pressure_scale;
    double previous_central       = TidalPyConstants::d_NAN;
    double previous_diff          = TidalPyConstants::d_NAN;
    double oldest_diff            = TidalPyConstants::d_NAN;
    // The latest central pressure on each side of the root, with its mismatch; the root lies between them once
    // both are known. `last_replaced` is the side the last pass replaced (-1 below, +1 above, 0 none), and
    // `last_bracket_step` whether the step to it replaced a secant step inside the bracket (false position or
    // bisection), for the Illinois halving of an end kept twice.
    double below_central          = TidalPyConstants::d_NAN;   // mismatch < 0 there
    double below_diff             = TidalPyConstants::d_NAN;
    double above_central          = TidalPyConstants::d_NAN;   // mismatch >= 0 there
    double above_diff             = TidalPyConstants::d_NAN;
    int    last_replaced          = 0;
    bool   last_bracket_step      = false;
    // Sizes of the last two steps: a run of unusable slopes doubles the last, and the one before it bounds a
    // bracketed step (Brent's safeguard).
    double last_step_abs          = 0.0;
    double step_two_ago_abs       = TidalPyConstants::d_INF;

    // Only the four structure variables are integrated; density, moduli, and viscosities follow from the
    // retained dense output, so the interpolant stays a plain polynomial evaluation.
    const size_t num_y     = integrate_temperature ? C_EOS_THERMAL_Y_VALUES : C_EOS_Y_VALUES;
    const size_t num_extra = 0;
    eos_solution_ptr->num_y_solved = num_y;

    size_t top_of_last_segment_index = 0;
    // CyRK takes the number of state variables from this vector, so it must hold exactly num_y of them.
    std::vector<double> y0_bysegment_vec(num_y);

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
    // Whether this pass keeps its dense output, and whether it is kept whatever it finds: the repeat of a
    // converged pass, or the pass after the iteration cap. Capturing roughly doubles the cost of a pass,
    // so the first pass does it only from a warm start, which usually converges there.
    bool capture_dense              = warm_start;
    bool final_pass                 = false;
    std::vector<std::unique_ptr<CySolverResult>> pass_results_vec;
    std::string integrator_failure_message;

    size_t last_solution_size = 0;
    std::unique_ptr<CySolverResult> integration_result_uptr = std::make_unique<CySolverResult>(integration_method);
    CySolverResult* integration_result_ptr = nullptr;

    while (true)
    {
        calculated_surf_pressure = TidalPyConstants::d_INF;
        pass_results_vec.clear();

        if (!final_pass)
        {
            iterations++;
        }

        for (size_t segment_i = 0; segment_i < num_segments; segment_i++)
        {
            const c_EOSSegment& segment = segment_vec[segment_i];
            radius_stop = segment.upper_radius;
            if (segment_i == 0)
            {
                radius_start = 0.0;
                for (size_t y_i = 0; y_i < num_y; y_i++)
                {
                    y0_bysegment_vec[y_i] = y0[y_i];
                }
            }
            else if (integrate_temperature)
            {
                // A segment that sets its own base temperature breaks the profile there, an isothermal
                // layer against its neighbor; the rest continue from the segment below.
                if (std::isfinite(segment.start_temperature))
                {
                    y0_bysegment_vec[4] = segment.start_temperature;
                }
                y0_bysegment_vec[5] = segment.start_heat_flow;
            }

            // Maximum step of one third of the segment thickness.
            max_step = 0.33 * (radius_stop - radius_start);

            eos_input_layer_ptr = &eos_input_bylayer_vec[segment.layer_index];
            eos_input_layer_ptr->temperature_kind = segment.temperature_kind;
            eos_input_layer_ptr->conduction_coeff = segment.conduction_coeff;
            eos_input_layer_ptr->adiabat_coeff    = segment.adiabat_coeff;

            // The integration needs only the density, so skip the moduli, viscosity, and melt models.
            eos_input_layer_ptr->update_bulk  = false;
            eos_input_layer_ptr->update_shear = false;
            use_dense_output = capture_dense;

            std::memcpy(args_vec_eos_ptr, eos_input_layer_ptr, sizeof(c_EOS_ODEInput));

            layer_eos_func = eos_function_bylayer_ptr_vec[segment.layer_index];

            if (!integration_result_uptr)
            {
                integration_result_uptr = std::make_unique<CySolverResult>(integration_method);
            }
            integration_result_ptr = integration_result_uptr.get();

            baseline_cysolve_ivp_noreturn(
                integration_result_ptr,
                diffeq,            // Differential equation [DiffeqFuncType]
                radius_start,      // Start radius for this segment
                radius_stop,       // Stop radius for this segment
                y0_bysegment_vec,  // y0 array vector<double>
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

            if (!failed)
            {
                top_of_last_segment_index = (num_extra + num_y) * (last_solution_size - 1);
                if (segment_i == num_segments - 1)
                {
                    // Pressure is the second state variable of the last step.
                    calculated_surf_pressure = integration_result_ptr->solution[top_of_last_segment_index + 1];
                }
                else
                {
                    // The next segment starts from the top of this one.
                    radius_start = segment.upper_radius;
                    std::memcpy(
                        y0_bysegment_vec.data(),
                        &integration_result_ptr->solution[top_of_last_segment_index],
                        sizeof(double) * num_y);
                }
                if (capture_dense)
                {
                    // Held until the pass is judged; a converged pass hands these to the solution.
                    pass_results_vec.push_back(std::move(integration_result_uptr));
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

        pressure_diff     = calculated_surf_pressure - surface_pressure;
        pressure_diff_abs = std::fabs(pressure_diff);
        const bool converged = (pressure_diff_abs <= pressure_tol_abs);

        if (final_pass || (converged && capture_dense))
        {
            // The whole result of every segment is kept for later dense calls.
            for (std::unique_ptr<CySolverResult>& result_uptr : pass_results_vec)
            {
                eos_solution_ptr->save_cyresult(std::move(result_uptr));
            }
            pass_results_vec.clear();
            break;
        }
        else
        {
            if (converged)
            {
                // Converged without its dense output, so run the same central pressure again and keep it.
                final_pass    = true;
                capture_dense = true;
            }
            else
            {
                // Record which side of the root this pass landed on. Two false-position passes in a row on one side
                // leave the other end stale, so its mismatch is halved (Illinois), which keeps the step moving.
                const double current_central = y0[1];
                if (pressure_diff < 0.0)
                {
                    if (last_bracket_step && (last_replaced == -1)) { above_diff *= 0.5; }
                    below_central = current_central;
                    below_diff    = pressure_diff;
                    last_replaced = -1;
                }
                else
                {
                    if (last_bracket_step && (last_replaced == +1)) { below_diff *= 0.5; }
                    above_central = current_central;
                    above_diff    = pressure_diff;
                    last_replaced = +1;
                }
                const bool bracketed = std::isfinite(below_central) && std::isfinite(above_central);

                // The secant step on the measured slope when it is positive (and, once bracketed, stays inside
                // the bracket). Otherwise a false-position step inside the bracket, or before there is one, a
                // unit-slope step that doubles on every pass the slope stays unusable, so a mismatch that first
                // moves the wrong way costs a few passes rather than one per residual's worth of pressure.
                double step      = TidalPyConstants::d_NAN;
                bool   secant_ok = false;
                last_bracket_step = false;
                if (std::isfinite(previous_central))
                {
                    const double slope = (pressure_diff - previous_diff) / (current_central - previous_central);
                    if (std::isfinite(slope) && (slope > 0.0))
                    {
                        step      = -pressure_diff / slope;
                        secant_ok = std::isfinite(step);
                    }
                }
                if (bracketed)
                {
                    const double bracket_lo = std::min(below_central, above_central);
                    const double bracket_hi = std::max(below_central, above_central);
                    double next_central = current_central + step;
                    if (!(secant_ok && (next_central > bracket_lo) && (next_central < bracket_hi)))
                    {
                        next_central = above_central
                            - above_diff * (above_central - below_central) / (above_diff - below_diff);
                        last_bracket_step = true;
                    }
                    // A converging step is under half the one two passes back. One that is not (a mismatch with a
                    // near-jump, from an EOS at the end of its range, stalls both the secant and false position)
                    // is replaced by a bisection, which shrinks the bracket steadily.
                    const bool stalled = std::fabs(next_central - current_central) > 0.5 * step_two_ago_abs;
                    if (stalled || !(std::isfinite(next_central) && (next_central > bracket_lo)
                                     && (next_central < bracket_hi)))
                    {
                        next_central      = 0.5 * (bracket_lo + bracket_hi);
                        last_bracket_step = true;
                    }
                    step = next_central - current_central;
                }
                else if (!secant_ok)
                {
                    const double direction = (pressure_diff > 0.0) ? -1.0 : 1.0;
                    step = direction * std::max(std::fabs(pressure_diff), 2.0 * last_step_abs);
                }
                else if ((last_step_abs > 0.0) && (std::fabs(step) > 4.0 * last_step_abs))
                {
                    // Before the root is bracketed a nearly flat mismatch can send the secant orders of magnitude
                    // past it; growing at most fourfold per pass still reaches a distant root geometrically.
                    step = std::copysign(4.0 * last_step_abs, step);
                }
                // The secant error model e[k+1] = M e[k] e[k-1], with M measured from the residuals in
                // hand, says whether the next pass should converge and so whether to capture its dense
                // output. With one residual there is no model yet, and a wrong guess costs more than a
                // repeated pass saves, so that pass runs without it.
                const double reference_diff = std::isfinite(oldest_diff) ? oldest_diff : previous_diff;
                const double predicted_diff = std::isfinite(reference_diff)
                    ? pressure_diff * pressure_diff / std::fabs(reference_diff) : TidalPyConstants::d_INF;
                capture_dense = (predicted_diff <= pressure_tol_abs);

                oldest_diff      = previous_diff;
                previous_central = y0[1];
                previous_diff    = pressure_diff;

                // A non-finite step has nowhere to go, and this pass kept no dense output to report.
                if (!std::isfinite(step))
                {
                    failed = true;
                    integrator_failure_message =
                        std::string("the secant step on the central pressure is not finite (surface pressure "
                                    "residual ") + std::to_string(pressure_diff) + std::string(" Pa)");
                    break;
                }
                // Keep the central pressure positive by halving an overshooting step.
                double next_central = y0[1] + step;
                while ((next_central <= 0.0) && (std::fabs(step) > TidalPyConstants::d_EPS * pressure_scale))
                {
                    step        *= 0.5;
                    next_central = y0[1] + step;
                }
                // Before a second step there is no step two passes back to bound the next one.
                step_two_ago_abs = (last_step_abs > 0.0) ? last_step_abs : TidalPyConstants::d_INF;
                last_step_abs    = std::fabs(step);
                y0[1] = next_central;
            }
        }

        if (!final_pass && (iterations >= static_cast<int>(max_iters)))
        {
            max_iters_hit = true;
            eos_solution_ptr->max_iters_hit = true;
            // Still produce output: one more pass, kept whatever it finds.
            final_pass    = true;
            capture_dense = true;
        }
    }

    eos_solution_ptr->iterations = iterations;

    // The pass the cap forces is kept for its diagnostics, but a structure whose surface pressure misses the
    // target is not hydrostatic and must not be reported as solved.
    const bool unconverged = max_iters_hit && !failed && !(pressure_diff_abs <= pressure_tol_abs);

    if (failed)
    {
        eos_solution_ptr->success = false;
        if (iterations > 1)
        {
            // A later pass fails where the search for the central pressure has taken it, usually far past any
            // structure the layers can hold.
            eos_solution_ptr->message =
                std::string("`c_solve_eos` found no hydrostatic structure: the structure integration failed at "
                            "iteration ") + std::to_string(iterations) + std::string(", at a central pressure of ") +
                std::to_string(y0[1]) + std::string(" (") + std::to_string(y0[1] / pressure_scale) +
                std::string(" times the uniform-sphere estimate), while searching for the central pressure that "
                            "meets the target surface pressure. The layers' equations of state may have no "
                            "hydrostatic solution at this radius and mass");
        }
        else
        {
            eos_solution_ptr->message =
                std::string("Warning in `c_solve_eos`: Integrator failed at iteration ") + std::to_string(iterations);
        }
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
        eos_solution_ptr->success = !unconverged;
        eos_solution_ptr->pressure_error = pressure_diff_abs;
        if (unconverged)
        {
            eos_solution_ptr->message =
                std::string("`c_solve_eos` found no hydrostatic structure: after ") + std::to_string(iterations) +
                std::string(" iterations the surface pressure misses its target by ") +
                std::to_string(pressure_diff_abs) + std::string(" Pa (tolerance ") +
                std::to_string(pressure_tol_abs) + std::string(" Pa). The layers' equations of state may have no "
                "hydrostatic solution at this radius and mass; otherwise raise `max_iters`.");
            if (verbose)
            {
                std::printf("%s\n", eos_solution_ptr->message.c_str());
            }
        }

        // Keep the layer EOS functions for on-demand evaluation, then sample onto the radius array.
        eos_solution_ptr->save_eos_functions(eos_function_bylayer_ptr_vec, eos_input_bylayer_vec);
        eos_solution_ptr->interpolate_full_planet();
    }

    if (integration_result_uptr)
    {
        integration_result_uptr.reset();
    }

    integration_result_ptr = nullptr;
}
