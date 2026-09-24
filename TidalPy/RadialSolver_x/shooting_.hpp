#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <memory>
#include <vector>
#include <complex>

// CyRK imports
#include "cysolution.hpp"
#include "c_events.hpp"
#include "cysolve.hpp"

// TidalPy imports
#include "../constants_.hpp"

// RadialSolver imports
#include "rs_constants_.hpp"
#include "rs_solution_.hpp"
#include "love_.hpp"
#include "starting/common_.hpp"
#include "starting/driver_.hpp"
#include "interfaces/interfaces_.hpp"
#include "interfaces/reversed_.hpp"
#include "derivatives/odes_.hpp"
#include "boundaries/boundaries_.hpp"
#include "boundaries/surface_bc_.hpp"

// Material imports
#include "../Material_x/eos/eos_solution_.hpp"


constexpr double d_EPS_DBL_10000 = 10000.0 * TidalPyConstants::d_EPS;

// A value in scientific notation for a status message; std::to_string prints fixed point, which shows a small
// radius or condition number as zero.
inline std::string c_format_scientific(double value)
{
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%.3e", value);
    return std::string(buffer);
}


// Inputs of the shooting method. The layer structure is set once per EOS solve; the rest are per-call knobs.
struct c_ShootingInputs {
    // layer_types: 0 = solid, 1 = liquid. bool[] because std::vector<bool> is bit-packed and the solver
    // wants a bool*.
    std::vector<int>        layer_types;
    std::unique_ptr<bool[]> is_static;
    std::unique_ptr<bool[]> is_incompressible;
    size_t                  num_layers = 0;

    // Per-layer slice partitioning over the (non-dim) radius grid.
    std::vector<size_t> first_slice_index_by_layer;
    std::vector<size_t> num_slices_by_layer;

    // Surface boundary conditions to solve for, in order; one block of radial functions per entry. The
    // independent solutions do not depend on the boundary condition, so n conditions cost one integration
    // and n surface solves rather than n integrations.
    std::vector<int> bc_models = {1};

    // Non-dim planet scalars.
    double planet_bulk_density = 0.0;
    double G                   = 0.0;
    int    degree_l            = 2;

    // Shooting-method knobs (per-call, set from the runtime config).
    bool      use_kamata          = false;
    double    starting_radius     = 0.0;          // non-dim
    double    start_radius_tol    = 1.0e-4;
    ODEMethod integration_method  = ODEMethod::DOP853;
    double    integration_rtol    = 1.0e-5;
    double    integration_atol    = 1.0e-7;
    bool      scale_rtols         = true;
    size_t    max_num_steps       = 500000;
    size_t    expected_size       = 500;
    size_t    max_ram_MB          = 500;
    double    max_step            = 0.0;
};


// What the collapse reads of one integrated layer, recorded as the integration finishes it so the surface-down
// pass repeats no dense or material evaluation.
struct c_LayerCollapseInputs {
    // Each independent solution's y at the top of the layer, [solution * C_MAX_NUM_Y + y]; unused entries NaN.
    std::array<std::complex<double>, C_MAX_NUM_SOL * C_MAX_NUM_Y> top_y;

    // Material at the layer's lower bound (the starting radius in the starting layer) and at its top slice.
    double gravity_lower = TidalPyConstants::d_NAN;
    double density_lower = TidalPyConstants::d_NAN;
    double gravity_upper = TidalPyConstants::d_NAN;
    double density_upper = TidalPyConstants::d_NAN;

    size_t num_sols          = 0;
    int    layer_type        = 0;
    bool   is_static         = false;
    bool   is_incompressible = false;
};


/// Solve the viscoelastic-gravitational problem with the shooting method.
inline int c_shooting_solver(
    c_RadialSolutionStorage* solution_storage_ptr,
    const c_ShootingInputs& inputs,
    double frequency,
    bool verbose) noexcept
{
    c_EOSSolution* eos_solution_storage_ptr = solution_storage_ptr->get_eos_solution_ptr();

    solution_storage_ptr->message = std::string("RadialSolver.ShootingMethod:: Starting integration\n");
    // The conditioning diagnostics describe this solve only; they stay at these values if it stops early.
    solution_storage_ptr->surface_amplification = 0.0;
    solution_storage_ptr->surface_rcond         = TidalPyConstants::d_NAN;
    if (verbose)
    {
        printf("%s", solution_storage_ptr->message.c_str());
    }

    const int*  layer_types_ptr                = inputs.layer_types.data();
    const bool* is_static_by_layer_ptr         = inputs.is_static.get();
    const bool* is_incompressible_by_layer_ptr = inputs.is_incompressible.get();
    const std::vector<size_t>& first_slice_index_by_layer_vec = inputs.first_slice_index_by_layer;
    const std::vector<size_t>& num_slices_by_layer_vec        = inputs.num_slices_by_layer;

    const double G_to_use        = inputs.G;
    const int    degree_l        = inputs.degree_l;
    const double integration_rtol = inputs.integration_rtol;
    const double integration_atol = inputs.integration_atol;
    const ODEMethod integration_method = inputs.integration_method;

    const double degree_l_dbl = static_cast<double>(degree_l);

    double* radius_array_ptr  = eos_solution_storage_ptr->radius_array_vec.data();

    const size_t num_layers   = eos_solution_storage_ptr->num_layers;

    const double planet_radius   = eos_solution_storage_ptr->radius;
    const double surface_gravity = eos_solution_storage_ptr->surface_gravity;

    // Surface boundary conditions per forcing type; tides follow (y2, y4, y6) = (0, 0, (2l+1)/R).
    const size_t num_ytypes = inputs.bc_models.size();

    // 15 = 5 (max solve_for entries) * 3 (surface conditions)
    double boundary_conditions[15];
    double* bc_pointer = &boundary_conditions[0];
    const int surface_bc_code = c_get_surface_bc(
        bc_pointer,
        inputs.bc_models.data(),
        num_ytypes,
        planet_radius,
        inputs.planet_bulk_density,
        degree_l_dbl
    );
    if (surface_bc_code != 0)
    {
        // A bad model would otherwise leave NaN boundary conditions for the surface solve to collapse onto.
        solution_storage_ptr->error_code = -14;
        solution_storage_ptr->success    = false;
        solution_storage_ptr->message    =
            std::string("RadialSolver.ShootingMethod:: Invalid surface boundary conditions (code ") +
            std::to_string(surface_bc_code) +
            std::string("): between 1 and 5 models are allowed, each free (0), tidal (1), or loading (2).\n");
        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }
        return solution_storage_ptr->error_code;
    }

    const size_t num_extra            = 0;
    const double first_step_size      = 0.0;
    // Every independent solution's dense CyRK result is retained so the collapsed y can be evaluated at any radius;
    // an empty t_eval keeps the adaptive dense segments with no fixed grid.
    const bool   capture_dense_output = true;
    std::vector<double> teval_empty;
    std::vector<Event> events_vec;  // unused

    // A zero max_step means one third of each layer's thickness (below).
    double max_step_to_use      = TidalPyConstants::d_NAN;
    bool   max_step_from_arrays = false;

    if (inputs.max_step == 0.0)
    {
        max_step_from_arrays = true;
    }
    else
    {
        max_step_to_use = inputs.max_step;
    }

    std::vector<double> rtols_vec{integration_rtol};
    std::vector<double> atols_vec{integration_atol};

    std::vector<size_t> num_solutions_by_layer_vec;
    num_solutions_by_layer_vec.resize(num_layers);
    size_t* num_solutions_by_layer_ptr = num_solutions_by_layer_vec.data();

    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        num_solutions_by_layer_ptr[current_layer_i] = c_find_num_shooting_solutions(
            layer_types_ptr[current_layer_i],
            is_static_by_layer_ptr[current_layer_i],
            is_incompressible_by_layer_ptr[current_layer_i]
        );
    }

    // Storage for the per-(layer, solution) dense interpolants, the collapse constants, and the per-layer metadata
    // needed to combine them at any radius.
    const std::complex<double> c_constant_NAN(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
    solution_storage_ptr->reset_interpolant_storage();
    solution_storage_ptr->p_uses_interpolants = true;
    solution_storage_ptr->p_frequency_solve   = frequency;
    solution_storage_ptr->p_interp_by_layer_sol.resize(num_layers);
    solution_storage_ptr->p_num_sols_by_layer.assign(num_solutions_by_layer_vec.begin(), num_solutions_by_layer_vec.end());
    solution_storage_ptr->p_layer_types.resize(num_layers);
    solution_storage_ptr->p_layer_is_static.resize(num_layers);
    solution_storage_ptr->p_upper_radii_solve.resize(num_layers);
    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        const size_t num_sols = num_solutions_by_layer_ptr[current_layer_i];
        solution_storage_ptr->p_interp_by_layer_sol[current_layer_i].resize(num_sols);
        solution_storage_ptr->p_layer_types[current_layer_i]     = layer_types_ptr[current_layer_i];
        solution_storage_ptr->p_layer_is_static[current_layer_i] = is_static_by_layer_ptr[current_layer_i] ? 1 : 0;
        solution_storage_ptr->p_upper_radii_solve[current_layer_i] =
            eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i];
    }
    // Collapse constants: [ytype][layer][<=3 solutions], filled NaN until the collapse phase.
    solution_storage_ptr->p_constants_by_ytype_layer.assign(
        num_ytypes,
        std::vector<std::array<std::complex<double>, 3>>(
            num_layers, std::array<std::complex<double>, 3>{c_constant_NAN, c_constant_NAN, c_constant_NAN})
        );

    // Stack buffers sized for the largest case: 6 ys x 3 solutions = 18 complex, 36 real.
    std::complex<double> initial_y[18];
    std::complex<double>* initial_y_ptr = &initial_y[0];
    double initial_y_only_real[36];
    double* initial_y_only_real_ptr = &initial_y_only_real[0];
    bool starting_y_check = false;

    // Filled layer by layer as the integration finishes each one; layers below the starting layer stay unset.
    std::vector<c_LayerCollapseInputs> collapse_inputs_vec(num_layers);

    // Scratch for the per-radius EOS reads below: gravity, density, and the complex moduli at this frequency.
    c_EOSMaterialState eos_material_state;
    std::unique_ptr<CySolverResult> integration_solution_uptr = std::make_unique<CySolverResult>(integration_method);
    CySolverResult* integration_solution_ptr = nullptr;

    // Extra diffeq arguments beyond y and r.
    std::vector<char> diffeq_args_vec(sizeof(c_RadialSolverArgs));
    c_RadialSolverArgs* diffeq_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(diffeq_args_vec.data());
    PreEvalFunc diffeq_preeval_ptr = nullptr;
    std::vector<double> y0_vec(C_MAX_NUM_Y_REAL);

    diffeq_args_ptr->degree_l         = degree_l_dbl;
    diffeq_args_ptr->lp1              = degree_l_dbl + 1.0;
    diffeq_args_ptr->lm1              = degree_l_dbl - 1.0;
    diffeq_args_ptr->llp1             = degree_l_dbl * (degree_l_dbl + 1.0);
    diffeq_args_ptr->G                = G_to_use;
    diffeq_args_ptr->grav_coeff       = 4.0 * TidalPyConstants::d_PI * G_to_use;
    diffeq_args_ptr->frequency        = frequency;
    diffeq_args_ptr->layer_index      = 0;
    diffeq_args_ptr->eos_solution_ptr = eos_solution_storage_ptr;

    // At most 3 solutions per layer, so the constant vectors live on the stack.
    std::complex<double> constant_vector[3];
    std::complex<double>* constant_vector_ptr = &constant_vector[0];
    std::complex<double> layer_above_constant_vector[3];
    std::complex<double>* layer_above_constant_vector_ptr = &layer_above_constant_vector[0];

    for (size_t i = 0; i < 3; ++i)
    {
        constant_vector_ptr[i]             = c_constant_NAN;
        layer_above_constant_vector_ptr[i] = c_constant_NAN;
    }

    // Surface linear-solve status; -999 means not yet called.
    int bc_solution_info = -999;

    // The integration cannot start at r = 0 (singularity); starting higher is more stable but skips more of the
    // planet. The automatic choice follows Martens' thesis and the LoadDef manual, capped by
    // config_x [numerical].max_start_radius_fraction.
    double starting_radius = inputs.starting_radius;
    if (starting_radius == 0.0)
    {
        starting_radius = planet_radius * std::pow(inputs.start_radius_tol, 1.0 / degree_l_dbl);
        starting_radius = std::fmin(starting_radius, tidalpy_config_ptr->d_MAX_START_RADIUS_FRAC * planet_radius);
    }

    // Find the layer holding the starting radius, [lower, upper), so a start on an interface begins in the layer
    // above it, and that layer's first slice at or above the start; lower layers are skipped. The starting layer
    // is integrated from the starting radius to its top slice, so that slice must lie above the start. A starting
    // radius at or below zero, at or above the surface, or NaN lies in no layer and fails the solve.
    bool   start_layer_found       = false;
    size_t start_layer_i           = 0;
    size_t start_first_slice_index = 0;  // first slice at or above the starting radius
    size_t start_layer_slices      = 0;  // slices from there to the top of the starting layer
    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        const double layer_upper_radius = eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i];
        const double layer_lower_radius = (current_layer_i == 0) ?
            0.0 : eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i - 1];

        if ((starting_radius > 0.0) && (layer_lower_radius <= starting_radius) &&
            (starting_radius < layer_upper_radius))
        {
            start_layer_found = true;
            start_layer_i     = current_layer_i;

            const size_t layer_first_slice = first_slice_index_by_layer_vec[current_layer_i];
            const size_t layer_end_slice   = layer_first_slice + num_slices_by_layer_vec[current_layer_i];
            size_t slice_i = layer_first_slice;
            while ((slice_i < layer_end_slice) && (radius_array_ptr[slice_i] < starting_radius))
            {
                ++slice_i;
            }
            start_first_slice_index = slice_i;
            start_layer_slices      = layer_end_slice - slice_i;
            break;
        }
    }

    if (!start_layer_found)
    {
        solution_storage_ptr->error_code = -5;
        solution_storage_ptr->success    = false;
        solution_storage_ptr->message    =
            std::string("RadialSolver.ShootingMethod:: The starting radius (") + c_format_scientific(starting_radius) +
            std::string(" in solve units) is not inside the planet, (0, ") + c_format_scientific(planet_radius) +
            std::string("). Use a starting radius inside the planet, or 0 for the automatic choice.\n");
        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }
        return solution_storage_ptr->error_code;
    }
    if ((start_layer_slices == 0) ||
        !(radius_array_ptr[start_first_slice_index + start_layer_slices - 1] > starting_radius))
    {
        solution_storage_ptr->error_code = -5;
        solution_storage_ptr->success    = false;
        solution_storage_ptr->message    =
            std::string("RadialSolver.ShootingMethod:: No radial slice of layer ") + std::to_string(start_layer_i) +
            std::string(" lies above the starting radius (") + c_format_scientific(starting_radius) +
            std::string(" in solve units), so there is nothing to integrate through. Lower the starting radius.\n");
        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }
        return solution_storage_ptr->error_code;
    }

    // Record the start info so the dense calling system NaNs any query below the starting radius.
    solution_storage_ptr->p_start_layer_i         = start_layer_i;
    solution_storage_ptr->p_starting_radius_solve = starting_radius;

    // =================================================================================================================
    // Main integration loop: bottom layer up; starting conditions, then the ODE for each independent solution
    // =================================================================================================================
    for (size_t current_layer_i = start_layer_i; current_layer_i < num_layers; ++current_layer_i)
    {
        // The starting layer begins at the slice at or above the starting radius.
        const bool   is_start_layer    = (current_layer_i == start_layer_i);
        const size_t first_slice_index =
            is_start_layer ? start_first_slice_index : first_slice_index_by_layer_vec[current_layer_i];
        const size_t layer_slices      =
            is_start_layer ? start_layer_slices : num_slices_by_layer_vec[current_layer_i];

        const size_t num_sols   = num_solutions_by_layer_ptr[current_layer_i];
        const size_t num_ys     = 2 * num_sols;
        const size_t num_ys_dbl = 2 * num_ys;
        double* layer_radius_ptr = &radius_array_ptr[first_slice_index];

        const int layer_type       = layer_types_ptr[current_layer_i];
        const bool layer_is_static = is_static_by_layer_ptr[current_layer_i];
        const bool layer_is_incomp = is_incompressible_by_layer_ptr[current_layer_i];

        c_LayerCollapseInputs& layer_inputs = collapse_inputs_vec[current_layer_i];
        layer_inputs.top_y.fill(c_constant_NAN);
        layer_inputs.num_sols          = num_sols;
        layer_inputs.layer_type        = layer_type;
        layer_inputs.is_static         = layer_is_static;
        layer_inputs.is_incompressible = layer_is_incomp;

        // The starting radius is generally not on a stored slice, so the EOS is evaluated there. Every layer asks
        // the solution for its base rather than reading the slice arrays: the two agree at a slice radius, but
        // going through call_material means the material-state provider is honoured, so a world solve takes its
        // interface values from the same models the integration uses.
        const double radius_lower = is_start_layer ? starting_radius : layer_radius_ptr[0];
        eos_solution_storage_ptr->call_material(current_layer_i, radius_lower, eos_material_state);
        const double gravity_lower               = eos_material_state.gravity;
        const double density_lower               = eos_material_state.density;
        const std::complex<double> shear_lower   = eos_material_state.shear_modulus;
        const std::complex<double> bulk_lower    = eos_material_state.bulk_modulus;
        layer_inputs.gravity_lower = gravity_lower;
        layer_inputs.density_lower = density_lower;

        const double radius_upper = layer_radius_ptr[layer_slices - 1];
        eos_solution_storage_ptr->call_material(current_layer_i, radius_upper, eos_material_state);
        layer_inputs.gravity_upper = eos_material_state.gravity;
        layer_inputs.density_upper = eos_material_state.density;

        if (max_step_from_arrays)
        {
            max_step_to_use = std::abs(0.33 * (radius_upper - radius_lower));
            max_step_to_use = std::fmax(max_step_to_use, d_EPS_DBL_10000);
        }

        if (inputs.scale_rtols)
        {
            // Two reals per complex y.
            rtols_vec.resize(num_ys * 2);
            atols_vec.resize(num_ys * 2);

            for (size_t y_i = 0; y_i < num_ys; ++y_i)
            {
                // TODO: Change up the tolerance scaling between real and imaginary?
                double layer_rtol_real = integration_rtol;
                double layer_rtol_imag = integration_rtol;
                double layer_atol_real = integration_atol;
                double layer_atol_imag = integration_atol;

                // Tighter rtols on the ys that drive instability. TODO: test these scales (Issue #44).
                if (layer_type == 0)
                {
                    // Solid: y2 and y3.
                    if ((y_i == 1) || (y_i == 2))
                    {
                        layer_rtol_real *= 0.1;
                        layer_rtol_imag *= 0.1;
                    }
                }
                else
                {
                    if (!layer_is_static)
                    {
                        // Dynamic liquid: y2.
                        if (y_i == 1)
                        {
                            layer_rtol_real *= 0.01;
                            layer_rtol_imag *= 0.01;
                        }
                    }
                }
                rtols_vec[2 * y_i]     = layer_rtol_real;
                rtols_vec[2 * y_i + 1] = layer_rtol_imag;
                atols_vec[2 * y_i]     = layer_atol_real;
                atols_vec[2 * y_i + 1] = layer_atol_imag;
            }
        }

        // NaN so uninitialized reads are visible (could be dropped for a small performance gain).
        for (size_t y_i = 0; y_i < 36; ++y_i)
        {
            if (y_i < 18)
            {
                initial_y_ptr[y_i] = c_constant_NAN;
            }
            initial_y_only_real_ptr[y_i] = TidalPyConstants::d_NAN;
        }

        if (is_start_layer)
        {
            c_find_starting_conditions(
                &solution_storage_ptr->success,
                solution_storage_ptr->message,
                layer_type,
                layer_is_static,
                layer_is_incomp,
                inputs.use_kamata,
                frequency,
                radius_lower,
                density_lower,
                bulk_lower,
                shear_lower,
                degree_l,
                G_to_use,
                C_MAX_NUM_Y,
                initial_y_ptr,
                starting_y_check
            );

            if (!solution_storage_ptr->success)
            {
                solution_storage_ptr->error_code = -10;
                if (verbose)
                {
                    printf("%s", solution_storage_ptr->message.c_str());
                }
                return solution_storage_ptr->error_code;
            }
        }
        else
        {
            c_LayerCollapseInputs& layer_below = collapse_inputs_vec[current_layer_i - 1];
            const c_InterfaceValues interface_values = c_interface_values(
                c_InterfaceSide{
                    layer_below.layer_type, layer_below.is_static, layer_below.gravity_upper,
                    layer_below.density_upper},
                c_InterfaceSide{layer_type, layer_is_static, gravity_lower, density_lower});

            c_solve_upper_y_at_interface(
                layer_below.top_y.data(),
                initial_y_ptr,
                layer_below.num_sols,
                num_sols,
                C_MAX_NUM_Y,
                layer_below.layer_type,
                layer_below.is_static,
                layer_below.is_incompressible,
                layer_type,
                layer_is_static,
                layer_is_incomp,
                interface_values.gravity,
                interface_values.liquid_density,
                G_to_use
            );
        }

        // The integrator works on real pairs.
        for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
        {
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
            {
                const std::complex<double> dcomplex_tmp = initial_y_ptr[solution_i * C_MAX_NUM_Y + y_i];
                initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL + 2 * y_i]     = dcomplex_tmp.real();
                initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL + 2 * y_i + 1] = dcomplex_tmp.imag();
            }
        }

        DiffeqFuncType layer_diffeq = c_find_layer_diffeq(layer_type, layer_is_static, layer_is_incomp);

        diffeq_args_ptr->layer_index = current_layer_i;

        for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
        {
            if (!integration_solution_uptr)
            {
                integration_solution_uptr = std::make_unique<CySolverResult>(integration_method);
            }
            integration_solution_ptr = integration_solution_uptr.get();

            y0_vec.resize(num_ys_dbl);
            std::memcpy(y0_vec.data(), &initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL], sizeof(double) * num_ys_dbl);

            baseline_cysolve_ivp_noreturn(
                integration_solution_ptr,  // Raw pointer to CySolverResult structure.
                layer_diffeq,              // Differential equation [DiffeqFuncType]
                radius_lower,              // Start radius
                radius_upper,              // End radius
                y0_vec,                    // y0 array vector[double]
                inputs.expected_size,      // Expected final integration size (0 = find good value) [size_t]
                num_extra,                 // Number of extra parameters captured during integration [size_t]
                diffeq_args_vec,           // Extra input args to diffeq vector[char]
                inputs.max_num_steps,      // Max number of steps (0 = find good value) [size_t]
                inputs.max_ram_MB,         // Max amount of RAM allowed [size_t]
                capture_dense_output,      // Use dense output [bool]
                teval_empty,               // No fixed eval grid; the dense interpolant is retained instead
                diffeq_preeval_ptr,        // Pre-eval function used in diffeq [PreEvalFunc]
                events_vec,                // Events vector[Event]
                rtols_vec,                 // Relative Tolerance (as array) vector[double]
                atols_vec,                 // Absolute Tolerance (as array) vector[double]
                max_step_to_use,           // Maximum step size [double]
                first_step_size,           // Initial step size (0 = find good value) [double]
                true,                      // Force retain solver [bool]
                nullptr                    // Analytic jacobian (null = numerical; used by implicit methods) [JacobianFuncType]
            );

            solution_storage_ptr->shooting_method_steps_taken_vec[(3 * current_layer_i) + solution_i] =
                integration_solution_ptr->steps_taken;

            if (!integration_solution_ptr->success)
            {
                solution_storage_ptr->error_code = -11;
                solution_storage_ptr->success    = false;
                solution_storage_ptr->message    =
                    std::string("RadialSolver.ShootingMethod:: Integration problem at layer ") +
                    std::to_string(current_layer_i) + std::string("; solution ") + std::to_string(solution_i) +
                    std::string(":\n\t") + integration_solution_ptr->message + std::string("\n");

                if (verbose)
                {
                    printf("%s", solution_storage_ptr->message.c_str());
                }
                return solution_storage_ptr->error_code;
            }

            // Top-of-layer y, for the next layer's interface condition and for the collapse.
            double interp_top[C_MAX_NUM_Y_REAL];
            if (!c_call_dense_checked(integration_solution_ptr, radius_upper, interp_top, 2 * num_ys))
            {
                solution_storage_ptr->error_code = -11;
                solution_storage_ptr->success    = false;
                solution_storage_ptr->message    =
                    std::string("RadialSolver.ShootingMethod:: Dense output unavailable at the top of layer ") +
                    std::to_string(current_layer_i) + std::string("; solution ") + std::to_string(solution_i) +
                    std::string(".\n");
                if (verbose)
                {
                    printf("%s", solution_storage_ptr->message.c_str());
                }
                return solution_storage_ptr->error_code;
            }
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
            {
                layer_inputs.top_y[solution_i * C_MAX_NUM_Y + y_i] =
                    std::complex<double>(interp_top[2 * y_i], interp_top[2 * y_i + 1]);
            }

            // Retain the dense interpolant and arm a fresh result for the next integration.
            solution_storage_ptr->p_interp_by_layer_sol[current_layer_i][solution_i] =
                std::move(integration_solution_uptr);
            integration_solution_uptr = std::make_unique<CySolverResult>(integration_method);
            integration_solution_ptr  = integration_solution_uptr.get();
        }
    }

    // =================================================================================================================
    // Collapse: surface boundary conditions, then the interface conditions from the surface down to the start
    // =================================================================================================================
    solution_storage_ptr->message = std::string("Integration completed for all layers. Beginning solution collapse.\n");

    c_LayerCollapseInputs& surface_layer = collapse_inputs_vec[num_layers - 1];

    // Rank of the surface system, which does not depend on the boundary condition: a singular system still hands
    // back finite, arbitrary constants that the amplification cannot flag.
    const double surface_rcond = c_estimate_surface_rcond(
        surface_layer.top_y.data(),
        surface_layer.num_sols,
        2 * surface_layer.num_sols,
        C_MAX_NUM_Y,
        surface_layer.layer_type,
        surface_layer.is_static);
    solution_storage_ptr->surface_rcond = surface_rcond;
    const double minimum_surface_rcond = tidalpy_config_ptr->d_MIN_SURFACE_RCOND;
    // The negated comparison also fails a NaN rcond; a NaN threshold (config unloaded) never fails.
    if (!(surface_rcond >= minimum_surface_rcond) && !std::isnan(minimum_surface_rcond))
    {
        solution_storage_ptr->error_code = -13;
        solution_storage_ptr->success    = false;
        solution_storage_ptr->message    =
            std::string("RadialSolver.ShootingMethod:: The surface boundary condition system is ") +
            std::string("singular to working precision (reciprocal condition number ") +
            c_format_scientific(surface_rcond) + std::string(" < [numerical] ") +
            std::string("minimum_surface_rcond ") + c_format_scientific(minimum_surface_rcond) +
            std::string("), so its solution constants are undetermined. A degree-1 solve for a ") +
            std::string("static body has a rigid-translation mode that no surface condition fixes; ") +
            std::string("otherwise try a larger or the automatic starting radius.\n");
        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }
        return solution_storage_ptr->error_code;
    }

    for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
    {
        solution_storage_ptr->message =
            std::string("Collapsing radial solutions for \"") +
            std::to_string(ytype_i) +
            std::string("\" solver.\n");

        bc_solution_info = -999;
        for (size_t i = 0; i < 3; ++i)
        {
            constant_vector_ptr[i] = c_constant_NAN;
        }

        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }

        // From the surface down to the starting layer.
        for (size_t layer_i = num_layers; layer_i-- > start_layer_i;)
        {
            c_LayerCollapseInputs& layer = collapse_inputs_vec[layer_i];
            const size_t num_sols = layer.num_sols;

            if (layer_i == num_layers - 1)
            {
                c_apply_surface_bc(
                    constant_vector_ptr,
                    &bc_solution_info,
                    bc_pointer,
                    layer.top_y.data(),
                    surface_gravity,
                    G_to_use,
                    num_sols,
                    C_MAX_NUM_Y,
                    ytype_i,
                    layer.layer_type,
                    layer.is_static,
                    layer.is_incompressible
                );

                if (bc_solution_info != 0)
                {
                    solution_storage_ptr->error_code = -12;
                    solution_storage_ptr->success    = false;
                    solution_storage_ptr->message    =
                        std::string(
                            "RadialSolver.ShootingMethod:: Error encountered while applying surface "
                            "boundary condition. Eigen LU decomp code: ") +
                            std::to_string(bc_solution_info) +
                            std::string("\nThe solutions may not be valid at the surface.\n");

                    if (verbose)
                    {
                        printf("%s", solution_storage_ptr->message.c_str());
                    }
                    return solution_storage_ptr->error_code;
                }

                // Worst-case error amplification of the surface solve across ytypes (large cancelling constants
                // amplify roundoff); cheap, so always recorded, and the wrappers decide whether to warn.
                solution_storage_ptr->surface_amplification = std::fmax(
                    solution_storage_ptr->surface_amplification,
                    c_estimate_surface_amplification(
                        constant_vector_ptr,
                        layer.top_y.data(),
                        num_sols,
                        2 * num_sols,
                        C_MAX_NUM_Y));
            }
            else
            {
                // Interior layer: constants follow from the layer above.
                const c_LayerCollapseInputs& layer_above = collapse_inputs_vec[layer_i + 1];
                c_top_to_bottom_interface_bc(
                    constant_vector_ptr,
                    layer_above_constant_vector_ptr,
                    layer.top_y.data(),
                    layer.gravity_upper,
                    layer_above.gravity_lower,
                    layer.density_upper,
                    layer_above.density_lower,
                    layer.layer_type,
                    layer_above.layer_type,
                    layer.is_static,
                    layer_above.is_static,
                    layer.is_incompressible,
                    layer_above.is_incompressible,
                    num_sols,
                    C_MAX_NUM_Y
                );
            }

            // The collapsed y is evaluated on demand from the interpolants and these constants
            // (c_RadialSolutionStorage::get_radial_solution_nondim); unused entries stay NaN.
            std::array<std::complex<double>, 3>& dest_constants =
                solution_storage_ptr->p_constants_by_ytype_layer[ytype_i][layer_i];
            for (size_t s = 0; s < num_sols; ++s)
                dest_constants[s] = constant_vector_ptr[s];

            for (size_t solution_i = 0; solution_i < 3; ++solution_i)
            {
                layer_above_constant_vector_ptr[solution_i] =
                    (solution_i < num_sols) ? constant_vector_ptr[solution_i] : c_constant_NAN;
            }
        }
    }

    if (solution_storage_ptr->error_code != 0)
    {
        solution_storage_ptr->success = false;
    }
    else
    {
        solution_storage_ptr->success = true;
        solution_storage_ptr->message = std::string("RadialSolver.ShootingMethod: Completed without any noted issues.");
    }

    return solution_storage_ptr->error_code;
}
