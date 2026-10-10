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
#include "layer_kind_.hpp"     // c_layer_layout
#include "start_radius_.hpp"   // c_auto_starting_radius, c_find_start_layer
#include "rs_solution_.hpp"
#include "love_.hpp"
#include "starting/common_.hpp"
#include "starting/driver_.hpp"
#include "interfaces/interfaces_.hpp"
#include "derivatives/odes_.hpp"
#include "orthonormalize_.hpp"
#include "boundaries/boundaries_.hpp"
#include "boundaries/surface_bc_.hpp"

// Material imports
#include "../Material/eos/eos_solution_.hpp"


constexpr double d_EPS_DBL_10000 = 10000.0 * TidalPyConstants::d_EPS;
// Radii across a dynamic liquid layer at which its stratification is read to bound its buoyancy frequency.
constexpr size_t C_STRATIFICATION_SAMPLES = 9;
// Bytes a kept integration segment holds beyond the arrays counted for it: the allocators' bookkeeping for its dozen
// or so allocations and the objects they point to. Measured as the private memory of earth_simple's dynamic core at
// 1000 and 3000 days less the arrays, about 3 KB a segment.
constexpr size_t C_KEPT_SEGMENT_OVERHEAD_BYTES = 3072;


/// Values each dense interpolant keeps per y: y itself and the method's interpolation data (CySolverDense), which is
/// polynomial coefficients for RK23, RK45, and DOP853, step history for the implicit methods, and barycentric node
/// values for Tsit5, Vern7, and Vern8.
inline size_t c_dense_values_per_y(ODEMethod method) noexcept
{
    switch (method)
    {
        case ODEMethod::RK23:   return 1 + RK23_len_Pcols;
        case ODEMethod::RK45:   return 1 + RK45_len_Pcols;
        case ODEMethod::DOP853: return 1 + DOP853_INTERPOLATOR_POWER;
        case ODEMethod::TSIT5:  return 1 + Tsit5_len_Pcols;
        case ODEMethod::VERN7:  return 1 + Vern7_len_Pcols;
        case ODEMethod::VERN8:  return 1 + Vern8_len_Pcols;
        case ODEMethod::BDF:    return 1 + BDF_MAX_ORDER;
        case ODEMethod::RADAU:  return 1 + RADAU_INTERPOLATOR_POWER;
        default:                return 1 + LSODA_MAX_ORDER;
    }
}

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
    int       starting_method     = 0;            // c_StartingMethod (starting/starting_method_.hpp)
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
    // Keep only what the Love numbers need: no dense output, so the radial functions below the surface are not
    // available afterwards (c_RadialSolutionStorage::p_love_only).
    bool      love_only           = false;
};


/// The y (`num_y` reals) of an integration's last stored step, when that step ended at `radius` to the layer
/// continuity tolerance: the top of a layer integrated without dense output. False, NaN-filled, otherwise.
inline bool c_last_step_at(const CySolverResult* result_ptr, double radius, double* y_out_ptr, size_t num_y) noexcept
{
    for (size_t y_i = 0; y_i < num_y; ++y_i) { y_out_ptr[y_i] = TidalPyConstants::d_NAN; }
    if (!result_ptr || !result_ptr->config_uptr || (result_ptr->size == 0) || (result_ptr->num_y < num_y)) {
        return false;
    }
    const size_t stride = result_ptr->config_uptr->capture_extra ? result_ptr->num_dy : result_ptr->num_y;
    const size_t last_i = result_ptr->size - 1;
    if ((result_ptr->time_domain_vec.size() <= last_i) || (result_ptr->solution.size() < (last_i + 1) * stride)) {
        return false;
    }
    const double end_rtol = tidalpy_config_ptr ? tidalpy_config_ptr->d_LAYER_CONTINUITY_RTOL : 0.0;
    if (!(std::abs(result_ptr->time_domain_vec[last_i] - radius) <= end_rtol * std::abs(radius))) { return false; }
    const double* last_y = result_ptr->solution.data() + last_i * stride;
    for (size_t y_i = 0; y_i < num_y; ++y_i) { y_out_ptr[y_i] = last_y[y_i]; }
    return true;
}


// What the collapse reads of one integrated layer, recorded as the integration finishes it so the surface-down
// pass repeats no dense or material evaluation.
struct c_LayerCollapseInputs {
    // Each independent solution's y at the top of the layer (the layer's last segment), [solution * C_MAX_NUM_Y + y];
    // unused entries NaN. top_y holds TS72's ys; top_y_integrated holds them as integrated, with a dynamic liquid's P
    // in y2's slot.
    std::array<std::complex<double>, C_MAX_NUM_SOL * C_MAX_NUM_Y> top_y;
    std::array<std::complex<double>, C_MAX_NUM_SOL * C_MAX_NUM_Y> top_y_integrated;

    // Material at the layer's lower bound (the starting radius in the starting layer) and at its top slice.
    double gravity_lower = TidalPyConstants::d_NAN;
    double density_lower = TidalPyConstants::d_NAN;
    double gravity_upper = TidalPyConstants::d_NAN;
    double density_upper = TidalPyConstants::d_NAN;

    size_t num_sols          = 0;
    int    layer_type        = 0;
    bool   is_static         = false;
    bool   is_incompressible = false;

    // The map from this layer's solutions to the layer above's, recorded by the upward pass at the interface on top
    // of this layer (c_solve_upper_y_at_interface); the collapse applies it downward (c_collapse_through_interface).
    std::array<std::complex<double>, C_INTERFACE_TRANSFER_SIZE> transfer_to_above{};
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
    solution_storage_ptr->surface_frame_residual = TidalPyConstants::d_NAN;
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

    double boundary_conditions[C_MAX_SURFACE_CONDITIONS];
    double* bc_pointer = &boundary_conditions[0];
    const int surface_bc_code = c_get_surface_bc(
        bc_pointer, inputs.bc_models.data(), num_ytypes, planet_radius, inputs.planet_bulk_density, degree_l_dbl);
    if (surface_bc_code != 0)
    {
        // A bad model would otherwise leave NaN boundary conditions for the surface solve to collapse onto.
        return solution_storage_ptr->fail(
            -14, "RadialSolver.ShootingMethod:: " + c_surface_bc_error_message(surface_bc_code), verbose);
    }

    const size_t num_extra            = 0;
    // A layer's solutions are integrated as one system, and every segment's dense CyRK result is retained so the
    // collapsed y can be evaluated at any radius; an empty t_eval keeps the adaptive dense segments with no fixed grid.
    // A Love-only solve needs the solutions at the surface alone, the last step CyRK stores, so it keeps no dense
    // output and reuses one result. With DOP853 that saves three of the fifteen right-hand-side calls of a step only
    // in a layer without the independence event (one solution, or re-orthonormalization off): CyRK builds every
    // step's interpolant to check an event.
    const bool love_only            = inputs.love_only;
    const bool capture_dense_output = !love_only;
    std::vector<double> teval_empty;

    // The factor by which an integration may reduce the normalized Gram determinant of a layer's solutions before they
    // are re-orthonormalized (orthonormalize_.hpp); zero never.
    const double independence_floor = tidalpy_config_ptr->d_MIN_SOLUTION_INDEPENDENCE;
    const bool use_orthonormalization = (independence_floor > 0.0);

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
    solution_storage_ptr->p_love_only         = love_only;
    solution_storage_ptr->p_frequency_solve   = frequency;
    solution_storage_ptr->p_pressure_unit     = c_pressure_unit(frequency, inputs.planet_bulk_density, G_to_use);
    solution_storage_ptr->p_segments_by_layer.resize(num_layers);
    solution_storage_ptr->shooting_method_orthonormalizations_vec.assign(num_layers, 0);
    solution_storage_ptr->p_num_sols_by_layer.assign(num_solutions_by_layer_vec.begin(), num_solutions_by_layer_vec.end());
    solution_storage_ptr->p_layer_types.resize(num_layers);
    solution_storage_ptr->p_layer_is_static.resize(num_layers);
    solution_storage_ptr->p_upper_radii_solve.resize(num_layers);
    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        solution_storage_ptr->p_layer_types[current_layer_i]     = layer_types_ptr[current_layer_i];
        solution_storage_ptr->p_layer_is_static[current_layer_i] = is_static_by_layer_ptr[current_layer_i] ? 1 : 0;
        solution_storage_ptr->p_upper_radii_solve[current_layer_i] =
            eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i];
    }
    // Collapse constants: [ytype][layer][segment][<=3 solutions], filled in the collapse phase.
    solution_storage_ptr->p_constants_by_ytype_layer.assign(
        num_ytypes, std::vector<std::vector<std::array<std::complex<double>, 3>>>(num_layers));

    // Stack buffer sized for the largest case: 6 ys x 3 solutions = 18 complex.
    std::complex<double> initial_y[C_MAX_NUM_SOL * C_MAX_NUM_Y];
    std::complex<double>* initial_y_ptr = &initial_y[0];
    bool starting_y_check = false;

    // Filled layer by layer as the integration finishes each one; layers below the starting layer stay unset.
    std::vector<c_LayerCollapseInputs> collapse_inputs_vec(num_layers);

    // Scratch for the per-radius EOS reads below: gravity, density, and the complex moduli at this frequency.
    c_EOSMaterialState eos_material_state;
    std::unique_ptr<CySolverResult> integration_solution_uptr = std::make_unique<CySolverResult>(integration_method);
    CySolverResult* integration_solution_ptr = nullptr;
    // Bytes of the dense results kept so far, which max_ram_MB bounds over the whole solve.
    size_t retained_bytes = 0;
    const size_t ram_limit_bytes = inputs.max_ram_MB * 1024 * 1024;

    // Extra diffeq arguments beyond y and r.
    std::vector<char> diffeq_args_vec(sizeof(c_RadialSolverArgs));
    c_RadialSolverArgs* diffeq_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(diffeq_args_vec.data());
    PreEvalFunc diffeq_preeval_ptr = nullptr;
    // A layer's solutions as one real vector, [solution * 2 * num_ys + real/imag pair], at a segment's start and top.
    std::vector<double> y0_vec(C_MAX_NUM_SOL * C_MAX_NUM_Y_REAL);
    std::vector<double> segment_top_vec(C_MAX_NUM_SOL * C_MAX_NUM_Y_REAL);

    diffeq_args_ptr->degree_l         = degree_l_dbl;
    diffeq_args_ptr->lp1              = degree_l_dbl + 1.0;
    diffeq_args_ptr->lm1              = degree_l_dbl - 1.0;
    diffeq_args_ptr->llp1             = degree_l_dbl * (degree_l_dbl + 1.0);
    diffeq_args_ptr->G                = G_to_use;
    diffeq_args_ptr->grav_coeff       = 4.0 * TidalPyConstants::d_PI * G_to_use;
    diffeq_args_ptr->frequency        = frequency;
    diffeq_args_ptr->pressure_unit    = c_pressure_unit(frequency, inputs.planet_bulk_density, G_to_use);
    diffeq_args_ptr->layer_index      = 0;
    diffeq_args_ptr->eos_solution_ptr = eos_solution_storage_ptr;
    diffeq_args_ptr->num_solutions    = 0;
    diffeq_args_ptr->num_ys           = 0;
    diffeq_args_ptr->independence_floor = independence_floor;

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

    const double starting_radius =
        c_auto_starting_radius(inputs.starting_radius, planet_radius, inputs.start_radius_tol, degree_l_dbl);

    // The layer holding the starting radius (c_find_start_layer), and that layer's first slice at or above the start;
    // lower layers are skipped. The starting layer is integrated from the starting radius to its top slice, so that
    // slice must lie above the start. A starting radius in no layer fails the solve.
    size_t start_layer_i           = 0;
    size_t start_first_slice_index = 0;  // first slice at or above the starting radius
    size_t start_layer_slices      = 0;  // slices from there to the top of the starting layer
    const bool start_layer_found = c_find_start_layer(
        eos_solution_storage_ptr->upper_radius_bylayer_vec, num_layers, starting_radius, start_layer_i);
    if (start_layer_found)
    {
        const size_t layer_first_slice = first_slice_index_by_layer_vec[start_layer_i];
        const size_t layer_end_slice   = layer_first_slice + num_slices_by_layer_vec[start_layer_i];
        size_t slice_i = layer_first_slice;
        while ((slice_i < layer_end_slice) && (radius_array_ptr[slice_i] < starting_radius))
        {
            ++slice_i;
        }
        start_first_slice_index = slice_i;
        start_layer_slices      = layer_end_slice - slice_i;
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
    // Main integration loop: bottom layer up; starting conditions, then the layer's solutions as one ODE system
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
        layer_inputs.top_y_integrated.fill(c_constant_NAN);
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

        // A compressible solid or dynamic liquid reads the bulk modulus, and a non-positive one (a Poisson ratio of
        // -1 or below with shear, an infinite compressibility without) solves without complaint to wrong Love
        // numbers. A zero is what a material that names no bulk modulus carries, so it is refused at both ends of
        // the layer. A static liquid never reads it.
        const bool layer_reads_bulk = (!layer_is_incomp) && ((layer_type == 0) || !layer_is_static);
        if (layer_reads_bulk &&
            (!(std::real(bulk_lower) > 0.0) || !(std::real(eos_material_state.bulk_modulus) > 0.0)))
        {
            solution_storage_ptr->error_code = -15;
            solution_storage_ptr->success    = false;
            solution_storage_ptr->message    =
                std::string("RadialSolver.ShootingMethod:: Layer ") + std::to_string(current_layer_i) +
                std::string(" is compressible but its bulk modulus is not positive (real part ") +
                c_format_scientific(std::real(bulk_lower)) + std::string(" at its base, ") +
                c_format_scientific(std::real(eos_material_state.bulk_modulus)) +
                std::string(" at its top, in solve units). Give the layer's material a bulk modulus, or mark ") +
                std::string("the layer incompressible (is_incompressible True).\n");
            if (verbose)
            {
                printf("%s", solution_storage_ptr->message.c_str());
            }
            return solution_storage_ptr->error_code;
        }

        if (max_step_from_arrays)
        {
            max_step_to_use = std::abs(0.33 * (radius_upper - radius_lower));
            max_step_to_use = std::fmax(max_step_to_use, d_EPS_DBL_10000);
        }

        if (inputs.scale_rtols)
        {
            // Two reals per complex y, the same pattern for each solution of the layer.
            rtols_vec.resize(num_sols * num_ys_dbl);
            atols_vec.resize(num_sols * num_ys_dbl);

            const c_LayerKindLayout& layout = c_layer_layout(layer_type, layer_is_static);
            for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
            {
                for (size_t y_i = 0; y_i < num_ys; ++y_i)
                {
                    // Tighter rtols on the ys that drive instability (c_LayerKindLayout::rtol_scale).
                    const size_t value_i = solution_i * num_ys_dbl + 2 * y_i;
                    rtols_vec[value_i]     = integration_rtol * layout.rtol_scale[y_i];
                    rtols_vec[value_i + 1] = integration_rtol * layout.rtol_scale[y_i];
                    atols_vec[value_i]     = integration_atol;
                    atols_vec[value_i + 1] = integration_atol;
                }
            }
        }

        // NaN so uninitialized reads are visible (could be dropped for a small performance gain).
        for (size_t y_i = 0; y_i < C_MAX_NUM_SOL * C_MAX_NUM_Y; ++y_i)
        {
            initial_y_ptr[y_i] = c_constant_NAN;
        }

        if (is_start_layer)
        {
            c_find_starting_conditions(
                &solution_storage_ptr->success,
                solution_storage_ptr->message,
                layer_type,
                layer_is_static,
                layer_is_incomp,
                inputs.starting_method,
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
                layer_type,
                layer_is_static,
                interface_values.gravity,
                interface_values.liquid_density,
                G_to_use,
                layer_below.transfer_to_above.data()
            );
        }

        // A dynamic liquid is integrated in P = y2 - rho g y1 + rho y5 (derivatives/odes_.hpp), from this layer's own
        // material at its base.
        const bool pressure_form = (layer_type != 0) && !layer_is_static;
        if (pressure_form)
        {
            for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
            {
                std::complex<double>* start_y = &initial_y_ptr[solution_i * C_MAX_NUM_Y];
                start_y[1] =
                    c_pressure_variable_from_y2(start_y[0], start_y[1], start_y[2], gravity_lower, density_lower,
                                                diffeq_args_ptr->pressure_unit);
            }
        }

        // The integrator works on real pairs, every solution of the layer in one system.
        y0_vec.resize(num_sols * num_ys_dbl);
        segment_top_vec.resize(num_sols * num_ys_dbl);
        for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
        {
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
            {
                const std::complex<double> dcomplex_tmp = initial_y_ptr[solution_i * C_MAX_NUM_Y + y_i];
                y0_vec[solution_i * num_ys_dbl + 2 * y_i]     = dcomplex_tmp.real();
                y0_vec[solution_i * num_ys_dbl + 2 * y_i + 1] = dcomplex_tmp.imag();
            }
        }

        DiffeqFuncType layer_diffeq = c_find_layer_diffeq(layer_type, layer_is_static, layer_is_incomp);

        diffeq_args_ptr->layer_index   = current_layer_i;
        diffeq_args_ptr->num_solutions = num_sols;
        diffeq_args_ptr->num_ys        = num_ys;
        diffeq_args_ptr->y_weights     =
            c_layer_y_weights(layer_type, layer_is_static, planet_radius, inputs.planet_bulk_density, G_to_use);
        // A stratified liquid's solutions are waves (N^2 > 0) or growth (N^2 < 0) on the scale omega / N, which carry P
        // (in its unit) about N / omega times y1, so the independence takes P at that larger size. Its largest |N^2|
        // across the layer sets it; a neutral liquid's is 0 and needs no change.
        if (pressure_form)
        {
            double buoyancy_bound = 0.0;
            for (size_t sample_i = 0; sample_i < C_STRATIFICATION_SAMPLES; ++sample_i)
            {
                const double sample_radius = radius_lower + (radius_upper - radius_lower)
                    * static_cast<double>(sample_i) / static_cast<double>(C_STRATIFICATION_SAMPLES - 1);
                const c_RadialMaterial material = c_read_eos(diffeq_args_ptr, sample_radius);
                const std::complex<double> stratification = layer_is_incomp ?
                    c_liquid_stratification<true>(material) : c_liquid_stratification<false>(material);
                const double buoyancy = std::abs(material.gravity * stratification.real() / material.density);
                if (std::isfinite(buoyancy)) { buoyancy_bound = std::max(buoyancy_bound, buoyancy); }
            }
            if (buoyancy_bound > frequency * frequency)
            {
                diffeq_args_ptr->y_weights[1] *= frequency / std::sqrt(buoyancy_bound);
            }
        }
        const double* y_weights_ptr = diffeq_args_ptr->y_weights.data();

        // The layer's solutions are re-orthonormalized where the integration makes them nearly dependent. A terminal
        // event ends a segment where their independence falls by the floor's factor from where the segment started
        // (to the floor itself after a re-orthonormalization, which starts at 1), and the next segment starts from an
        // orthonormal basis of their span (orthonormalize_.hpp). The solutions a layer starts from are otherwise taken
        // as given, since the starting conditions near the center can be nearly dependent in the leading powers of r,
        // and orthonormalizing them only makes the integrator resolve what the surface solve does not need. The event's
        // threshold must stay well above the roundoff of the independence (about machine epsilon), or the root finder
        // chases noise, so solutions too dependent for that are orthonormalized first.
        const bool orthonormalize = use_orthonormalization && (num_sols > 1);
        c_BasisChange basis_change = c_identity_basis_change();
        size_t& layer_orthonormalizations =
            solution_storage_ptr->shooting_method_orthonormalizations_vec[current_layer_i];
        if (pressure_form && (num_sols == 2))
        {
            c_split_pressure_solutions(y0_vec.data(), num_ys, y_weights_ptr, basis_change);
        }
        double start_independence = c_solution_independence(y0_vec.data(), num_sols, num_ys, y_weights_ptr);
        if (orthonormalize && !(independence_floor * std::fmin(1.0, start_independence) > d_EPS_DBL_10000))
        {
            c_BasisChange start_change = c_identity_basis_change();
            c_orthonormalize_solutions(y0_vec.data(), num_sols, num_ys, y_weights_ptr, start_change);
            basis_change = c_compose_basis_changes(basis_change, start_change, num_sols);
            ++layer_orthonormalizations;
            start_independence = 1.0;
        }
        std::vector<Event> layer_events_vec;
        if (orthonormalize)
        {
            layer_events_vec.emplace_back(c_solution_independence_event, 1, -1);
        }
        diffeq_args_ptr->independence_floor = independence_floor * std::fmin(1.0, start_independence);
        // Every solution shares the layer's steps, recorded as each segment finishes so a failure reports them too.
        const auto record_layer_steps = [&](size_t layer_steps)
        {
            for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
            {
                solution_storage_ptr->shooting_method_steps_taken_vec[(3 * current_layer_i) + solution_i] = layer_steps;
            }
        };

        std::vector<c_ShootingSegment>& layer_segments = solution_storage_ptr->p_segments_by_layer[current_layer_i];
        layer_segments.clear();
        double segment_lower   = radius_lower;
        double first_step_size = 0.0;   // found by CyRK for a layer's first segment
        size_t layer_steps     = 0;
        size_t expected_size   = inputs.expected_size;   // a restart takes about its predecessor's steps
        while (true)
        {
            if (!integration_solution_uptr)
            {
                integration_solution_uptr = std::make_unique<CySolverResult>(integration_method);
            }
            integration_solution_ptr = integration_solution_uptr.get();

            // The step cap holds for the whole layer, whatever its segments. Zero asks for CyRK's own, the steps
            // max_ram_MB can store.
            const size_t layer_step_cap = (inputs.max_num_steps > 0) ? inputs.max_num_steps :
                find_max_num_steps(num_sols * num_ys_dbl, num_extra, 0, inputs.max_ram_MB).max_num_steps;
            if (layer_steps >= layer_step_cap)
            {
                return solution_storage_ptr->fail(
                    -11,
                    std::string("RadialSolver.ShootingMethod:: Integration problem at layer ") +
                        std::to_string(current_layer_i) + std::string(": its ") + std::to_string(layer_steps) +
                        std::string(" steps over ") + std::to_string(layer_segments.size()) +
                        std::string(" re-orthonormalized segments reached max_num_steps.\n"),
                    verbose);
            }
            const size_t steps_left = layer_step_cap - layer_steps;

            baseline_cysolve_ivp_noreturn(
                integration_solution_ptr,  // Raw pointer to CySolverResult structure.
                layer_diffeq,              // Differential equation [DiffeqFuncType]
                segment_lower,             // Start radius
                radius_upper,              // End radius
                y0_vec,                    // y0 array vector[double]
                expected_size,             // Expected final integration size (0 = find good value) [size_t]
                num_extra,                 // Number of extra parameters captured during integration [size_t]
                diffeq_args_vec,           // Extra input args to diffeq vector[char]
                steps_left,                // Max number of steps (0 = find good value) [size_t]
                inputs.max_ram_MB,         // Max amount of RAM allowed [size_t]
                capture_dense_output,      // Use dense output [bool]
                teval_empty,               // No fixed eval grid; the dense interpolant is retained instead
                diffeq_preeval_ptr,        // Pre-eval function used in diffeq [PreEvalFunc]
                layer_events_vec,          // Events vector[Event]
                rtols_vec,                 // Relative Tolerance (as array) vector[double]
                atols_vec,                 // Absolute Tolerance (as array) vector[double]
                max_step_to_use,           // Maximum step size [double]
                first_step_size,           // Initial step size (0 = find good value) [double]
                false,                     // Force retain solver: its interpolants do not need it [bool]
                nullptr                    // Analytic jacobian (null = numerical; used by implicit methods) [JacobianFuncType]
            );
            layer_steps += integration_solution_ptr->steps_taken;
            record_layer_steps(layer_steps);

            if (!integration_solution_ptr->success)
            {
                solution_storage_ptr->error_code = -11;
                solution_storage_ptr->success    = false;
                solution_storage_ptr->message    =
                    std::string("RadialSolver.ShootingMethod:: Integration problem at layer ") +
                    std::to_string(current_layer_i) + std::string(":\n\t") + integration_solution_ptr->message +
                    std::string("\n");

                if (verbose)
                {
                    printf("%s", solution_storage_ptr->message.c_str());
                }
                return solution_storage_ptr->error_code;
            }

            // The segment ends at the top of the layer, or where the event stopped it; CyRK stores the last step there.
            // An event within roundoff of the top has reached the top. Any wider margin would read a fast-growing
            // solution short of the top as its value there.
            const size_t num_stored     = integration_solution_ptr->size;
            const double stopped_radius = integration_solution_ptr->time_domain_vec[num_stored - 1];
            const double top_slack      = TidalPyConstants::d_EPS_100 * std::abs(radius_upper);
            const bool   restart        =
                integration_solution_ptr->event_terminated && ((radius_upper - stopped_radius) > top_slack);
            const double segment_upper  = restart ? stopped_radius : radius_upper;
            if (!c_last_step_at(integration_solution_ptr, segment_upper, segment_top_vec.data(), num_sols * num_ys_dbl))
            {
                return solution_storage_ptr->fail(
                    -11,
                    std::string("RadialSolver.ShootingMethod:: No solution at the top of layer ") +
                        std::to_string(current_layer_i) + std::string(".\n"),
                    verbose);
            }
            if (restart && !(segment_upper > segment_lower))
            {
                return solution_storage_ptr->fail(
                    -11,
                    std::string("RadialSolver.ShootingMethod:: The solutions of layer ") +
                        std::to_string(current_layer_i) + std::string(" lost their independence faster than the ") +
                        std::string("integrator could step at radius ") + c_format_scientific(segment_lower) +
                        std::string(" (solve units). Lower [numerical] minimum_solution_independence.\n"),
                    verbose);
            }

            c_ShootingSegment segment;
            segment.radius_lower = segment_lower;
            segment.radius_upper = segment_upper;
            segment.basis_change = basis_change;
            // The next segment takes the last whole step the stopped one accepted, rather than a fresh guess.
            first_step_size = 0.0;
            if (restart && (num_stored >= 3))
            {
                const double last_step = integration_solution_ptr->time_domain_vec[num_stored - 2] -
                    integration_solution_ptr->time_domain_vec[num_stored - 3];
                const double step_cap  = std::fmin(max_step_to_use, 0.5 * (radius_upper - segment_upper));
                first_step_size = (last_step > 0.0) ? std::fmin(last_step, step_cap) : 0.0;
            }
            // Retain the dense interpolant, shrunk to what it holds, and arm a fresh result for the next integration; a
            // Love-only solve reuses its one result.
            expected_size = std::max<size_t>(8, 2 * num_stored);
            if (!love_only)
            {
                CySolverResult& kept = *integration_solution_ptr;
                kept.time_domain_vec.shrink_to_fit();
                kept.solution.shrink_to_fit();
                kept.time_domain_vec_sorted.shrink_to_fit();
                kept.interp_time_vec.shrink_to_fit();
                kept.dense_vec.shrink_to_fit();
                // Each interpolant holds y and the method's interpolation coefficients, and the problem's
                // configuration stays with the result.
                const ProblemConfig& config = *kept.config_uptr;
                const size_t dense_doubles = c_dense_values_per_y(kept.integrator_method) * kept.num_y;
                retained_bytes += sizeof(CySolverResult) + C_KEPT_SEGMENT_OVERHEAD_BYTES
                    + sizeof(double) * (kept.time_domain_vec.capacity()
                    + kept.solution.capacity() + kept.time_domain_vec_sorted.capacity()
                    + kept.interp_time_vec.capacity() + config.y0_vec.capacity() + config.rtols.capacity()
                    + config.atols.capacity()) + config.args_vec.capacity()
                    + config.events_vec.capacity() * sizeof(Event)
                    + kept.dense_vec.capacity() * (sizeof(CySolverDense) + sizeof(double) * dense_doubles);
                segment.result = std::move(integration_solution_uptr);
                if ((ram_limit_bytes > 0) && (retained_bytes > ram_limit_bytes))
                {
                    layer_segments.push_back(std::move(segment));
                    return solution_storage_ptr->fail(
                        -11,
                        std::string("RadialSolver.ShootingMethod:: The dense solution kept through layer ") +
                            std::to_string(current_layer_i) + std::string(" (") +
                            std::to_string(layer_segments.size()) +
                            std::string(" re-orthonormalized segments in this layer) exceeds max_ram_MB. ") +
                            std::string("Raise max_ram_MB, or solve with love_only.\n"),
                        verbose);
                }
            }
            layer_segments.push_back(std::move(segment));
            if (!restart)
            {
                break;
            }

            y0_vec = segment_top_vec;
            c_orthonormalize_solutions(y0_vec.data(), num_sols, num_ys, y_weights_ptr, basis_change);
            ++layer_orthonormalizations;
            diffeq_args_ptr->independence_floor = independence_floor;
            segment_lower = segment_upper;
        }

        // Top-of-layer y, for the next layer's interface condition and for the collapse; a dynamic liquid's P back to
        // y2 from its material at the top.
        for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
        {
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
            {
                const size_t value_i = solution_i * num_ys_dbl + 2 * y_i;
                layer_inputs.top_y_integrated[solution_i * C_MAX_NUM_Y + y_i] =
                    std::complex<double>(segment_top_vec[value_i], segment_top_vec[value_i + 1]);
            }
        }
        layer_inputs.top_y = layer_inputs.top_y_integrated;
        if (pressure_form)
        {
            for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
            {
                std::complex<double>* top = &layer_inputs.top_y[solution_i * C_MAX_NUM_Y];
                top[1] = c_y2_from_pressure_variable(
                    top[0], top[1], top[2], layer_inputs.gravity_upper, layer_inputs.density_upper,
                    diffeq_args_ptr->pressure_unit);
            }
        }
    }
    // The surface layer's solutions at the surface as integrated, which a Love-only solve's find_love collapses.
    if (love_only)
    {
        solution_storage_ptr->p_surface_top_y.assign(
            collapse_inputs_vec[num_layers - 1].top_y_integrated.begin(),
            collapse_inputs_vec[num_layers - 1].top_y_integrated.end());
    }

    // =================================================================================================================
    // Collapse: surface boundary conditions, then the interface conditions from the surface down to the start
    // =================================================================================================================
    solution_storage_ptr->message = std::string("Integration completed for all layers. Beginning solution collapse.\n");

    c_LayerCollapseInputs& surface_layer = collapse_inputs_vec[num_layers - 1];

    // A rigid translation of the whole body (y1 = y3 = u, y2 = y4 = 0, y5 = g u, y6 = 0) solves the static
    // equations at degree 1 for any structure and meets every surface condition, so the degree-1 Love numbers are
    // fixed only once a reference frame is chosen (Farrell 1972; Blewitt 2003). A loading solve fixes the frame of the
    // body's own center of mass (k' = 0) in the surface solve itself (c_apply_surface_bc), static or dynamic: with
    // inertia the exact solution is already in that frame (Saito 1974 App. 2), but the system it would be found from
    // is singular as the frequency squared. The solution is shifted to another frame afterwards
    // (c_RadialSolutionStorage::find_love). Tidal and free solves have no degree-1 answer to fix.
    bool degree1_frame = (degree_l == 1);
    for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
    {
        degree1_frame = degree1_frame && (inputs.bc_models[ytype_i] == 2);
    }

    // Rank of the surface system, which does not depend on the boundary condition: a singular system still hands
    // back finite, arbitrary constants that the amplification cannot flag.
    const double surface_rcond = c_estimate_surface_rcond(
        surface_layer.top_y.data(),
        surface_layer.num_sols,
        2 * surface_layer.num_sols,
        C_MAX_NUM_Y,
        surface_layer.layer_type,
        surface_layer.is_static,
        degree1_frame);
    solution_storage_ptr->surface_rcond = surface_rcond;

    // Without a frame row, the static degree-1 system is singular in exact arithmetic. Integration error can hold the
    // computed rcond far above machine precision (1e-11 at rtol 1e-6 below a static liquid core), so this is decided
    // from the structure rather than from the rcond threshold below.
    if ((degree_l == 1) && !degree1_frame)
    {
        bool all_layers_static = true;
        for (size_t current_layer_i = start_layer_i; current_layer_i < num_layers; ++current_layer_i)
        {
            if (!is_static_by_layer_ptr[current_layer_i])
            {
                all_layers_static = false;
                break;
            }
        }
        if (all_layers_static)
        {
            solution_storage_ptr->error_code = -13;
            solution_storage_ptr->success    = false;
            solution_storage_ptr->message    =
                std::string("RadialSolver.ShootingMethod:: A degree-1 tidal or free solve in which every ") +
                std::string("integrated layer is static is singular: a rigid translation of the body satisfies the ") +
                std::string("equations and every surface condition. Only loading has degree-1 Love numbers; solve ") +
                std::string("for loading alone, which fixes the reference frame.\n");
            if (verbose)
            {
                printf("%s", solution_storage_ptr->message.c_str());
            }
            return solution_storage_ptr->error_code;
        }
    }

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
            std::string("), so its solution constants are undetermined. Keep [numerical] ") +
            std::string("minimum_solution_independence above 0, which re-orthonormalizes nearly dependent ") +
            std::string("solutions; a system that stays singular is singular by construction.\n");
        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }
        return solution_storage_ptr->error_code;
    }

    // A static liquid surface layer integrates only y5 and y7, but at the free surface y2 is the boundary
    // condition and y2 = rho (g y1 - y5) (S74 Eq. 20), so y1, and with it h, is defined there; y3, and l, are not.
    solution_storage_ptr->p_static_surface_y2.clear();
    if ((surface_layer.layer_type != 0) && surface_layer.is_static)
    {
        solution_storage_ptr->p_static_surface_density = surface_layer.density_upper;
        solution_storage_ptr->p_static_surface_gravity = surface_gravity;
        solution_storage_ptr->p_static_surface_y2.resize(num_ytypes);
        for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
        {
            solution_storage_ptr->p_static_surface_y2[ytype_i] =
                std::complex<double>(bc_pointer[ytype_i * 3 + 0], 0.0);
        }
    }

    for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
    {
        solution_storage_ptr->message =
            std::string("Collapsing radial solutions for \"") +
            std::to_string(ytype_i) +
            std::string("\" solver.\n");

        bc_solution_info = -999;
        double frame_residual = TidalPyConstants::d_NAN;
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
                    layer.is_incompressible,
                    degree1_frame,
                    &frame_residual
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

                // How far the condition the frame row replaced is from met, worst across ytypes (NaN unless
                // degree1_frame). A static body meets it to roundoff; inertia in some layers but not others leaves
                // a residual of about omega^2 R / g.
                if (degree1_frame)
                {
                    const double previous = solution_storage_ptr->surface_frame_residual;
                    solution_storage_ptr->surface_frame_residual =
                        std::isnan(previous) ? frame_residual : std::fmax(previous, frame_residual);
                }

                // Worst-case error amplification of the surface solve across ytypes (large cancelling constants
                // amplify roundoff); cheap, so always recorded, and the wrappers decide whether to warn.
                const std::array<double, 6> surface_weights = c_layer_y_weights(
                    layer.layer_type, layer.is_static, planet_radius, inputs.planet_bulk_density, G_to_use);
                solution_storage_ptr->surface_amplification = std::fmax(
                    solution_storage_ptr->surface_amplification,
                    c_estimate_surface_amplification(
                        constant_vector_ptr,
                        layer.top_y.data(),
                        num_sols,
                        2 * num_sols,
                        C_MAX_NUM_Y,
                        surface_weights.data()));
            }
            else
            {
                // Interior layer: constants follow from the layer above through the map the upward pass recorded.
                const c_LayerCollapseInputs& layer_above = collapse_inputs_vec[layer_i + 1];
                c_collapse_through_interface(
                    constant_vector_ptr,
                    layer_above_constant_vector_ptr,
                    layer.transfer_to_above.data(),
                    num_sols,
                    layer_above.num_sols);
            }

            // The constants hold for the layer's last segment. Each segment keeps its own, and those of the solutions
            // it continues are R^-1 times them (orthonormalize_.hpp), down to the layer's starting solutions, which
            // the interface below maps from. The collapsed y is evaluated on demand from the interpolants and these
            // constants (c_RadialSolutionStorage::get_radial_solution_nondim); unused entries stay NaN.
            const std::vector<c_ShootingSegment>& layer_segments = solution_storage_ptr->p_segments_by_layer[layer_i];
            std::vector<std::array<std::complex<double>, 3>>& dest_constants =
                solution_storage_ptr->p_constants_by_ytype_layer[ytype_i][layer_i];
            dest_constants.assign(
                layer_segments.size(),
                std::array<std::complex<double>, 3>{c_constant_NAN, c_constant_NAN, c_constant_NAN});
            for (size_t segment_i = layer_segments.size(); segment_i-- > 0;)
            {
                for (size_t s = 0; s < num_sols; ++s)
                    dest_constants[segment_i][s] = constant_vector_ptr[s];
                c_constants_before_basis_change(
                    layer_segments[segment_i].basis_change, constant_vector_ptr, constant_vector_ptr, num_sols);
            }

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
