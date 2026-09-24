#pragma once

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
#include "constants_.hpp"

// RadialSolver imports
#include "rs_constants_.hpp"
#include "rs_solution_.hpp"
#include "love_.hpp"
#include "starting/common_.hpp"
#include "starting/driver_.hpp"
#include "interfaces/interfaces_.hpp"
#include "interfaces/reversed_.hpp"
#include "derivatives/odes_.hpp"
#include "collapse/collapse_.hpp"
#include "boundaries/boundaries_.hpp"
#include "boundaries/surface_bc_.hpp"

// Material imports
#include "../Material_x/eos/eos_solution_.hpp"


constexpr double d_EPS_DBL_10000 = 10000.0 * TidalPyConstants::d_EPS;

inline size_t c_find_num_shooting_solutions(
    int layer_type,
    bool is_static,
    bool is_incompressible
) noexcept
{
    /** Number of independent shooting solutions a layer needs: 3 solid, 2 dynamic liquid, 1 static liquid.
    layer_type is 0 for solid and 1 for liquid. */
    size_t num_sols = 0;

    if (layer_type == 0)
    {
        if (is_static)
        {
            if (is_incompressible)
            {
                // TODO: Confirm
                num_sols = 3;
            }
            else
            {
                num_sols = 3;
            }
        }
        else
        {
            if (is_incompressible)
            {
                // TODO: Confirm
                num_sols = 3;
            }
            else
            {
                num_sols = 3;
            }
        }
    }
    else
    {
        if (is_static)
        {
            if (is_incompressible)
            {
                // TODO: Confirm
                num_sols = 1;
            }
            else
            {
                num_sols = 1;
            }
        }
        else
        {
            if (is_incompressible)
            {
                // TODO: Confirm
                num_sols = 2;
            }
            else
            {
                num_sols = 2;
            }
        }
    }
    return num_sols;
}

int c_shooting_solver(
    c_RadialSolutionStorage* solution_storage_ptr,
    double frequency,
    double planet_bulk_density,
    int* layer_types_ptr,
    bool* is_static_by_layer_ptr,
    bool* is_incompressible_by_layer_ptr,
    std::vector<size_t>& first_slice_index_by_layer_vec,
    std::vector<size_t>& num_slices_by_layer_vec,
    size_t num_bc_models,
    int* bc_models_ptr,
    double G_to_use,
    int degree_l,
    bool use_kamata,
    double starting_radius,
    double start_radius_tolerance,
    ODEMethod integration_method,
    double integration_rtol,
    double integration_atol,
    bool scale_rtols_by_layer_type,
    size_t max_num_steps,
    size_t expected_size,
    size_t max_ram_MB,
    double max_step,
    bool verbose,
    bool warnings
) noexcept
{
    /** Solve the viscoelastic-gravitational problem with the shooting method. */
    c_EOSSolution* eos_solution_storage_ptr = solution_storage_ptr->get_eos_solution_ptr();

    solution_storage_ptr->message = std::string("RadialSolver.ShootingMethod:: Starting integration\n");
    if (verbose)
    {
        printf("%s", solution_storage_ptr->message.c_str());
    }

    const double degree_l_dbl = static_cast<double>(degree_l);

    double* radius_array_ptr  = eos_solution_storage_ptr->radius_array_vec.data();

    const size_t num_layers   = eos_solution_storage_ptr->num_layers;
    const size_t total_slices = eos_solution_storage_ptr->radius_array_size;

    const double planet_radius   = eos_solution_storage_ptr->radius;
    const double surface_gravity = eos_solution_storage_ptr->surface_gravity;

    // Surface boundary conditions per forcing type; tides follow (y2, y4, y6) = (0, 0, (2l+1)/R).
    const size_t num_ytypes = num_bc_models;

    // 15 = 5 (max solve_for entries) * 3 (surface conditions)
    double boundary_conditions[15];
    double* bc_pointer = &boundary_conditions[0];
    c_get_surface_bc(
        bc_pointer,
        bc_models_ptr,
        num_ytypes,
        planet_radius,
        planet_bulk_density,
        degree_l_dbl
    );

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

    if (max_step == 0.0)
    {
        max_step_from_arrays = true;
    }
    else
    {
        max_step_to_use = max_step;
    }

    std::vector<double> rtols_vec{integration_rtol};
    std::vector<double> atols_vec{integration_atol};

    std::vector<size_t> num_solutions_by_layer_vec;
    num_solutions_by_layer_vec.resize(num_layers);
    size_t* num_solutions_by_layer_ptr = num_solutions_by_layer_vec.data();

    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        const int layer_type       = layer_types_ptr[current_layer_i];
        const bool layer_is_static = is_static_by_layer_ptr[current_layer_i];
        const bool layer_is_incomp = is_incompressible_by_layer_ptr[current_layer_i];

        const size_t num_sols = c_find_num_shooting_solutions(
            layer_type,
            layer_is_static,
            layer_is_incomp
        );
        num_solutions_by_layer_ptr[current_layer_i] = num_sols;
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
    solution_storage_ptr->p_layer_is_incomp.resize(num_layers);
    solution_storage_ptr->p_upper_radii_solve.resize(num_layers);
    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        const size_t num_sols = num_solutions_by_layer_ptr[current_layer_i];
        solution_storage_ptr->p_interp_by_layer_sol[current_layer_i].resize(num_sols);
        solution_storage_ptr->p_layer_types[current_layer_i]     = layer_types_ptr[current_layer_i];
        solution_storage_ptr->p_layer_is_static[current_layer_i] = is_static_by_layer_ptr[current_layer_i] ? 1 : 0;
        solution_storage_ptr->p_layer_is_incomp[current_layer_i] = is_incompressible_by_layer_ptr[current_layer_i] ? 1 : 0;
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
    std::complex<double> uppermost_y_per_solution[18];
    std::complex<double>* uppermost_y_per_solution_ptr = &uppermost_y_per_solution[0];

    for (size_t i = 0; i < 18; ++i)
    {
        uppermost_y_per_solution_ptr[i] = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
    }

    std::complex<double> initial_y[18];
    std::complex<double>* initial_y_ptr = &initial_y[0];
    double initial_y_only_real[36];
    double* initial_y_only_real_ptr = &initial_y_only_real[0];
    bool starting_y_check = false;

    std::complex<double>* solution_ptr = reinterpret_cast<std::complex<double>*>(solution_storage_ptr->full_solution_vec.data());

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
    std::complex<double> surface_solutions[6];
    std::complex<double>* surface_solutions_ptr = &surface_solutions[0];

    for (size_t i = 0; i < 6; ++i)
    {
        if (i < 3)
        {
            constant_vector_ptr[i]             = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
            layer_above_constant_vector_ptr[i] = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        }
        surface_solutions_ptr[i] = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
    }

    // Surface linear-solve status; -999 means not yet called.
    int bc_solution_info = -999;

    // The integration cannot start at r = 0 (singularity); starting higher is more stable but skips more of the
    // planet. The automatic choice follows Martens' thesis and the LoadDef manual, capped by
    // config_x [numerical].max_start_radius_fraction.
    if (starting_radius == 0.0)
    {
        starting_radius = planet_radius * std::pow(start_radius_tolerance, 1.0 / degree_l_dbl);
        starting_radius = std::fmin(starting_radius, tidalpy_config_ptr->d_MAX_START_RADIUS_FRAC * planet_radius);
    }

    // Find the layer holding the starting radius; lower layers are skipped.
    size_t start_layer_i           = 0;
    size_t last_index_before_start = 0;
    size_t start_index_in_layer    = 0;
    double last_radius_check       = 0.0;
    double layer_upper_radius      = TidalPyConstants::d_INF;
    double last_layer_upper_radius = TidalPyConstants::d_INF;
    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        layer_upper_radius = eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i];
        if (current_layer_i == 0)
        {
            last_layer_upper_radius = 0.0;
        }
        else
        {
            last_layer_upper_radius = eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i - 1];
        }

        if (last_layer_upper_radius < starting_radius && starting_radius <= layer_upper_radius)
        {
            start_layer_i = current_layer_i;
            const size_t first_slice_index = first_slice_index_by_layer_vec[current_layer_i];

            // Find the last slice before the starting radius.
            start_index_in_layer = 0;
            for (size_t slice_i = first_slice_index; slice_i < first_slice_index + num_slices_by_layer_vec[current_layer_i]; ++slice_i)
            {
                const double radius_check = radius_array_ptr[slice_i];
                if (last_radius_check < starting_radius && starting_radius <= radius_check)
                {
                    if (slice_i == 0)
                    {
                        last_index_before_start = 0;
                    }
                    else
                    {
                        last_index_before_start = slice_i - 1;
                    }
                    break;
                }
                else
                {
                    start_index_in_layer += 1;
                    last_radius_check     = radius_check;
                }
            }
            break;
        }
        else
        {
            last_radius_check = last_layer_upper_radius;
        }
    }

    // Record the start info so the dense calling system NaNs any query below the starting radius.
    solution_storage_ptr->p_start_layer_i         = start_layer_i;
    solution_storage_ptr->p_starting_radius_solve = starting_radius;

    // NaN the gridded solution below the starting radius.
    for (size_t slice_i = 0; slice_i < last_index_before_start + 1; ++slice_i)
    {
        for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
        {
            for (size_t y_i = 0; y_i < C_MAX_NUM_Y; ++y_i)
            {
                solution_ptr[slice_i * C_MAX_NUM_Y * num_ytypes + ytype_i * C_MAX_NUM_Y + y_i] = 
                    std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
            }
        }
    }

    // =================================================================================================================
    // Main integration loop: bottom layer up; starting conditions, then the ODE for each independent solution
    // =================================================================================================================
    size_t layer_below_num_sols  = 0;
    int layer_below_type         = -1;
    bool layer_below_is_static   = false;
    bool layer_below_is_incomp   = false;
    double interface_gravity        = TidalPyConstants::d_NAN;
    double static_liquid_density    = TidalPyConstants::d_NAN;
    double last_layer_upper_gravity = TidalPyConstants::d_NAN;
    double last_layer_upper_density = TidalPyConstants::d_NAN;
    
    // Physical parameters at the starting radius, reused during the collapse.
    double starting_gravity;
    double starting_density;
    std::complex<double> starting_shear{TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    std::complex<double> starting_bulk{TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};

    for (size_t current_layer_i = start_layer_i; current_layer_i < num_layers; ++current_layer_i)
    {
        size_t layer_slices = num_slices_by_layer_vec[current_layer_i];
        size_t first_slice_index;
        if (current_layer_i == start_layer_i)
        {
            // The starting layer begins at the slice at or above the starting radius.
            first_slice_index = last_index_before_start + 1;
            layer_slices -= start_index_in_layer;
        }
        else
        {
            first_slice_index = first_slice_index_by_layer_vec[current_layer_i];
        }

        const size_t num_sols   = num_solutions_by_layer_ptr[current_layer_i];
        const size_t num_ys     = 2 * num_sols;
        const size_t num_ys_dbl = 2 * num_ys;
        double* layer_radius_ptr = &radius_array_ptr[first_slice_index];

        double radius_lower{starting_radius};
        double gravity_lower{TidalPyConstants::d_NAN};
        double density_lower{TidalPyConstants::d_NAN};
        std::complex<double> shear_lower{starting_shear};
        std::complex<double> bulk_lower{starting_bulk};
        if (current_layer_i == start_layer_i)
        {
            // The starting radius is generally not on a stored slice, so evaluate the EOS there.
            eos_solution_storage_ptr->call_material(current_layer_i, starting_radius, eos_material_state);

            starting_gravity = eos_material_state.gravity;
            // TODO: At very small r the interpolated g can come back negative (likely an EOS artifact);
            // clamp it to a small positive floor so the shooting start is well defined.
            if (starting_gravity < TidalPyConstants::d_EPS)
            {
                starting_gravity = TidalPyConstants::d_EPS;
            }
            starting_density = eos_material_state.density;
            starting_shear   = eos_material_state.shear_modulus;
            starting_bulk    = eos_material_state.bulk_modulus;

            radius_lower  = starting_radius;
            gravity_lower = starting_gravity;
            density_lower = starting_density;
            shear_lower   = starting_shear;
            bulk_lower    = starting_bulk;
        }
        else
        {
            // Ask the solution for this layer's base rather than reading the slice arrays. The two agree at a
            // slice radius, but going through call_material means the material-state provider is honoured, so a
            // world solve takes its interface values from the same models the integration uses.
            radius_lower = layer_radius_ptr[0];
            eos_solution_storage_ptr->call_material(current_layer_i, radius_lower, eos_material_state);
            gravity_lower = eos_material_state.gravity;
            density_lower = eos_material_state.density;
            shear_lower   = eos_material_state.shear_modulus;
            bulk_lower    = eos_material_state.bulk_modulus;
        }

        const double radius_upper = layer_radius_ptr[layer_slices - 1];
        eos_solution_storage_ptr->call_material(current_layer_i, radius_upper, eos_material_state);
        const double gravity_upper = eos_material_state.gravity;
        const double density_upper = eos_material_state.density;

        if (max_step_from_arrays)
        {
            max_step_to_use = std::abs(0.33 * (radius_upper - radius_lower));
            max_step_to_use = std::fmax(max_step_to_use, d_EPS_DBL_10000);
        }

        const int layer_type       = layer_types_ptr[current_layer_i];
        const bool layer_is_static = is_static_by_layer_ptr[current_layer_i];
        const bool layer_is_incomp = is_incompressible_by_layer_ptr[current_layer_i];

        if (scale_rtols_by_layer_type)
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
                initial_y_ptr[y_i] = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
            }
            initial_y_only_real_ptr[y_i] = TidalPyConstants::d_NAN;
        }

        if (current_layer_i == start_layer_i)
        {
            c_find_starting_conditions(
                &solution_storage_ptr->success,
                solution_storage_ptr->message,
                layer_type,
                layer_is_static,
                layer_is_incomp,
                use_kamata,
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
                break;
            }
        }
        else
        {
            layer_below_type      = layer_types_ptr[current_layer_i - 1];
            layer_below_is_static = is_static_by_layer_ptr[current_layer_i - 1];
            layer_below_is_incomp = is_incompressible_by_layer_ptr[current_layer_i - 1];

            // Interface gravity: mean of the bottom of this layer and the top of the one below.
            interface_gravity = 0.5 * (gravity_lower + last_layer_upper_gravity);

            // Liquid density at the interface: the liquid side of a solid-liquid pair, the static side of a
            // liquid-liquid pair, NaN when neither side needs it.
            if ((layer_type == 0) && (layer_below_type == 0))
            {
                static_liquid_density = TidalPyConstants::d_NAN;
            }
            else if (!(layer_type == 0) && (layer_below_type == 0))
            {
                static_liquid_density = density_lower;
            }
            else if ((layer_type == 0) && !(layer_below_type == 0))
            {
                static_liquid_density = last_layer_upper_density;
            }
            else
            {
                if (layer_is_static && layer_below_is_static)
                {
                    // TODO: Not sure what to do here so just use this layer's density.
                    static_liquid_density = density_lower;
                }
                else if (layer_is_static && !layer_below_is_static)
                {
                    static_liquid_density = density_lower;
                }
                else if (!layer_is_static && layer_below_is_static)
                {
                    static_liquid_density = last_layer_upper_density;
                }
                else
                {
                    static_liquid_density = TidalPyConstants::d_NAN;
                }
            }

            c_solve_upper_y_at_interface(
                uppermost_y_per_solution_ptr,
                initial_y_ptr,
                layer_below_num_sols,
                num_sols,
                C_MAX_NUM_Y,
                layer_below_type,
                layer_below_is_static,
                layer_below_is_incomp,
                layer_type,
                layer_is_static,
                layer_is_incomp,
                interface_gravity,
                static_liquid_density,
                G_to_use
            );
        }

        for (size_t i = 0; i < 18; ++i)
        {
            uppermost_y_per_solution_ptr[i] = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
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
                expected_size,             // Expected final integration size (0 = find good value) [size_t]
                num_extra,                 // Number of extra parameters captured during integration [size_t]
                diffeq_args_vec,           // Extra input args to diffeq vector[char]
                max_num_steps,             // Max number of steps (0 = find good value) [size_t]
                max_ram_MB,                // Max amount of RAM allowed [size_t]
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

            // Top-of-layer y for the next layer's interface condition.
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
                uppermost_y_per_solution_ptr[solution_i * C_MAX_NUM_Y + y_i] =
                    std::complex<double>(interp_top[2 * y_i], interp_top[2 * y_i + 1]);
            }

            // Retain the dense interpolant and arm a fresh result for the next integration.
            solution_storage_ptr->p_interp_by_layer_sol[current_layer_i][solution_i] =
                std::move(integration_solution_uptr);
            integration_solution_uptr = std::make_unique<CySolverResult>(integration_method);
            integration_solution_ptr  = integration_solution_uptr.get();
        }
        if (solution_storage_ptr->error_code != 0)
        {
            solution_storage_ptr->success = false;
            return solution_storage_ptr->error_code;
        }

        layer_below_num_sols     = num_sols;
        last_layer_upper_gravity = gravity_upper;
        last_layer_upper_density = density_upper;
    }



    if (solution_storage_ptr->error_code < 0 || !solution_storage_ptr->success)
    {
        solution_storage_ptr->success = false;
        if (integration_solution_ptr)
        {
            solution_storage_ptr->message =
                std::string("RadialSolver.ShootingMethod:: Integration failed:\n\t") +
                integration_solution_ptr->message + std::string("\n");
        }

        if (verbose)
        {
            printf("%s", solution_storage_ptr->message.c_str());
        }
        return solution_storage_ptr->error_code;
    }
    else
    {
        solution_storage_ptr->message = std::string("Integration completed for all layers. Beginning solution collapse.\n");

        double layer_above_lower_gravity = TidalPyConstants::d_NAN;
        double layer_above_lower_density = TidalPyConstants::d_NAN;
        int layer_above_type             = 9;
        bool layer_above_is_static       = false;
        bool layer_above_is_incomp       = false;

        // The storage is reused across solves.
        solution_storage_ptr->surface_amplification = 0.0;

        for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
        {
            solution_storage_ptr->message =
                std::string("Collapsing radial solutions for \"") +
                std::to_string(ytype_i) +
                std::string("\" solver.\n");

            bc_solution_info            = -999;
            constant_vector_ptr[0]      = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
            constant_vector_ptr[1]      = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
            constant_vector_ptr[2]      = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);

            if (verbose)
            {
                printf("%s", solution_storage_ptr->message.c_str());
            }

            // Collapse from the surface down to the starting layer.
            for (size_t current_layer_i = 0; current_layer_i < num_layers - start_layer_i; ++current_layer_i)
            {
                const size_t layer_i_reversed = num_layers - (current_layer_i + 1);

                size_t first_slice_index = 0;
                size_t layer_slices      = num_slices_by_layer_vec[layer_i_reversed];
                if (layer_i_reversed == start_layer_i)
                {
                    first_slice_index = last_index_before_start + 1;
                    layer_slices -= start_index_in_layer;
                }
                else
                {
                    first_slice_index = first_slice_index_by_layer_vec[layer_i_reversed];
                }

                const size_t num_sols = num_solutions_by_layer_ptr[layer_i_reversed];
                const size_t num_ys   = 2 * num_sols;

                double* layer_radius_ptr = &radius_array_ptr[first_slice_index];

                double radius_lower{TidalPyConstants::d_NAN};
                double gravity_lower{TidalPyConstants::d_NAN};
                double density_lower{TidalPyConstants::d_NAN};
                std::complex<double> shear_lower{TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
                std::complex<double> bulk_lower{TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};

                if (layer_i_reversed == start_layer_i)
                {
                    radius_lower  = starting_radius;
                    gravity_lower = starting_gravity;
                    density_lower = starting_density;
                    shear_lower   = starting_shear;
                    bulk_lower    = starting_bulk;
                }
                else
                {
                    // Through call_material(), as above, so the provider supplies the interface values.
                    radius_lower = layer_radius_ptr[0];
                    eos_solution_storage_ptr->call_material(
                        layer_i_reversed, radius_lower, eos_material_state);
                    gravity_lower = eos_material_state.gravity;
                    density_lower = eos_material_state.density;
                    shear_lower   = eos_material_state.shear_modulus;
                    bulk_lower    = eos_material_state.bulk_modulus;
                }

                const double radius_upper = layer_radius_ptr[layer_slices - 1];
                eos_solution_storage_ptr->call_material(layer_i_reversed, radius_upper, eos_material_state);
                const double density_upper = eos_material_state.density;
                const double gravity_upper = eos_material_state.gravity;

                const int layer_type      = layer_types_ptr[layer_i_reversed];
                const bool layer_is_static = is_static_by_layer_ptr[layer_i_reversed];
                const bool layer_is_incomp = is_incompressible_by_layer_ptr[layer_i_reversed];

                // Evaluate each independent solution's y-values at the top of the layer from its dense interpolant.
                for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
                {
                    double interp_top[C_MAX_NUM_Y_REAL];
                    c_call_dense_checked(
                        solution_storage_ptr->p_interp_by_layer_sol[layer_i_reversed][solution_i].get(),
                        radius_upper,
                        interp_top,
                        2 * num_ys);
                    for (size_t y_i = 0; y_i < num_ys; ++y_i)
                    {
                        uppermost_y_per_solution_ptr[solution_i * C_MAX_NUM_Y + y_i] =
                            std::complex<double>(interp_top[2 * y_i], interp_top[2 * y_i + 1]);
                    }
                }

                if (current_layer_i == 0)
                {
                    c_apply_surface_bc(
                        constant_vector_ptr,
                        &bc_solution_info,
                        bc_pointer,
                        uppermost_y_per_solution_ptr,
                        surface_gravity,
                        G_to_use,
                        num_sols,
                        C_MAX_NUM_Y,
                        ytype_i,
                        layer_type,
                        layer_is_static,
                        layer_is_incomp
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
                            uppermost_y_per_solution_ptr,
                            num_sols,
                            num_ys,
                            C_MAX_NUM_Y));
                }
                else
                {
                    // Interior layer: constants follow from the layer above.
                    c_top_to_bottom_interface_bc(
                        constant_vector_ptr,
                        layer_above_constant_vector_ptr,
                        uppermost_y_per_solution_ptr,
                        gravity_upper,
                        layer_above_lower_gravity,
                        density_upper,
                        layer_above_lower_density,
                        layer_type,
                        layer_above_type,
                        layer_is_static,
                        layer_above_is_static,
                        layer_is_incomp,
                        layer_above_is_incomp,
                        num_sols,
                        C_MAX_NUM_Y
                    );
                }

                // The collapsed y is evaluated on demand from the interpolants and these constants
                // (c_RadialSolutionStorage::get_radial_solution_nondim); unused entries stay NaN.
                std::array<std::complex<double>, 3>& dest_constants =
                    solution_storage_ptr->p_constants_by_ytype_layer[ytype_i][layer_i_reversed];
                for (size_t s = 0; s < num_sols; ++s)
                    dest_constants[s] = constant_vector_ptr[s];

                layer_above_lower_gravity = gravity_lower;
                layer_above_lower_density = density_lower;
                layer_above_type          = layer_type;
                layer_above_is_static     = layer_is_static;
                layer_above_is_incomp     = layer_is_incomp;

                if (num_sols == 1)
                {
                    layer_above_constant_vector_ptr[0] = constant_vector_ptr[0];
                    layer_above_constant_vector_ptr[1] =
                        std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
                    layer_above_constant_vector_ptr[2] =
                        std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
                }
                else if (num_sols == 2)
                {
                    layer_above_constant_vector_ptr[0] = constant_vector_ptr[0];
                    layer_above_constant_vector_ptr[1] = constant_vector_ptr[1];
                    layer_above_constant_vector_ptr[2] =
                        std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
                }
                else if (num_sols == 3)
                {
                    layer_above_constant_vector_ptr[0] = constant_vector_ptr[0];
                    layer_above_constant_vector_ptr[1] = constant_vector_ptr[1];
                    layer_above_constant_vector_ptr[2] = constant_vector_ptr[2];
                }
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