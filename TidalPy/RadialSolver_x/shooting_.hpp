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
    /** Determine number of solutions required for layer based on assumptions.
    
    Parameters
    ----------
    layer_type : int
        - 0: Layer is solid 
        - 1: Layer is liquid
    is_static : bool
        Use static (True) or dynamic (False) assumption.
    is_incompressible : bool
        Use incompressible (True) or compressible (False) assumption.

    Returns
    -------
    num_sols : int
        Number of solutions required for layer.

    */

    // Initialize
    size_t num_sols = 0;

    if (layer_type == 0)
    {
        // Solid
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
            // Dynamic
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
        // Liquid
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
            // Dynamic
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
    bool verbose
) noexcept
{
    /** Solves the viscoelastic-gravitational problem for planets using a shooting method.
    */

    // Get raw pointer of radial solver storage and eos storage
    c_EOSSolution* eos_solution_storage_ptr = solution_storage_ptr->get_eos_solution_ptr();

    solution_storage_ptr->message = std::string("RadialSolver.ShootingMethod:: Starting integration\n");
    if (verbose)
    {
        printf("%s", solution_storage_ptr->message.c_str());
    }

    // Type conversions
    const double degree_l_dbl = static_cast<double>(degree_l);

    // Alias pointers to EOS properties
    double* radius_array_ptr  = eos_solution_storage_ptr->radius_array_vec.data();
    double* gravity_array_ptr = eos_solution_storage_ptr->gravity_array_vec.data();
    // double* pressure_array_ptr = eos_solution_storage_ptr->pressure_array_vec;   // Unused
    double* density_array_ptr = eos_solution_storage_ptr->density_array_vec.data();
    std::complex<double>* complex_shear_array_ptr = eos_solution_storage_ptr->complex_shear_array_vec.data();
    std::complex<double>* complex_bulk_array_ptr  = eos_solution_storage_ptr->complex_bulk_array_vec.data();

    // Pull out key information
    const size_t num_layers   = eos_solution_storage_ptr->num_layers;
    const size_t total_slices = eos_solution_storage_ptr->radius_array_size;
    const size_t top_slice_i  = total_slices - 1;
    
    // Pull out any constants now that arrays have had dimensional protocol applied to them.
    const double planet_radius   = eos_solution_storage_ptr->radius;
    const double surface_gravity = eos_solution_storage_ptr->surface_gravity;

    // Find boundary condition at the top of the planet -- this is dependent on the forcing type.
    //     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    // The [5] represents the maximum number of solvers that can be invoked with a single call to radial_solver
    const size_t num_ytypes = num_bc_models;

    // 15 = 5 (max_num_solutions) * 3 (number of surface conditions)
    double boundary_conditions[15];
    double* bc_pointer = &boundary_conditions[0];
    c_get_surface_bc(
        bc_pointer,  // Changed parameter
        bc_models_ptr,
        num_ytypes,
        planet_radius,
        planet_bulk_density,
        degree_l_dbl
    );

    // Integration information
    // Number of extra parameters captured during integration
    const size_t num_extra            = 0;
    const double first_step_size      = 0.0;
    const bool   capture_dense_output = false;
    // Events (empty — not used by RadialSolver)
    std::vector<Event> events_vec;

    // Max step size
    double max_step_to_use      = TidalPyConstants::d_NAN;
    bool   max_step_from_arrays = false;
    
    if (max_step == 0.0)
    {
        // If max_step is zero use the array information to determine max_step_size
        max_step_from_arrays = true;
    }
    else
    {
        // Otherwise use user input.
        max_step_to_use = max_step;
    }

    // Setup tolerance arrays
    // For simplicity just make these all as large as the maximum number of ys.
    // Maximum number of ys = 6. Then 2x for conversion from complex to real
    std::vector<double> rtols_vec{integration_rtol};
    std::vector<double> atols_vec{integration_atol};

    // Create storage for flags and information about each layer.
    std::vector<size_t> num_solutions_by_layer_vec;
    num_solutions_by_layer_vec.resize(num_layers);
    size_t* num_solutions_by_layer_ptr = num_solutions_by_layer_vec.data();

    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        // Pull out information on this layer
        const int layer_type            = layer_types_ptr[current_layer_i];
        const bool layer_is_static      = is_static_by_layer_ptr[current_layer_i];
        const bool layer_is_incomp      = is_incompressible_by_layer_ptr[current_layer_i];
        const double layer_upper_radius = eos_solution_storage_ptr->upper_radius_bylayer_vec[current_layer_i];

        // Find number of solutions based on this layer's assumptions
        const size_t num_sols = c_find_num_shooting_solutions(
            layer_type,
            layer_is_static,
            layer_is_incomp
        );
        num_solutions_by_layer_ptr[current_layer_i] = num_sols;
    }

    // We have all the size information needed to build storage pointers
    // Main storage pointer is setup like [current_layer_i][solution_i][y_i + r_i]
    std::vector<std::vector<std::vector<std::complex<double>>>> main_storage_vec;

    try
    {
        main_storage_vec.resize(num_layers);

        for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
        {
            const size_t num_sols     = num_solutions_by_layer_ptr[current_layer_i];
            const size_t layer_slices = num_slices_by_layer_vec[current_layer_i];
            
            // Number of ys = 2x num sols
            const size_t num_ys = 2 * num_sols;

            // Resize the middle vector (storage by solution)
            main_storage_vec[current_layer_i].resize(num_sols);

            for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
            {
                // Resize the inner vector (storage by y). 
                // Note: std::vector automatically zero-initializes std::complex elements
                main_storage_vec[current_layer_i][solution_i].resize(layer_slices * num_ys);
            }
        }
    }
    catch (const std::bad_alloc& e)
    {
        throw std::runtime_error("Failed to allocate memory for: main_storage_vec (radial_solver; init)");
    }

    // Create storage for uppermost ys for each solution. We don't know how many solutions or ys per layer so assume the
    //  worst.
    std::complex<double> uppermost_y_per_solution[18];
    std::complex<double>* uppermost_y_per_solution_ptr = &uppermost_y_per_solution[0];

    for (size_t i = 0; i < 18; ++i)
    {
        uppermost_y_per_solution_ptr[i] = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
    }

    // Layer specific pointers; set the size based on the layer with the most slices.
    std::vector<double> layer_radius_vec(eos_solution_storage_ptr->radius_array_vec.size());

    // Starting solutions (initial conditions / lower boundary conditions)
    // Allocate memory for the initial value arrays now. We don't know the number of solutions or ys. But the max number
    //  is not all that different from the min. So it is more efficient to assume the largest size.
    // Largest size = 6 (ys) x 3 (sols) = 18
    // For the "only_real" this is further multiplied by 2 since we convert the number of _complex_ ys to 2x _real_ ys
    std::complex<double> initial_y[18];
    std::complex<double>* initial_y_ptr = &initial_y[0];
    double initial_y_only_real[36];
    double* initial_y_only_real_ptr = &initial_y_only_real[0];
    bool starting_y_check = false;

    // No matter the results of the integration, we know the shape and size of the final solution.
    // The number of rows will depend on if the user wants to simultaneously calculate loading Love numbers.
    const size_t num_output_ys = C_MAX_NUM_Y * num_ytypes;

    // Get a reference pointer to solution array
    // Cast the solution pointer from double to double complex
    std::complex<double>* solution_ptr = reinterpret_cast<std::complex<double>*>(solution_storage_ptr->full_solution_vec.data());

    // Layer's differential equation will vary by layer type
    double eos_interp_array[9];  // The "9" here is the number of variables (dependent + extra) found in the EOS solver.
    // See TidalPy.Material.eos.solver.pyx for details.
    double* eos_interp_array_ptr = &eos_interp_array[0];
    std::unique_ptr<CySolverResult> integration_solution_uptr = std::make_unique<CySolverResult>(integration_method);
    CySolverResult* integration_solution_ptr = nullptr;

    // The radial solver diffeq's require additional inputs other than "y" and "r". Build a structure that stores these
    // extra arguments will be passed to the diffeq's while solving the ODE.
    std::vector<char> diffeq_args_vec(sizeof(c_RadialSolverArgs));
    c_RadialSolverArgs* diffeq_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(diffeq_args_vec.data());
    PreEvalFunc diffeq_preeval_ptr = nullptr;
    std::vector<double> y0_vec(C_MAX_NUM_Y_REAL);

    // The diffeq needs to be able to call the EOS solution which is different for each layer.
    // We will load the EOS solution into this argument structure for each layer. For now, set it to null.

    // Set diffeq inputs that do not change with layer
    diffeq_args_ptr->degree_l         = degree_l_dbl;
    diffeq_args_ptr->lp1              = degree_l_dbl + 1.0;
    diffeq_args_ptr->lm1              = degree_l_dbl - 1.0;
    diffeq_args_ptr->llp1             = degree_l_dbl * (degree_l_dbl + 1.0);
    diffeq_args_ptr->G                = G_to_use;
    diffeq_args_ptr->grav_coeff       = 4.0 * TidalPyConstants::d_PI * G_to_use;
    diffeq_args_ptr->frequency        = frequency;
    diffeq_args_ptr->layer_index      = 0;
    diffeq_args_ptr->eos_solution_ptr = eos_solution_storage_ptr;

    // The constant vectors are the same size as the number of solutions in the layer. But since the largest they can
    //  ever be is 3, it is more efficient to just preallocate them on the stack rather than dynamically allocate them
    //  on the heap.
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

    // Variables used to solve the linear equation at the planet's surface.
    // Info = flag set by the solver. Set equal to -999. This will indicate that the solver has not been called yet.
    int bc_solution_info = -999;

    // The radial solver can not start at r=0 (singularity) and the higher in the planet the better the stability of
    // the solution. The flip side is that the more of the planet you skip the less accurate your results will be. 
    // Below we determine a reasonable starting radius for the radial solver solution. 
    if (starting_radius == 0.0)
    {
        // User did not provide a starting radius.
        // Use a model involving planet radius and degree l to determine a reasonable choice. 
        // This formulism is based on Hilary Martens thesis and LoadDef manual.
        starting_radius = planet_radius * std::pow(start_radius_tolerance, 1.0 / degree_l_dbl);

        // Ensure the starting radius is not too close to the surface of the planet.
        starting_radius = std::fmin(starting_radius, 0.95 * planet_radius);
    }

    // Determine which layer this starting radius resides in. We will skip the lower layers
    size_t start_layer_i           = 0;
    size_t last_index_before_start = 0;
    size_t start_index_in_layer    = 0;
    double last_radius_check       = 0.0;
    double layer_upper_radius      = TidalPyConstants::d_INF;
    double last_layer_upper_radius = TidalPyConstants::d_INF;
    for (size_t current_layer_i = 0; current_layer_i < num_layers; ++current_layer_i)
    {
        // Check if the radius is in this layer
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
            // It is!
            start_layer_i = current_layer_i;
            const size_t first_slice_index = first_slice_index_by_layer_vec[current_layer_i];
            
            // Now find the last radial slice before the starting radius
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
                        // We found that the starting radius is in-between this slice and the last slice.
                        // We want to set the last index to the slice before this one.
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

    // Step through the solution vector and NAN out data below the starting radius
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
    // Main Diffeq Loop!
    // Steps through each layer starting from the bottom; finds starting conditions; solves the ode and stores results.
    // =================================================================================================================
    // Initialize parameters for layer below
    size_t layer_below_num_sols  = 0;
    int layer_below_type         = -1;
    bool layer_below_is_static   = false;
    bool layer_below_is_incomp   = false;
    double interface_gravity        = TidalPyConstants::d_NAN;
    double static_liquid_density    = TidalPyConstants::d_NAN;
    double last_layer_upper_gravity = TidalPyConstants::d_NAN;
    double last_layer_upper_density = TidalPyConstants::d_NAN;
    
    // Physical parameters at the very bottom of the planet.
    double starting_gravity;
    double starting_density;
    std::complex<double> starting_shear{TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    std::complex<double> starting_bulk{TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};

    for (size_t current_layer_i = start_layer_i; current_layer_i < num_layers; ++current_layer_i)
    {
        // Get layer's index information
        size_t layer_slices = num_slices_by_layer_vec[current_layer_i];
        size_t first_slice_index;
        if (current_layer_i == start_layer_i)
        {
            // The first slice index is not going to actually be the bottom-most index
            // Instead it is the one right at or above our starting radius.
            first_slice_index = last_index_before_start + 1;

            // When we loop through slices we only want to loop between the starting slice and the top of the layer
            layer_slices -= start_index_in_layer;
        }
        else
        {
            first_slice_index = first_slice_index_by_layer_vec[current_layer_i];
        }

        // Get solution and y information
        const size_t num_sols   = num_solutions_by_layer_ptr[current_layer_i];
        const size_t num_ys     = 2 * num_sols;
        const size_t num_ys_dbl = 2 * num_ys;

        // Setup pointer array slices starting at the start of this layer (either base or at starting index)
        double* layer_radius_ptr    = &radius_array_ptr[first_slice_index];
        double* layer_density_ptr   = &density_array_ptr[first_slice_index];
        double* layer_gravity_ptr   = &gravity_array_ptr[first_slice_index];
        std::complex<double>* layer_shear_mod_ptr = &complex_shear_array_ptr[first_slice_index];
        std::complex<double>* layer_bulk_mod_ptr  = &complex_bulk_array_ptr[first_slice_index];

        // Populate our t_eval vector with this layer's radius values.
        layer_radius_vec.resize(layer_slices);
        std::memcpy(layer_radius_vec.data(), layer_radius_ptr, sizeof(double) * layer_slices);
        
        // Get physical parameters at the bottom of the layer
        double radius_lower{starting_radius};
        double gravity_lower{TidalPyConstants::d_NAN};
        double density_lower{TidalPyConstants::d_NAN};
        std::complex<double> shear_lower{starting_shear};
        std::complex<double> bulk_lower{starting_bulk};
        if (current_layer_i == start_layer_i)
        {
            // In the first layer we can not use the physical properties at the bottom of the arrays
            // because the starting radius may not be at the bottom.
            // Even worse, it may not be at any of the slice indices that are stored.
            // To get the most accurate result we need to perform an interpolation to find various properties at
            // this starting radius. 
            eos_solution_storage_ptr->call(current_layer_i, starting_radius, eos_interp_array_ptr);

            // Save the values, look at the "TidalPy.Material.eos.eos_solution_.hpp" to see how these are saved. 
            // We are storing these in function-global variables because they will be used again during collapse
            starting_gravity = eos_interp_array_ptr[0];
            // TODO: Sometimes at very small r the g can be negative. Probably an issue with the EOS but for now just put a force check in?
            if (starting_gravity < TidalPyConstants::d_EPS)
            {
                starting_gravity = TidalPyConstants::d_EPS;
            }
            starting_density = eos_interp_array_ptr[4];
            starting_shear   = std::complex<double>(eos_interp_array_ptr[5], eos_interp_array_ptr[6]);
            starting_bulk    = std::complex<double>(eos_interp_array_ptr[7], eos_interp_array_ptr[8]);
            
            // Now set local variables used in this loop
            radius_lower  = starting_radius;
            gravity_lower = starting_gravity;
            density_lower = starting_density;
            shear_lower   = starting_shear;
            bulk_lower    = starting_bulk;
        }
        else
        {
            // Otherwise we can just use the values at the base of the layer
            radius_lower  = layer_radius_ptr[0];
            gravity_lower = layer_gravity_ptr[0];
            density_lower = layer_density_ptr[0];
            shear_lower   = layer_shear_mod_ptr[0];
            bulk_lower    = layer_bulk_mod_ptr[0];
        }

        const double radius_upper  = layer_radius_ptr[layer_slices - 1];
        const double gravity_upper = layer_gravity_ptr[layer_slices - 1];
        const double density_upper = layer_density_ptr[layer_slices - 1];
        const std::complex<double> shear_upper   = layer_shear_mod_ptr[layer_slices - 1];
        const std::complex<double> bulk_upper    = layer_bulk_mod_ptr[layer_slices - 1];

        // Determine max step size (if not provided by user)
        if (max_step_from_arrays)
        {
            // Maximum step size during integration can not exceed ~1/3 of the layer size.
            max_step_to_use = std::abs(0.33 * (radius_upper - radius_lower));
            max_step_to_use = std::fmax(max_step_to_use, d_EPS_DBL_10000);
        }

        // Get assumptions for layer
        const int layer_type       = layer_types_ptr[current_layer_i];
        const bool layer_is_static = is_static_by_layer_ptr[current_layer_i];
        const bool layer_is_incomp = is_incompressible_by_layer_ptr[current_layer_i];

        // Determine rtols and atols for this layer.
        // Scale rtols by layer type
        if (scale_rtols_by_layer_type)
        {
            // The x2 accounts for size of double complex being stored in double vector
            rtols_vec.resize(num_ys * 2);
            atols_vec.resize(num_ys * 2);
            
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
            {
                // Default is that each layer's rtol and atol equal user input.
                // TODO: Change up the tolerance scaling between real and imaginary?
                double layer_rtol_real = integration_rtol;
                double layer_rtol_imag = integration_rtol;
                double layer_atol_real = integration_atol;
                double layer_atol_imag = integration_atol;

                // Certain layer assumptions can affect solution stability so use additional scales on the relevant rtols
                // TODO test that these scales are the best. See Issue #44
                if (layer_type == 0)
                {
                    // Solid layer
                    // Scale y2 and y3 by 0.1
                    if ((y_i == 1) || (y_i == 2))
                    {
                        // Scale both the real and imaginary portions by the same amount.
                        layer_rtol_real *= 0.1;
                        layer_rtol_imag *= 0.1;
                    }
                }
                else
                {
                    // Liquid layer
                    if (!layer_is_static)
                    {
                        // Scale dynamic liquid layer's y2 by additional 0.01.
                        if (y_i == 1)
                        {
                            // Scale both the real and imaginary portions by the same amount.
                            layer_rtol_real *= 0.01;
                            layer_rtol_imag *= 0.01;
                        }
                    }
                }
                // Populate rtol and atol pointers.
                rtols_vec[2 * y_i]     = layer_rtol_real;
                rtols_vec[2 * y_i + 1] = layer_rtol_imag;
                atols_vec[2 * y_i]     = layer_atol_real;
                atols_vec[2 * y_i + 1] = layer_atol_imag;
            }
        }

        // Set initial y to nan for now.
        // OPT: this is for debugging purposes. could likely be commented out in the future for small performance gain.
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
            // In the first layer. Use initial condition function to find initial conditions
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
                initial_y_ptr,  // Modified Variable
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
            // Not in the first layer. Use interface function to find initial conditions
            layer_below_type      = layer_types_ptr[current_layer_i - 1];
            layer_below_is_static = is_static_by_layer_ptr[current_layer_i - 1];
            layer_below_is_incomp = is_incompressible_by_layer_ptr[current_layer_i - 1];

            // Find gravity at the base interface using bottom of this layer and top of previous.
            interface_gravity = 0.5 * (gravity_lower + last_layer_upper_gravity);

            // Find the density needed for some initial conditions.
            if ((layer_type == 0) && (layer_below_type == 0))
            {
                // Both layers are solid. A liquid interface density is not needed.
                static_liquid_density = TidalPyConstants::d_NAN;
            }
            else if (!(layer_type == 0) && (layer_below_type == 0))
            {
                // Layer below is solid, this layer is liquid. Use its density.
                static_liquid_density = density_lower;
            }
            else if ((layer_type == 0) && !(layer_below_type == 0))
            {
                // Layer below is liquid, this layer is solid. Use layer below's density.
                static_liquid_density = last_layer_upper_density;
            }
            else
            {
                // Both layers are liquid. Choose the layer's density which is static.
                if (layer_is_static && layer_below_is_static)
                {
                    // Both layers are static.
                    // TODO: Not sure what to do here so just use this layer's density.
                    static_liquid_density = density_lower;
                }
                else if (layer_is_static && !layer_below_is_static)
                {
                    // This layer is static, one below is not. Use this layer's density
                    static_liquid_density = density_lower;
                }
                else if (!layer_is_static && layer_below_is_static)
                {
                    // This layer is dynamic, layer below is static. Use layer below's density
                    static_liquid_density = last_layer_upper_density;
                }
                else
                {
                    // Both layers are dynamic. Static liquid density is not needed.
                    static_liquid_density = TidalPyConstants::d_NAN;
                }
            }

            // Find the starting values for this layer using the results a the top of the previous layer + an interface
            //  function.
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

        // Reset the uppermost y value array
        for (size_t i = 0; i < 18; ++i)
        {
            uppermost_y_per_solution_ptr[i] = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
        }

        // Change initial conditions into 2x real values instead of complex for integration
        for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
        {
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
            {
                const std::complex<double> dcomplex_tmp = initial_y_ptr[solution_i * C_MAX_NUM_Y + y_i];
                initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL + 2 * y_i]     = dcomplex_tmp.real();
                initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL + 2 * y_i + 1] = dcomplex_tmp.imag();
            }
        }

        // Find correct diffeq
        DiffeqFuncType layer_diffeq = c_find_layer_diffeq(layer_type, layer_is_static, layer_is_incomp);

        // Set the layer index to current layer in the additional diffeq arg struct
        diffeq_args_ptr->layer_index = current_layer_i;

        // Solve for each solution
        for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
        {
            // Ensure the integration solution object is initialized.
            if (!integration_solution_uptr)
            {
                integration_solution_uptr = std::make_unique<CySolverResult>(integration_method);
            }
            integration_solution_ptr = integration_solution_uptr.get();

            // Copy over initial y values into our vector
            y0_vec.resize(num_ys_dbl);
            std::memcpy(y0_vec.data(), &initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL], sizeof(double) * num_ys_dbl);

            // ###### Integrate! #######
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
                layer_radius_vec,          // Interpolate at this layer's radius array vector[double]
                diffeq_preeval_ptr,        // Pre-eval function used in diffeq [PreEvalFunc]
                events_vec,                // Events vector[Event]
                rtols_vec,                 // Relative Tolerance (as array) vector[double]
                atols_vec,                 // Absolute Tolerance (as array) vector[double]
                max_step_to_use,           // Maximum step size [double]
                first_step_size,           // Initial step size (0 = find good value) [double]
                true                       // Force retain solver
            );
            // #########################

            // Store diagnostic data
            solution_storage_ptr->shooting_method_steps_taken_vec[(3 * current_layer_i) + solution_i] = 
                integration_solution_ptr->steps_taken;

            // Check for problems
            if (!integration_solution_ptr->success)
            {
                // Problem with integration.
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

            // If no problems, store results.
            // TODO: Eventually we will want to store the full integration_solution_uptr so we can use dense outputs.
            //    When we do that we will need to make sure we std::move(integration_solution_uptr) the unique pointer.
            double* integrator_data_ptr = integration_solution_ptr->solution.data();
            
            // Need to make a copy because the solver pointers will be reallocated during the next solution.
            // Get storage pointer for this solution
            std::vector<std::complex<double>>& storage_by_y_vec = main_storage_vec[current_layer_i][solution_i];

            for (size_t slice_i = 0; slice_i < layer_slices; ++slice_i)
            {
                for (size_t y_i = 0; y_i < num_ys; ++y_i)
                {
                    // Convert 2x real ys to 1x complex ys
                    const std::complex<double> dcomplex_tmp = std::complex<double>(
                        integrator_data_ptr[num_ys_dbl * slice_i + (2 * y_i)],
                        integrator_data_ptr[num_ys_dbl * slice_i + (2 * y_i) + 1]
                    );
                    storage_by_y_vec[num_ys * slice_i + y_i] = dcomplex_tmp;

                    if (slice_i == (layer_slices - 1))
                    {
                        // Store top most result for initial condition for the next layer
                        uppermost_y_per_solution_ptr[solution_i * C_MAX_NUM_Y + y_i] = dcomplex_tmp;
                    }
                }
            }
        }
        if (solution_storage_ptr->error_code != 0)
        {
            // Error was encountered during integration
            solution_storage_ptr->success = false;
            return solution_storage_ptr->error_code;
        }

        // Prepare for next layer
        layer_below_num_sols     = num_sols;
        last_layer_upper_gravity = gravity_upper;
        last_layer_upper_density = density_upper;
    } // Layer Loop (Main integration loop)



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
        // No errors. Proceed with collapsing all sub-solutions into final full solution.
        solution_storage_ptr->message = std::string("Integration completed for all layers. Beginning solution collapse.\n");

        // Tracker variables
        double layer_above_lower_gravity   = TidalPyConstants::d_NAN;
        double layer_above_lower_density   = TidalPyConstants::d_NAN;
        double liquid_density_at_interface = TidalPyConstants::d_NAN;
        int layer_above_type               = 9;
        bool layer_above_is_static         = false;
        bool layer_above_is_incomp         = false;

        for (size_t ytype_i = 0; ytype_i < num_ytypes; ++ytype_i)
        {
            solution_storage_ptr->message =
                std::string("Collapsing radial solutions for \"") +
                std::to_string(ytype_i) +
                std::string("\" solver.\n");

            // Reset variables for this solver
            bc_solution_info            = -999;
            constant_vector_ptr[0]      = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
            constant_vector_ptr[1]      = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);
            constant_vector_ptr[2]      = std::complex<double>(TidalPyConstants::d_NAN, TidalPyConstants::d_NAN);

            if (verbose)
            {
                printf("%s", solution_storage_ptr->message.c_str());
            }
            // Collapse the multiple solutions for each layer into one final combined solution.

            // Work from the surface down to the starting layer.
            for (size_t current_layer_i = 0; current_layer_i < num_layers - start_layer_i; ++current_layer_i)
            {
                const size_t layer_i_reversed = num_layers - (current_layer_i + 1);

                // Pull out layer information.
                size_t first_slice_index = 0;
                size_t layer_slices      = num_slices_by_layer_vec[layer_i_reversed];
                if (layer_i_reversed == start_layer_i)
                {
                    // The first slice index is not going to actually be the bottom-most index
                    // Instead it is the one right at or above our starting radius.
                    first_slice_index = last_index_before_start + 1;

                    // When we loop through slices we only want to loop between the starting slice and the top of the layer
                    layer_slices -= start_index_in_layer;
                }
                else
                {
                    first_slice_index = first_slice_index_by_layer_vec[layer_i_reversed];
                }

                // Get solution and y information
                const size_t num_sols = num_solutions_by_layer_ptr[layer_i_reversed];
                const size_t num_ys   = 2 * num_sols;

                // Setup pointer array slices starting at this layer's beginning
                double* layer_radius_ptr    = &radius_array_ptr[first_slice_index];
                double* layer_density_ptr   = &density_array_ptr[first_slice_index];
                double* layer_gravity_ptr   = &gravity_array_ptr[first_slice_index];
                std::complex<double>* layer_bulk_mod_ptr  = &complex_bulk_array_ptr[first_slice_index];
                std::complex<double>* layer_shear_mod_ptr = &complex_shear_array_ptr[first_slice_index];

                // Get physical parameters at the top and bottom of the layer
                double radius_lower{TidalPyConstants::d_NAN};
                double gravity_lower{TidalPyConstants::d_NAN};
                double density_lower{TidalPyConstants::d_NAN};
                std::complex<double> shear_lower{TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
                std::complex<double> bulk_lower{TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};

                if (layer_i_reversed == start_layer_i)
                {
                    // In the starting layer.
                    // Found earlier via interpolation at the starting radius
                    radius_lower  = starting_radius;
                    gravity_lower = starting_gravity;
                    density_lower = starting_density;
                    shear_lower   = starting_shear;
                    bulk_lower    = starting_bulk;
                }
                else
                {
                    radius_lower  = layer_radius_ptr[0];
                    density_lower = layer_density_ptr[0];
                    gravity_lower = layer_gravity_ptr[0];
                    bulk_lower    = layer_bulk_mod_ptr[0];
                    shear_lower   = layer_shear_mod_ptr[0];
                }

                const double radius_upper  = layer_radius_ptr[layer_slices - 1];
                const double density_upper = layer_density_ptr[layer_slices - 1];
                const double gravity_upper = layer_gravity_ptr[layer_slices - 1];

                // Get assumptions for layer
                const int layer_type      = layer_types_ptr[layer_i_reversed];
                const bool layer_is_static = is_static_by_layer_ptr[layer_i_reversed];
                const bool layer_is_incomp = is_incompressible_by_layer_ptr[layer_i_reversed];

                // Get radial solution values at the top of the layer
                for (size_t solution_i = 0; solution_i < num_sols; ++solution_i)
                {
                    std::vector<std::complex<double>>& storage_by_y_vec = 
                        main_storage_vec[layer_i_reversed][solution_i];
                    for (size_t y_i = 0; y_i < num_ys; ++y_i)
                    {
                        uppermost_y_per_solution_ptr[solution_i * C_MAX_NUM_Y + y_i] = 
                            storage_by_y_vec[(layer_slices - 1) * num_ys + y_i];
                    }
                }

                if (current_layer_i == 0)
                {
                    // Working on surface (uppermost) layer -- Apply surface boundary conditions.
                    c_apply_surface_bc(
                        constant_vector_ptr,  // Modified Variable
                        &bc_solution_info,    // Modified Variable
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

                    // Check that the boundary condition was successfully applied.
                    if (bc_solution_info != 0)
                    {
                        solution_storage_ptr->error_code = -12;
                        solution_storage_ptr->success    = false;
                        solution_storage_ptr->message    = 
                            std::string(
                                "RadialSolver.ShootingMethod:: Error encountered while applying surface"
                                "boundary condition. ZGESV code: ") + 
                                std::to_string(bc_solution_info) + 
                                std::string("\nThe solutions may not be valid at the surface.\n");
                        
                        if (verbose)
                        {
                            printf("%s", solution_storage_ptr->message.c_str());
                        }
                        return solution_storage_ptr->error_code;
                    }
                }
                else
                {
                    // Working on interior layers. Will need to find the constants of integration based on the layer above.
                    c_top_to_bottom_interface_bc(
                        constant_vector_ptr,  // Modified Variable
                        layer_above_constant_vector_ptr,
                        uppermost_y_per_solution_ptr,
                        gravity_upper, layer_above_lower_gravity,
                        density_upper, layer_above_lower_density,
                        layer_type, layer_above_type,
                        layer_is_static, layer_above_is_static,
                        layer_is_incomp, layer_above_is_incomp,
                        num_sols, C_MAX_NUM_Y
                    );
                }

                // Use constant vectors to find the full y from all of the solutions in this layer
                c_collapse_layer_solution(
                    solution_ptr,  // Modified Variable
                    constant_vector_ptr,
                    main_storage_vec[layer_i_reversed],
                    layer_radius_ptr,
                    layer_density_ptr,
                    layer_gravity_ptr,
                    frequency,
                    first_slice_index,
                    layer_slices,
                    num_sols,
                    C_MAX_NUM_Y,
                    num_ys,
                    num_output_ys,
                    ytype_i,
                    layer_type,
                    layer_is_static,
                    layer_is_incomp
                );

                // Setup for next layer
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

            // Ready for next y-type
        }
    }

    // Update solution status and return
    if (solution_storage_ptr->error_code != 0)
    {
        solution_storage_ptr->success = false;
    }
    else
    {
        solution_storage_ptr->success = true;
        solution_storage_ptr->message = std::string("RadialSolver.ShootingMethod: Completed without any noted issues.");
    }

    // Done!
    return solution_storage_ptr->error_code;
}