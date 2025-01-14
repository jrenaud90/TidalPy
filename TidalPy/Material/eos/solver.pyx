# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.stdio cimport sprintf, printf
from libc.string cimport strcpy
from libcpp.memory cimport make_shared

from CyRK cimport CySolverResult, DiffeqFuncType, cysolve_ivp_noreturn

from TidalPy.constants cimport d_G, d_PI_DBL, d_INF_DBL, d_EPS_DBL_100
from TidalPy.Material.eos.eos_solution cimport EOS_Y_VALUES, EOS_EXTRA_VALUES
from TidalPy.Material.eos.ode cimport eos_diffeq


cdef void solve_eos(
        EOSSolutionCC* eos_solution_ptr,
        vector[PreEvalFunc] eos_function_bylayer_ptr_vec,
        vector[EOS_ODEInput] eos_input_bylayer_vec,
        double planet_bulk_density,
        double surface_pressure = 0.0,
        double G_to_use = d_G,
        unsigned int integration_method = 2,
        double rtol = 1.0e-6,
        double atol = 1.0e-10,
        double pressure_tol = 1.0e-3,
        unsigned int max_iters = 100,
        cpp_bool verbose = True
        ) noexcept nogil:

    # Set the message assuming success, it will be updated if we run into failure 
    strcpy(eos_solution_ptr.message_ptr, "Equation of state solver finished without issue.")

    # Don't use rtol or atol arrays
    # cdef double[2] rtols_arr = [1.0e-8, 1.0e-12]
    # cdef double[2] atols_arr = [1.0e-10, 1.0e-15]
    cdef double* rtols_ptr = NULL #&rtols_arr[0]
    cdef double* atols_ptr = NULL #&atols_arr[0]
    
    # Determine planetary properties
    cdef size_t len_radius_array  = eos_solution_ptr.radius_array_vec.size()
    cdef double planet_radius     = eos_solution_ptr.radius_array_vec.back()
    cdef double r0_gravity        = 0.0
    cdef double r0_pressure_guess = (2. / 3.) * d_PI_DBL * G_to_use * planet_radius**2 * planet_bulk_density**2
    r0_pressure_guess += surface_pressure
    
    cdef double r0_mass = 0.0
    cdef double r0_moi  = 0.0
    
    # Setup bound variables
    cdef double[2] radial_span
    cdef double* radial_span_ptr = &radial_span[0]
    cdef DiffeqFuncType diffeq = eos_diffeq

    # We need the centeral pressure of the planet. Use the global bulk modulus to calculate this.
    cdef double[4] y0   = [r0_gravity, r0_pressure_guess, r0_mass, r0_moi]
    cdef double* y0_ptr = &y0[0]

    # y information
    cdef size_t num_y     = EOS_Y_VALUES
    cdef size_t num_extra = EOS_EXTRA_VALUES

    # Layer information
    cdef size_t layer_i
    cdef size_t top_of_last_layer_index
    cdef double[4] y0_bylayer
    cdef double* y0_bylayer_ptr = &y0_bylayer[0]

    # EOS functions and inputs
    cdef char* args_ptr                    = NULL
    cdef EOS_ODEInput* eos_input_layer_ptr = NULL
    cdef size_t sizeof_args                = sizeof(EOS_ODEInput)
    
    # Other integration information
    cdef double max_step
    cdef size_t max_num_steps  = 10_000
    cdef size_t max_ram_MB     = 500
    cdef bint use_dense_output = False
    cdef double first_step     = 0.0
    cdef size_t expected_size  = 200
    cdef double* t_eval_ptr    = NULL
    cdef size_t len_t_eval     = 0
    
    # Pressure convergence variables
    cdef size_t surface_pressure_index   = 0
    cdef double calculated_surf_pressure = d_INF_DBL
    cdef double pressure_diff            = d_INF_DBL
    cdef int iterations                  = 0
    cdef cpp_bool failed                 = False
    cdef cpp_bool max_iters_hit          = False
    cdef cpp_bool final_run              = False
    
    # Integration solution variables
    cdef size_t last_solution_size = 0
    cdef shared_ptr[CySolverResult] integration_result_sptr
    cdef CySolverResult* integration_result_ptr = NULL

    # Loop variables
    cdef size_t y_i
    cdef size_t num_layers = eos_solution_ptr.num_layers

    # Solve the equation of state in a convergence loop based on the surface pressure.
    while True:
        
        # Reset calculated pressure 
        calculated_surf_pressure = d_INF_DBL

        if not final_run:
            iterations += 1

        # Step through each macro layer of the planet and solve the equation of state starting from bottom to top
        for layer_i in range(num_layers):
            # Setup bounds and initial conditions for next layer's integration
            radial_span_ptr[1] = eos_solution_ptr.upper_radius_bylayer_vec[layer_i]
            if layer_i == 0:
                radial_span_ptr[0] = 0.0
                # Set y0 for bottom-most layer equal to the global y0
                for y_i in range(num_y):
                    y0_bylayer_ptr[y_i] = y0_ptr[y_i]

            # Set the maximum step size equal to 1/3 the layer's thickness
            max_step = 0.33 * (radial_span_ptr[1] - radial_span_ptr[0])
                
            # Get eos function and inputs for this layer
            eos_input_layer_ptr = &eos_input_bylayer_vec[layer_i]
            
            if final_run:
                # We now want to make sure that all final calculations are performed.
                eos_input_layer_ptr.update_bulk  = True
                eos_input_layer_ptr.update_shear = True
                eos_input_layer_ptr.final_solve  = True
                # Capture extra outputs and store interpolators
                num_extra = EOS_EXTRA_VALUES
                use_dense_output = True
            else:
                # During the iterations we do not need to update the complex bulk or shear
                eos_input_layer_ptr.update_bulk  = False
                eos_input_layer_ptr.update_shear = False
                # We also are not at the final call step. 
                eos_input_layer_ptr.final_solve  = False
                num_extra = 0
                use_dense_output = False
            
            # Convert input pointer to void pointer (required by cysolve)
            args_ptr = <char*>eos_input_layer_ptr

            ###### Radial Integrate the EOS Through the Planet ######
            integration_result_sptr = make_shared[CySolverResult](
                num_y,
                num_extra,
                expected_size,
                radial_span_ptr[1],
                True,
                use_dense_output,
                False)

            cysolve_ivp_noreturn(
                integration_result_sptr,
                diffeq,                        # Differential equation [DiffeqFuncType]
                radial_span_ptr,               # Radial span [const double*]
                y0_bylayer_ptr,                # y0 array [const double*]
                num_y,                         # Number of dependent y values [size_t]
                integration_method,            # Integration method [int]
                rtol,                          # Relative Tolerance (as scalar) [double]
                atol,                          # Absolute Tolerance (as scalar) [double]
                args_ptr,                      # Extra input args to diffeq [char*]
                sizeof_args,                   # Size of argument structure in bytes [size_t]
                num_extra,                     # Number of extra outputs tracked [size_t]
                max_num_steps,                 # Max number of steps (0 = find good value) [size_t]
                max_ram_MB,                    # Max amount of RAM allowed [size_t]
                use_dense_output,              # Use dense output [bint]
                t_eval_ptr,                    # Interpolate at radius array [double*]
                len_t_eval,                    # Size of interpolation array [size_t]
                eos_function_bylayer_ptr_vec[layer_i],  # Pre-eval function used in diffeq [PreEvalFunc]
                rtols_ptr,                     # Relative Tolerance (as array) [double*]
                atols_ptr,                     # Absolute Tolerance (as array) [double*]
                max_step,                      # Maximum step size [double]
                first_step,                    # Initial step size (0 = find good value) [double]
                expected_size                  # Expected final integration size (0 = find good value) [size_t]
                )

            integration_result_ptr = integration_result_sptr.get()
            #########################################################
            last_solution_size = integration_result_ptr.size
            eos_solution_ptr.save_steps_taken(integration_result_ptr.steps_taken)

            if not integration_result_ptr.success:
                failed = True
            
            if final_run and (not failed):
                eos_solution_ptr.save_cyresult(integration_result_sptr)
                # eos_solution_ptr.cysolver_results_sptr_bylayer_vec.push_back(integration_result)
                # eos_solution_ptr.current_layers_saved += 1
            elif (layer_i == num_layers - 1) and (not failed):
                # Find planet surface pressure for this iteration
                surface_pressure_index   = (last_solution_size * num_y) - (num_y - 2) - 1    # (Total number of slices) - (num_y - location of pressure) - 1
                calculated_surf_pressure = integration_result_ptr.solution[surface_pressure_index]
            
            # Prepare for next layer
            if (num_layers > 1) and (not failed):
                # Bottom radius value equals top of lower layer's radius
                radial_span_ptr[0] = eos_solution_ptr.upper_radius_bylayer_vec[layer_i]
                top_of_last_layer_index = (num_extra + num_y) * (last_solution_size - 1)

                # y0 for this layer equals the top most result of the lower layer
                if integration_result_ptr:
                    for y_i in range(num_y):
                        y0_bylayer_ptr[y_i] = integration_result_ptr.solution[top_of_last_layer_index + y_i]
                else:
                    # Not sure why that would be null but in any case we are in a fail state.
                    failed = True

            # Clear pointers
            integration_result_ptr = NULL
            integration_result_sptr.reset()

            if failed:
                break

        if failed:
            break
        
        if final_run:
            # We are done!
            break
        else:
            # Update the central pressure using the error at the surface as the correction factor
            pressure_diff     = surface_pressure - calculated_surf_pressure
            pressure_diff_abs = pressure_diff
            if pressure_diff < 0.0:
                pressure_diff_abs = -pressure_diff

            # Calculate percent difference to use in convergence check.
            if surface_pressure > d_EPS_DBL_100:
                pressure_diff_abs /= surface_pressure
            
            # Check if we are done next iteration
            if pressure_diff_abs <= pressure_tol:
                final_run = True
            else:
                y0_ptr[1] += pressure_diff
        
        if iterations >= max_iters:
            max_iters_hit = True
            eos_solution_ptr.max_iters_hit = True
            # To ensure that there is some output we will go ahead and do a final run.
            final_run = True
    
    # Done with convergence loop.
    eos_solution_ptr.iterations = iterations
    
    # Display any warnings
    if max_iters_hit:
        strcpy(eos_solution_ptr.message_ptr, "Warning in `solve_eos`: Maximum number of iterations hit without convergence.")
        if verbose:
            printf(eos_solution_ptr.message_ptr)

    if failed:
        eos_solution_ptr.success = False
        if integration_result_ptr:
            sprintf(eos_solution_ptr.message_ptr, "Warning in `solve_eos`: Integrator failed at iteration %d. Message: %s", iterations, integration_result_ptr.message_ptr)
        else:
            sprintf(eos_solution_ptr.message_ptr, "Warning in `solve_eos`: Integrator failed at iteration %d.", iterations)
        if verbose:
            printf(eos_solution_ptr.message_ptr)
    else:
        # Set feedback attributes
        eos_solution_ptr.success = True

        # Set other final parameters
        eos_solution_ptr.pressure_error = pressure_diff_abs

        # Tell the eos solution to perform a full planet interpolation and store the results. Including surface results 
        eos_solution_ptr.interpolate_full_planet()
    
    if integration_result_sptr:
        integration_result_sptr.reset()
    
    if integration_result_ptr:
        integration_result_ptr = NULL