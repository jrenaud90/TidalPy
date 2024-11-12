# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.stdio cimport printf, sprintf
from libc.string cimport memcpy, strcpy

from CyRK cimport cysolve_ivp, CySolveOutput

from TidalPy.constants cimport d_G, d_PI_DBL, d_INF_DBL
from TidalPy.Material.eos.ode cimport eos_diffeq

cdef EOSSolutionVec solve_eos(
        cpp_bool* success_ptr,
        char* message_ptr,
        double* radius_array_ptr,
        size_t len_radius_array,
        double* layer_upper_radii,
        unsigned int num_layers,
        PreEvalFunc* eos_function_bylayer_ptrs,
        EOS_ODEInput** eos_input_bylayer_ptrs,
        double planet_bulk_density,
        double surface_pressure = 0.0,
        double G_to_use = d_G,
        unsigned int integration_method = 2,
        double rtol = 1.0e-8,
        double atol = 1.0e-16,
        double pressure_tol = 1.0e-3,
        unsigned int max_iters = 100,
        cpp_bool verbose = True
        ) noexcept nogil:

    # Feedback
    cdef char[256] local_message
    cdef char* local_message_ptr = &local_message[0]

    # Don't use rtol or atol arrays
    # cdef double[2] rtols_arr = [1.0e-8, 1.0e-12]
    # cdef double[2] atols_arr = [1.0e-10, 1.0e-15]
    cdef double* rtols_ptr = NULL #&rtols_arr[0]
    cdef double* atols_ptr = NULL #&atols_arr[0]
    
    # Determine planetary properties
    cdef double planet_radius = radius_array_ptr[len_radius_array - 1]
    cdef double r0_gravity = 0.0
    cdef double r0_pressure_guess = \
        surface_pressure + \
        (2. / 3.) * d_PI_DBL * G_to_use * planet_radius * planet_radius * planet_bulk_density * planet_bulk_density
    
    # Setup bound variables
    cdef double[2] radial_span
    cdef double* radial_span_ptr = &radial_span[0]

    # We need the centeral pressure of the planet. Use the global bulk modulus to calculate this.
    cdef double[2] y0 = [r0_gravity, r0_pressure_guess]
    cdef double* y0_ptr = &y0[0]

    # y information
    cdef unsigned int num_y = 2
    cdef unsigned int num_extra = 0

    # Layer information
    cdef unsigned int layer_i
    cdef size_t top_of_last_layer_index
    cdef double[2] y0_layer
    cdef double* y0_layer_ptr = &y0_layer[0]

    # EOS functions and inputs
    cdef void* args_ptr = NULL
    cdef PreEvalFunc eos_function_ptr = NULL
    cdef EOS_ODEInput* eos_input_layer_ptr = NULL
    
    # Other integration information
    cdef size_t max_num_steps  = 10_000
    cdef size_t max_ram_MB     = 500
    cdef bint use_dense_output = False
    cdef double max_step       = 0.3 * planet_radius / num_layers
    cdef double first_step     = 0.0
    cdef size_t expected_size  = 50
    cdef double* t_eval_ptr    = NULL
    cdef size_t len_t_eval     = 0
    
    # Pressure convergence variables
    cdef size_t surface_pressure_index   = 0
    cdef double calculated_surf_pressure = d_INF_DBL
    cdef double pressure_diff            = d_INF_DBL
    cdef unsigned int iterations         = 0
    cdef cpp_bool failed                 = False
    cdef cpp_bool max_iters_hit          = False
    cdef cpp_bool final_run              = False
    
    # Integration solution variables
    cdef size_t last_solution_size = 0
    cdef CySolveOutput integration_result
    cdef CySolverResult* integration_result_ptr = NULL
    cdef EOSSolutionVec layer_solutions = EOSSolutionVec(0)
    layer_solutions.reserve(num_layers)
    
    while True:
        
        if not final_run:
            iterations += 1
        # Step through each macro layer of the planet and solve for density and gravity
        for layer_i in range(num_layers):
            # Setup bounds and initial conditions for next layer's integration
            radial_span_ptr[1] = layer_upper_radii[layer_i]
            if layer_i == 0:
                radial_span_ptr[0] = 0.0
                # Set y0 for bottom-most layer equal to the global y0
                y0_layer_ptr[0] = y0_ptr[0]
                y0_layer_ptr[1] = y0_ptr[1]
            else:
                radial_span_ptr[0] = layer_upper_radii[layer_i - 1]
                top_of_last_layer_index = (num_extra + num_y) * (last_solution_size - 1)
                # y0 for this layer equals result of last layer
                if integration_result_ptr:
                    y0_layer_ptr[0] = integration_result_ptr.solution[top_of_last_layer_index]
                    y0_layer_ptr[1] = integration_result_ptr.solution[top_of_last_layer_index + 1]
                else:
                    # Not sure why that would be null but in any case we are in a fail state.
                    failed = True
                    break
                
            # Get eos function and inputs for this layer
            eos_input_layer_ptr = eos_input_bylayer_ptrs[layer_i]
            eos_function_ptr    = eos_function_bylayer_ptrs[layer_i]
            
            if final_run:
                # We now want to make sure that all final calculations are performed.
                eos_input_layer_ptr.update_bulk  = True
                eos_input_layer_ptr.update_shear = True
                eos_input_layer_ptr.final_solve  = True
                # Capture extra outputs and store interpolators
                num_extra = 5
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
            args_ptr = <void*>eos_input_layer_ptr

            ###### Radial Integrate the EOS Through the Planet ######
            integration_result = cysolve_ivp(
                eos_diffeq,          # Differential equation [DiffeqFuncType]
                radial_span_ptr,     # Radial span [const double*]
                y0_layer_ptr,        # y0 array [const double*]
                num_y,               # Integration method [unsigned int]
                integration_method,  # Integration method [unsigned int]
                rtol,                # Relative Tolerance (as scalar) [double]
                atol,                # Absolute Tolerance (as scalar) [double]
                args_ptr,            # Extra input args to diffeq [void*]
                num_extra,           # Number of extra outputs tracked [unsigned int]
                max_num_steps,       # Max number of steps (0 = find good value) [size_t]
                max_ram_MB,          # Max amount of RAM allowed [size_t]
                use_dense_output,    # Use dense output [bint]
                t_eval_ptr,          # Interpolate at radius array [double*]
                len_t_eval,          # Size of interpolation array [size_t]
                eos_function_ptr,    # Pre-eval function used in diffeq [PreEvalFunc]
                rtols_ptr,           # Relative Tolerance (as array) [double*]
                atols_ptr,           # Absolute Tolerance (as array) [double*]
                max_step,            # Maximum step size [double]
                first_step,          # Initial step size (0 = find good value) [doub;e]
                expected_size        # Expected final integration size (0 = find good value) [size_t]
                )
            integration_result_ptr = integration_result.get()
            #########################################################
            last_solution_size = integration_result_ptr.size

            if not integration_result_ptr.success:
                failed = True
                break
            
            if final_run:
                layer_solutions.push_back(integration_result_ptr)

        if failed:
            break
        
        if final_run:
            # We are done!
            break
        else:
            surface_pressure_index = num_y * last_solution_size - 1
            calculated_surf_pressure = integration_result_ptr.solution[surface_pressure_index]

            # Update the centeral pressure using the error at the surface as the correction factor
            pressure_diff = surface_pressure - calculated_surf_pressure
            pressire_diff_abs = pressure_diff
            if pressure_diff < 0.0:
                pressire_diff_abs = -pressure_diff
                
            y0_ptr[1] += pressure_diff

            # Calculate percent difference to use in convergence check.
            if surface_pressure != 0.0:
                pressire_diff_abs /= surface_pressure
                
            # Check if we are done next iteration
            if pressire_diff_abs <= pressure_tol:
                final_run = True
        
        if iterations >= max_iters:
            max_iters_hit = True
            # To ensure that there is some output we will go ahead and do a final run.
            final_run = True
    
    # Display any warnings
    if max_iters_hit:
        strcpy(local_message_ptr, "Warning in `solve_eos`: Maximum number of iterations hit without convergence.\n")
        if verbose:
            printf(local_message_ptr)
    
    if failed:
        success_ptr[0] = False
        if integration_result_ptr:
            sprintf(local_message_ptr, "Warning in `solve_eos`: Integrator failed at iteration %d. Message: %s\n", iterations, integration_result_ptr.message_ptr)
        else:
            sprintf(local_message_ptr, "Warning in `solve_eos`: Integrator failed at iteration %d.\n", iterations)
        if verbose:
            printf(local_message_ptr)
    else:
        success_ptr[0] = True
    strcpy(message_ptr, local_message_ptr)

    return layer_solutions
