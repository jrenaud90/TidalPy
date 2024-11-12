from libc.stdio cimport printf
from libc.string cimport memcpy
from libcpp import bool as cpp_bool

from CyRK cimport cysolve_ivp, CySolverResult

from TidalPy.Material.eos.ode cimport eos_diffeq
from TidalPy.utilities.math.complex cimport cf_build_dblcmplx
from TidalPy.utilities.constants_x cimport G, PI_DBL, INF_DBL, NAN_DBL

cdef EOSSolutionVec solve_eos(
        double* radius_array_ptr,
        size_t len_radius_array,
        double* layer_upper_radii,
        unsigned int num_layers,
        PreEvalFunc* eos_function_bylayer_ptrs,
        EOS_ODEInput** eos_input_bylayer_ptrs,
        double planet_bulk_density,
        double surface_pressure = 0.0,
        double G_to_use = G,
        unsigned int integration_method = 2,
        double rtol = 1.0e-8,
        double atol = 1.0e-16,
        double pressure_tol = 1.0e-3,
        unsigned int max_iters = 100
        ) noexcept nogil:

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
        (2. / 3.) * PI_DBL * G_to_use * planet_radius * planet_radius * planet_bulk_density * planet_bulk_density
    
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
    cdef double calculated_surf_pressure = INF_DBL
    cdef double pressure_diff            = INF_DBL
    cdef unsigned int iterations         = 0
    cdef cpp_bool failed                 = False
    cdef cpp_bool max_iters_hit          = False
    cdef cpp_bool final_run              = False
    
    # Integration solution variables
    cdef size_t last_solution_size = 0
    cdef CySolveOutput solution
    cdef CySolverResult* solution_ptr = NULL
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
                y0_layer_ptr[0] = solution_ptr.solution[top_of_last_layer_index]
                y0_layer_ptr[1] = solution_ptr.solution[top_of_last_layer_index + 1]
                
            # Get eos function and inputs for this layer
            eos_input_layer_ptr = eos_input_bylayer_ptrs[layer_i]
            eos_function_ptr    = eos_function_bylayer_ptrs[layer_i]
            
            # Convert input pointer to void pointer (required by cysolve)
            args_ptr = <void*>eos_input_layer_ptr
            
            if final_run:
                # We now want to make sure that all final calculations are performed.
                eos_input_layer_ptr.update_bulk = True
                eos_input_layer_ptr.update_shear = True
                eos_input_layer_ptr.final_solve = True
                # Capture extra outputs and store interpolators
                num_extra = 5
                use_dense_output = True
            else:
                # During the iterations we do not need to update the complex bulk or shear
                eos_input_layer_ptr.update_bulk = False
                eos_input_layer_ptr.update_shear = False
                # We also are not at the final call step. 
                eos_input_layer_ptr.final_solve = False
                num_extra = 0
                use_dense_output = False
            
            solution = cysolve_ivp(
                eos_solution, radial_span_ptr, y0_layer_ptr, num_y,
                integration_method, rtol, atol, args_ptr, num_extra,
                max_num_steps, max_ram_MB, use_dense_output, t_eval_ptr, len_t_eval, eos_function_ptr,
                rtols_ptr, atols_ptr, max_step, first_step, expected_size)
            solution_ptr = solution.get()
            last_solution_size = solution_ptr.size

            if not solution_ptr.success:
                failed = True
                break
            
            if final_run:
                layer_solutions.push_back(solution)

        if failed:
            break
        
        if final_run:
            # We are done!
            break
        else:
            surface_pressure_index = num_y * last_solution_size - 1
            calculated_surf_pressure = solution_ptr.solution[surface_pressure_index]

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
        printf("Warning in `solve_eos`: Maximum number of iterations hit without convergence.\n")
    
    if failed:
        if solution_ptr:
            printf("Warning in `solve_eos`: Integrator failed at iteration %d. Message: %s\n", iterations, solution_ptr.message_ptr)
        else:
            printf("Warning in `solve_eos`: Integrator failed at iteration %d.\n", iterations)
    
    return layer_solutions
