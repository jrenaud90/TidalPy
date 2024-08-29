from libc.stdio cimport printf
from libc.string cimport memcpy
from libcpp import bool as cpp_bool

from CyRK cimport cysolve_ivp, CySolveOutput, PreEvalFunc, CySolverResult

from TidalPy.Material.eos.ode cimport eos_solution
from TidalPy.utilities.math.complex cimport cf_build_dblcmplx
from TidalPy.utilities.constants_x cimport G, PI_DBL, INF_DBL, NAN_DBL

cdef CySolveOutput solve_eos(
        double* radius_array_ptr,
        size_t len_radius_array,
        double* layer_upper_radii,
        unsigned int num_layers,
        PreEvalFunc* eos_function_ptrs,
        EOS_ODEInput* eos_input_ptrs,
        double planet_bulk_density,
        double surface_pressure = 0.0,
        double G_to_use = G,
        unsigned int integration_method = 2,
        double rtol = 1.0e-6,
        double atol = 1.0e-8,
        double pressure_tol = 1.0e-3,
        unsigned int max_iters = 50
        ) noexcept nogil:

    # Don't use rtol or atol arrays
    cdef double* rtols_ptr = NULL
    cdef double* atols_ptr = NULL

    # Determine bounds
    cdef double planet_radius
    cdef double[2] radial_span
    cdef double* radial_span_ptr = &radial_span[0]
    radial_span_ptr[0] = radius_array_ptr[0]
    radial_span_ptr[1] = radius_array_ptr[len_radius_array - 1]
    planet_radius = radial_span_ptr[1]

    # We need the centeral pressure of the planet. Use the global bulk modulus to calculate this.
    cdef double centeral_gravity = 0.0
    cdef double centeral_pressure_guess = (4. / 3.) * PI_DBL * G_to_use * planet_radius * planet_radius * planet_bulk_density * planet_bulk_density
    cdef double[2] y0 = [centeral_gravity, centeral_pressure_guess]
    cdef double* y0_ptr = &y0[0]

    # y information
    cdef unsigned int num_y = 2
    cdef unsigned int num_extra = 0

    # Modify additional arguments to work with integrator
    cdef EOS_ODEInput eos_input
    cdef EOS_ODEInput* eos_input_ptr = &eos_input
    cdef void* args_ptr = <void*>eos_input

    # Other integration information
    cdef size_t max_num_steps = 10_000
    cdef size_t max_ram_MB = 500
    cdef bint use_dense_output = False
    cdef double max_step = 0.25 * planet_radius
    cdef double first_step = 0.0
    cdef size_t expected_size = 50

    # Solution information
    cdef double* t_eval_ptr = NULL
    cdef size_t len_t_eval = 0
    cdef CySolveOutput solution
    cdef CySolverResult* solution_ptr = NULL
    cdef double calculated_surf_pressure = INF_DBL
    cdef double pressure_diff = INF_DBL
    cdef size_t surface_pressure_index = 0
    cdef size_t iterations = 0
    cdef cpp_bool failed = False
    cdef cpp_bool max_iters_hit = False

    # Layer information
    cdef unsigned int layer_i
    cdef double[2] y0_layer = [centeral_gravity, centeral_pressure_guess]
    cdef double* y0_layer_ptr = &y0_layer[0]


    # Begin iterating until we find convergence on the surface pressure.
    while (pressure_diff > pressure_tol) and (iterations < max_iters):

        # Step through each macro layer of the planet and solve for density and gravity
        for layer_i in range(num_layers):

            radial_span_ptr[1] = layer_upper_radii[layer_i]

            if layer_i == 0:
                # Set y0 equal to the global y0
                y0_layer_ptr[0] = y0_ptr[0]
                y0_layer_ptr[1] = y0_ptr[1]
            
            # Copy over the contents of the eos arg input so values can be manipulated
            memcpy(eos_input_ptr, <void*>eos_input_ptrs[layer_i], sizeof(EOS_ODEInput))
            # During the iterations we do not need to update the complex bulk or shear
            eos_input_ptr.update_bulk = False
            eos_input_ptr.update_shear = False
            # We also are not at the final call step. 
            eos_input_ptr.final_solve = False

            solution = cysolve_ivp(
                eos_solution, radial_span_ptr, y0_ptr, num_y, integration_method, rtol, atol, args_ptr, num_extra,
                max_num_steps, max_ram_MB, use_dense_output, t_eval_ptr, len_t_eval, eos_function_ptrs[layer_i],
                rtols_ptr, atols_ptr, max_step, first_step, expected_size)
            solution_ptr = solution.get()

            radial_span_ptr[0] = radial_span_ptr[1]

            if not solution_ptr.success:
                failed = True
                break

        if failed:
            break

        surface_pressure_index = 2 * solution_ptr.size
        calculated_surf_pressure = solution_ptr.solution[surface_pressure_index]

        # Update the centeral pressure using the error at the surface as the correction factor
        pressure_diff = surface_pressure - calculated_surf_pressure
        y0_ptr[1] += pressure_diff

        # Calculate percent difference to use in loop check.
        if pressure_diff < 0.0:
            pressure_diff = -pressure_diff

        if surface_pressure != 0.0:
            pressure_diff /= surface_pressure

        iterations += 1

        if iterations > max_iters:
            max_iters_hit = True
            break

    # Display any warnings
    if max_iters_hit:
        printf("Warning in `solve_eos`: Maximum number of iterations hit without convergence.")
    
    if failed:
        if solution_ptr:
            printf("Warning in `solve_eos`: Integrator failed at iteration %d. Message: %s", iterations, solution_ptr.message_ptr)
        else:
            printf("Warning in `solve_eos`: Integrator failed at iteration %d.", iterations)

    # Solve one last time using the provided radius array as our target steps
    radial_span_ptr[0] = radius_array_ptr[0]

    # Parameters to set t_eval
    cdef double radius_check
    cdef cpp_bool found_start = False
    cdef cpp_bool found_end = False
    cdef size_t slice_i = 0
    cdef double start_radius = 0.0
    cdef double end_radius = 0.0
    cdef size_t t_eval_start_index = 0
    cdef size_t num_t_eval = 0

    # Also save additional data calculated during integration
    num_extra = 5

    # Step through each macro layer of the planet. Store solutions in a vector
    cdef vector[CySolveOutput] layer_solutions = vector[CySolveOutput](0)
    layer_solutions.reserve(num_layers)

    for layer_i in range(num_layers):

        if layer_i == 0:
            start_radius = 0.0
        else:
            start_radius = layer_upper_radii[layer_i - 1]
        if layer_i == (num_layers - 1):
            end_radius = planet_radius
        else:
            end_radius = layer_upper_radii[layer_i + 1]
        
        found_start = False
        found_end = False
        num_t_eval = 0
        for slice_i in range(len_radius_array):
            if found_start and found_end:
                break
            radius_check = radius_array_ptr[slice_i]

            if radius_check >= start_radius:
                
                if not found_start:
                    t_eval_ptr  = &radius_array_ptr[slice_i]
                    found_start = True
                
                if radius_check < end_radius:
                    # Still in layer
                    num_t_eval += 1
                elif radius_check == end_radius:
                    # At very top of layer. If the last layer then save this data point. Otherwise push it to next layer.
                    if layer_i == (num_layers - 1):
                        num_t_eval += 1
                        found_end = True
                    else:
                        found_end = True
                else:
                    found_end = True

        radial_span_ptr[1] = layer_upper_radii[layer_i]

        if layer_i == 0:
            # Set y0 equal to the global y0
            y0_layer_ptr[0] = y0_ptr[0]
            y0_layer_ptr[1] = y0_ptr[1]
        
        # Copy over the contents of the eos arg input so values can be manipulated
        memcpy(eos_input_ptr, <void*>eos_input_ptrs[layer_i], sizeof(EOS_ODEInput))

        # We now want to make sure that all final calculations are performed.
        eos_input_ptr.final_solve = True

        solution = cysolve_ivp(
            eos_solution, radial_span_ptr, y0_ptr, num_y, integration_method, rtol, atol, args_ptr, num_extra,
            max_num_steps, max_ram_MB, use_dense_output, t_eval_ptr, len_t_eval, eos_function_ptrs[layer_i],
            rtols_ptr, atols_ptr, max_step, first_step, expected_size)

        layer_solutions.push_back(solution)

        radial_span_ptr[0] = radial_span_ptr[1]
    
    return layer_solutions
