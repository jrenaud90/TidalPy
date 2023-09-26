# distutils: language = c++

from libcpp cimport bool as bool_cpp_t
from libc.math cimport NAN, fmax

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

import numpy as np
cimport numpy as np

from scipy.constants import G as G_
from TidalPy.radial_solver.nondimensional import re_dimensionalize_radial_func
from TidalPy.radial_solver.numerical.collapse import collapse_solutions
from TidalPy.radial_solver.numerical.initial import find_initial_guess
from TidalPy.radial_solver.numerical.interfaces import find_interface_func

# Import cythonized functions
from TidalPy.radial_solver.nondimensional_x import non_dimensionalize_physicals_x
from TidalPy.radial_solver.numerical.interfaces.interfaces_x cimport find_solution_num
from TidalPy.radial_solver.numerical.derivatives.ode_base_x cimport RadialSolverBase
from TidalPy.radial_solver.numerical.derivatives.odes_x cimport build_solver



# Setup globals
cdef double G
G = G_

cdef Py_ssize_t MAX_NUM_SOLS    = 3
cdef Py_ssize_t MAX_NUM_YS      = 6
cdef Py_ssize_t MAX_NUM_YS_REAL = 2 * MAX_NUM_YS

def radial_solver_x(
        double[:] radius_array,
        double[:] density_array,
        double[:] gravity_array,
        double[:] bulk_modulus_array,
        double complex[:] complex_shear_modulus_array,
        double frequency,
        double planet_bulk_density,
        tuple is_solid_by_layer,
        tuple is_static_by_layer,
        tuple is_incompressible_by_layer,
        tuple upper_radius_by_layer,
        unsigned int degree_l = 2,
        double complex[:] surface_boundary_conditions = None,
        bool_cpp_t solve_load_numbers = False,
        bool_cpp_t use_kamata = False,
        int integration_method = 1,
        double integration_rtol = 1.0e-4,
        double integration_atol = 1.0e-12,
        bool_cpp_t scale_rtols_by_layer_type = True,
        Py_ssize_t max_num_steps = 500_000,
        Py_ssize_t expected_size = 250,
        double max_step = 0,
        bool_cpp_t limit_solution_to_radius = True,
        bool_cpp_t verbose = False,
        bool_cpp_t nondimensionalize = True,
        ):
    # General indexing
    cdef Py_ssize_t i, j, k, m

    # Pull out key information
    cdef double radius_planet
    cdef Py_ssize_t num_layers, num_interfaces, total_slices

    total_slices   = len(radius_array)
    radius_planet  = radius_array[total_slices - 1]
    num_layers     = len(is_solid_by_layer)
    num_interfaces = num_layers - 1

    # Ensure there are enough slices for radial interpolations.
    if total_slices <= (3 * num_layers):
        raise AttributeError('Radial solver requires at least 3 radial slices per layer (ideally >= 10 per).')

    # Ensure there is at least one layer.
    if num_layers <= 0:
        raise AttributeError('Radial solver requires at least one layer.')

    # Build main input pointers
    cdef double* radius_array_ptr = <double *> PyMem_Malloc(total_slices * sizeof(double))
    if not radius_array_ptr:
        raise MemoryError()

    cdef double* density_array_ptr = <double *> PyMem_Malloc(total_slices * sizeof(double))
    if not density_array_ptr:
        raise MemoryError()

    cdef double* gravity_array_ptr = <double *> PyMem_Malloc(total_slices * sizeof(double))
    if not gravity_array_ptr:
        raise MemoryError()

    cdef double* bulk_array_ptr = <double *> PyMem_Malloc(total_slices * sizeof(double))
    if not bulk_array_ptr:
        raise MemoryError()

    cdef double complex* cmplx_shear_array_ptr = <double complex *> PyMem_Malloc(total_slices * sizeof(double complex))
    if not cmplx_shear_array_ptr:
        raise MemoryError()

    # Populate the arrays (making a copy of values)
    for i in range(total_slices):
        radius_array_ptr[i]      = radius_array[i]
        density_array_ptr[i]     = density_array[i]
        gravity_array_ptr[i]     = gravity_array[i]
        bulk_array_ptr[i]        = bulk_modulus_array[i]
        cmplx_shear_array_ptr[i] = complex_shear_modulus_array[i]


    # Non-dimensionalize inputs
    cdef double G_to_use, frequency_to_use
    if nondimensionalize:
        non_dimensionalize_physicals_x(
            total_slices,
            frequency,
            radius_planet,
            planet_bulk_density,
            radius_array_ptr,
            density_array_ptr,
            gravity_array_ptr,
            bulk_array_ptr,
            cmplx_shear_array_ptr,
            &frequency_to_use,
            &G_to_use
            )
    else:
        # Leave inputs alone.
        G_to_use         = G
        frequency_to_use = frequency

    # Find boundary condition at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    cdef double complex* bc_pointer = <double complex *> PyMem_Malloc(3 * sizeof(double complex))
    if not bc_pointer:
        raise MemoryError()

    if surface_boundary_conditions is None:
        # Assume tides
        for i in range(3):
            if i == 2:
                if nondimensionalize:
                    bc_pointer[i] == (2. * degree_l + 1.) / 1.
                else:
                    bc_pointer[i] == (2. * degree_l + 1.) / radius_planet
    else:
        # Use user input
        if len(surface_boundary_conditions) != 3:
            raise AttributeError('Unexpected number of user-provided surface boundary conditions.')
        for i in range(3):
            bc_pointer[i] = surface_boundary_conditions[i]

    # Integration information
    # Max step size
    cdef double max_step_touse
    cdef bool_cpp_t max_step_from_arrays
    max_step_from_arrays = False
    if max_step == 0:
        # If max_step is zero use the array information to determine max_step_size
        max_step_from_arrays = True
    else:
        # Otherwise use user input.
        max_step_touse = max_step

    # Tolerances
    cdef double* rtols_by_layer_ptr = <double *> PyMem_Malloc(num_layers * MAX_NUM_YS_REAL * sizeof(double))
    if not rtols_by_layer_ptr:
        raise MemoryError()

    cdef double* atols_by_layer_ptr = <double *> PyMem_Malloc(num_layers * MAX_NUM_YS_REAL * sizeof(double))
    if not atols_by_layer_ptr:
        raise MemoryError()

    # Find number of solutions per layer
    cdef Py_ssize_t* num_solutions_by_layer_ptr = \
        <Py_ssize_t*> PyMem_Malloc(num_layers * sizeof(Py_ssize_t))
    if not num_solutions_by_layer_ptr:
        raise MemoryError()

    cdef Py_ssize_t* start_index_by_layer_ptr = \
        <Py_ssize_t*> PyMem_Malloc(num_layers * sizeof(Py_ssize_t))
    if not start_index_by_layer_ptr:
        raise MemoryError()

    cdef Py_ssize_t* num_slices_by_layer_ptr = \
        <Py_ssize_t*> PyMem_Malloc(num_layers * sizeof(Py_ssize_t))
    if not num_slices_by_layer_ptr:
        raise MemoryError()

    cdef bool_cpp_t layer_is_solid, layer_is_static, layer_is_incomp
    cdef bool_cpp_t layer_below_is_solid, layer_below_is_static, layer_below_is_incomp
    cdef double layer_upper_radius, radius_check
    cdef double layer_rtol_real, layer_rtol_imag, layer_atol_real, layer_atol_imag
    cdef Py_ssize_t layer_slices, num_sols, num_ys, max_slices

    # Track the largest number of slices in a layer.
    # Used to more efficiently allocate memory for other storage arrays
    max_slices = -1

    for i in range(num_layers):
        # Pull out information on this layer
        layer_is_solid     = is_solid_by_layer[i]
        layer_is_static    = is_static_by_layer[i]
        layer_is_incomp    = is_incompressible_by_layer[i]
        layer_upper_radius = upper_radius_by_layer[i]

        # Find number of solutions based on this layer's assumptions
        num_sols = <Py_ssize_t>find_solution_num(
            layer_is_solid,
            layer_is_static,
            layer_is_incomp
            )
        num_ys = 2 * num_sols
        num_solutions_by_layer_ptr[i] = num_sols

        # Scale rtols by layer type
        for j in range(num_ys):
            # Default is that each layer's rtol and atol equal user input.
            layer_rtol_real = integration_rtol
            layer_rtol_imag = integration_rtol
            layer_atol_real = integration_atol
            layer_atol_imag = integration_atol

            if scale_rtols_by_layer_type:
                # Certain layer assumptions can affect solution stability so use additional scales on the relevant rtols
                # TODO test that these scales are the best.
                if layer_is_solid:
                    # Scale y2 and y3 by 0.1
                    if (j == 1) or (j == 2):
                        # Scale both the real and imaginary portions by the same amount.
                        layer_rtol_real *= 0.1
                        layer_rtol_imag *= 0.1
                else:
                    if not layer_is_static:
                        # Scale dynamic liquid layer's y2 by additional 0.01.
                        if (j == 1):
                            # Scale both the real and imaginary portions by the same amount.
                            layer_rtol_real *= 0.01
                            layer_rtol_imag *= 0.01
            # Populate rtol and atol pointers.
            rtols_by_layer_ptr[i * MAX_NUM_YS_REAL + (2 * j)]     = layer_rtol_real
            rtols_by_layer_ptr[i * MAX_NUM_YS_REAL + (2 * j + 1)] = layer_rtol_imag
            atols_by_layer_ptr[i * MAX_NUM_YS_REAL + (2 * j)]     = layer_atol_real
            atols_by_layer_ptr[i * MAX_NUM_YS_REAL + (2 * j + 1)] = layer_atol_imag

        # Determine how many slices are in this layer
        if i == 0:
            # First layer starts at 0.
            start_index_by_layer_ptr[0] = 0
        else:
            # Not first layer. Starting point is based on number of slices in previous layer.
            start_index_by_layer_ptr[i] = start_index_by_layer_ptr[i - 1] + num_slices_by_layer_ptr[i - 1]

        layer_slices = 0
        for j in range(start_index_by_layer_ptr[i], total_slices):
            radius_check = radius_array[j]
            if radius_check > layer_upper_radius:
                # We have passed this layer.
                break
            layer_slices += 1
        if layer_slices <= 3:
            raise ValueError('At least three layer slices per layer are required\n\tTry using more slices in the'
                             'input arrays.')
        num_slices_by_layer_ptr[i] = layer_slices
        if max_slices < layer_slices:
            max_slices = layer_slices

    # We have all the size information needed to build storage pointers
    # Main storage pointer is setup like [layer_i][solution_i][y_i + r_i]
    cdef double complex*** main_storage = <double complex ***> PyMem_Malloc(num_layers * sizeof(double complex**))
    if not main_storage:
        raise MemoryError()

    cdef double complex** storage_by_solution
    cdef double complex* storage_by_y

    for i in range(num_layers):
        num_sols     = num_solutions_by_layer_ptr[i]
        layer_slices = num_slices_by_layer_ptr[i]
        # Number of ys = 2x num sols
        num_ys = 2 * num_sols

        storage_by_solution = <double complex**> PyMem_Malloc(num_sols * sizeof(double complex*))
        if not storage_by_solution:
            raise MemoryError()

        for j in range(num_sols):
            storage_by_y = <double complex*> PyMem_Malloc(layer_slices * num_ys * sizeof(double complex))
            if not storage_by_y:
                raise MemoryError()

            storage_by_solution[j] = storage_by_y
        main_storage[i] = storage_by_solution

    # Storage for solutions at the top of a layer.
    cdef double complex[:, ::1] last_layer_top_solutions_view

    # Interface information
    cdef object interface_func
    cdef tuple interface_inputs

    # Layer specific pointers; set the size based on the layer with the most slices.
    cdef double* layer_radius_ptr
    cdef double* layer_density_ptr
    cdef double* layer_gravity_ptr
    cdef double* layer_bulk_mod_ptr
    cdef double complex* layer_shear_mod_ptr
    cdef double* layer_rtols_ptr
    cdef double* layer_atols_ptr

    # Properties at top and bottom of layer
    cdef double radius_lower, radius_upper
    cdef double density_lower, density_upper
    cdef double gravity_lower, gravity_upper
    cdef double bulk_lower
    cdef double complex shear_lower
    cdef tuple radial_span

    # Properties at interfaces between layers
    cdef double static_liquid_density
    cdef double interface_gravity
    cdef double last_layer_upper_gravity, last_layer_upper_density

    # Starting solutions (initial conditions / lower boundary conditions)
    cdef double complex[:, ::1] initial_solutions_view
    cdef double[:, ::1] initial_solutions_real_view

    # Solver class
    cdef RadialSolverBase solver

    # Feedback
    cdef str feedback_str
    cdef bool_cpp_t error
    feedback_str = 'Starting integration'
    if verbose:
        print(feedback_str)
    error = False

    cdef Py_ssize_t start_index, end_index
    for i in range(num_layers):
        # Get layer's index information
        start_index  = start_index_by_layer_ptr[i]
        layer_slices = num_slices_by_layer_ptr[i]
        end_index    = start_index + (layer_slices - 1)

        # Get solution and y information
        num_sols = num_solutions_by_layer_ptr[i]
        num_ys   = 2 * num_sols

        # Setup pointer array slices starting at this layer's beginning
        layer_radius_ptr    = &radius_array_ptr[start_index]
        layer_density_ptr   = &density_array_ptr[start_index]
        layer_gravity_ptr   = &gravity_array_ptr[start_index]
        layer_bulk_mod_ptr  = &bulk_array_ptr[start_index]
        layer_shear_mod_ptr = &cmplx_shear_array_ptr[start_index]

        # Setup other pointers
        layer_rtols_ptr = &rtols_by_layer_ptr[i * MAX_NUM_YS_REAL]
        layer_atols_ptr = &atols_by_layer_ptr[i * MAX_NUM_YS_REAL]

        # Get physical parameters at the top and bottom of the layer
        radius_lower  = layer_radius_ptr[0]
        density_lower = layer_density_ptr[0]
        gravity_lower = layer_gravity_ptr[0]
        bulk_lower    = layer_bulk_mod_ptr[0]
        shear_lower   = layer_shear_mod_ptr[0]
        radius_upper  = layer_radius_ptr[layer_slices - 1]
        density_upper = layer_density_ptr[layer_slices - 1]
        gravity_upper = layer_gravity_ptr[layer_slices - 1]

        # Determine max step size (if not provided by user)
        if max_step_from_arrays:
            # Maximum step size during integration can not exceed the average radial slice size.
            max_step_touse = (radius_upper - radius_lower) / <double>layer_slices

        # Get assumptions for layer
        layer_is_solid  = is_solid_by_layer[i]
        layer_is_static = is_static_by_layer[i]
        layer_is_incomp = is_incompressible_by_layer[i]

        # Find initial conditions for each solution at the base of this layer.
        radial_span = (radius_lower, radius_upper)
        if i == 0:
            # Use initial condition function
            initial_solutions_view = find_initial_guess(
                layer_is_solid,
                layer_is_static,
                layer_is_incomp,
                use_kamata,
                radius_lower,
                shear_lower,
                bulk_lower,
                density_lower,
                frequency_to_use,
                order_l=degree_l,
                G_to_use=G_to_use
                )
        else:
            layer_below_is_solid  = is_solid_by_layer[i - 1]
            layer_below_is_static = is_static_by_layer[i - 1]
            layer_below_is_incomp = is_incompressible_by_layer[i - 1]

            # Find gravity at the base interface using bottom of this layer and top of previous.
            interface_gravity = 0.5 * (gravity_lower + last_layer_upper_gravity)

            # Find the density needed for some initial conditions.
            if layer_is_solid and layer_below_is_solid:
                # Both layers are solid. A liquid interface density is not needed.
                static_liquid_density = NAN
            elif not layer_is_solid and layer_below_is_solid:
                # Layer below is solid, this layer is liquid. Use its density.
                static_liquid_density = density_lower
            elif layer_is_solid and not layer_below_is_solid:
                # Layer below is liquid, this layer is solid. Use layer below's density.
                static_liquid_density = last_layer_upper_density
            else:
                # Both layers are liquid. Choose the layer's density which is static.
                if layer_is_static and layer_below_is_static:
                    # Both layers are static.
                    # TODO: Not sure what to do here so just use this layer's density.
                    static_liquid_density = density_lower
                elif layer_is_static and not layer_below_is_static:
                    # This layer is static, one below is not. Use this layer's density
                    static_liquid_density = density_lower
                elif not layer_is_static and layer_below_is_static:
                    # This layer is dynamic, layer below is static. Use layer below's density
                    static_liquid_density = last_layer_upper_density
                else:
                    # Both layers are dynamic. Static liquid density is not needed.
                    static_liquid_density = NAN

            # Use results at the top of the layer below and a interface function
            interface_func, interface_inputs = find_interface_func(
                layer_below_is_solid,
                layer_below_is_static,
                layer_is_solid,
                layer_is_static,
                static_liquid_density,
                interface_gravity,
                G_to_use=G_to_use
                )

            # Run interface function
            initial_solutions_view = interface_func(
                last_layer_top_solutions_view
                *interface_inputs
                )

        # Change initial conditions into 2x real values instead of complex for integration
        # TODO: change to pointers
        initial_solutions_real_arr  = np.empty((num_sols, num_ys * 2), dtype=np.float64, order='C')
        initial_solutions_real_view = initial_solutions_real_arr
        for j in range(num_sols):
            for k in range(num_ys):
                initial_solutions_real_view[j, (2 * k)]     = np.real(initial_solutions_view[j, k])
                initial_solutions_real_view[j, (2 * k) + 1] = np.imag(initial_solutions_view[j, k])

        # Build solver instance
        solver = build_solver(
            layer_is_solid,
            layer_is_static,
            layer_is_incomp,
            layer_slices,
            layer_radius_ptr,
            layer_density_ptr,
            layer_gravity_ptr,
            layer_bulk_mod_ptr,
            layer_shear_mod_ptr,
            frequency_to_use,
            degree_l,
            G_to_use,
            radial_span,
            initial_solutions_real_view[0, :],
            layer_atols_ptr,
            layer_rtols_ptr,
            integration_method,
            max_step_touse,
            max_num_steps,
            expected_size,
            limit_solution_to_radius
            )

        # Get storage pointer for this layer
        storage_by_solution = main_storage[i]

        # TODO: make this pointers
        last_layer_top_solutions_arr  = np.empty((num_sols, num_ys), dtype=np.complex128, order='C')
        last_layer_top_solutions_view = last_layer_top_solutions_arr

        # Solve for each solution
        for j in range(num_sols):

            if j > 0:
                # Reset solver with new initial condition (this is already done for j==0)
                solver.change_y0(initial_solutions_real_view[j, :], auto_reset_state=False)

            # Integrate!
            solver._solve(reset=True)

            # Check for problems
            if not solver.success:
                # Problem with integration.
                feedback_str = f'Integration problem at layer {i}; solution {j}:\n\t{solver.message}'
                if verbose:
                    print(feedback_str)
                error = True
                break

            # If no problems, store results.
            # Need to make a copy because the solver pointers will be reallocated during the next solution.
            # Get storage pointer for this solution
            storage_by_y = storage_by_solution[j]
            for m in range(layer_slices):
                for k in range(num_ys):

                    # Convert 2x real ys to 1x complex ys
                    storage_by_y[m * num_ys + k] = \
                        (solver.solution_y_ptr[m * num_ys + (2 * k)] +
                         1.0j * solver.solution_y_ptr[m * num_ys + (2 * k) + 1])

                    # Store top most result for initial condition for the next layer
                    if m == (layer_slices - 1):
                        last_layer_top_solutions_view[j, k] = storage_by_y[m * num_ys + k]

        if error:
            # Error was encountered during integration
            break

        # Prepare for next layer
        last_layer_upper_gravity = gravity_upper
        last_layer_upper_density = density_upper

    cdef np.ndarray[np.complex128_t, ndim=2] full_solution_arr
    cdef double complex[:, ::1] full_solution_view
    # No matter the results of the integration, we know the shape and size of the final solution
    if solve_load_numbers:
        full_solution_arr = np.empty((6, total_slices), dtype=np.complex128, order='C')
    else:
        full_solution_arr = np.empty((12, total_slices), dtype=np.complex128, order='C')
    full_solution_view = full_solution_arr

    if error:
        feedback_str = 'Integration failed.'
        if verbose:
            print(feedback_str)
        for i in range(12):
            for j in range(total_slices):
                full_solution_view[i, j] = NAN
    else:
        feedback_str = 'Integration completed for all layers. Beginning solution collapse.'
        if verbose:
            print(feedback_str)
        # Collapse the multiple solutions for each layer into one final combined solution.

    # Free memory
    del solver

    # Release planet pointers
    PyMem_Free(radius_array_ptr)
    PyMem_Free(density_array_ptr)
    PyMem_Free(gravity_array_ptr)
    PyMem_Free(bulk_array_ptr)
    PyMem_Free(cmplx_shear_array_ptr)

    # Release integration-related pointers
    PyMem_Free(rtols_by_layer_ptr)
    PyMem_Free(atols_by_layer_ptr)

    # Release BC-related pointers
    PyMem_Free(bc_pointer)

    # Deconstruct main solution pointer
    # Main storage pointer is setup like [layer_i][solution_i][y_i + r_i]
    for i in range(num_layers):
        num_sols = num_solutions_by_layer_ptr[i]
        for j in range(num_sols):
            PyMem_Free(main_storage[i][j])
        PyMem_Free(main_storage[i])
    PyMem_Free(main_storage)

    # Release layer information pointers
    PyMem_Free(num_solutions_by_layer_ptr)
    PyMem_Free(start_index_by_layer_ptr)
    PyMem_Free(num_slices_by_layer_ptr)

    return full_solution_arr
