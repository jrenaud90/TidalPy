# distutils: language = c++

from libcpp cimport bool as bool_cpp_t
from libc.math cimport NAN, fmax

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

import numpy as np
cimport numpy as np

from scipy.constants import G as G_
from TidalPy.radial_solver.nondimensional import non_dimensionalize_physicals, re_dimensionalize_radial_func
from TidalPy.radial_solver.numerical.collapse import collapse_solutions
from TidalPy.radial_solver.numerical.initial import find_initial_guess
from TidalPy.radial_solver.numerical.interfaces import find_interface_func

# Import cythonized functions
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
        bool_cpp_t scale_rtols_by_layer_ptr_type = True,
        Py_ssize_t max_num_steps = 500_000,
        Py_ssize_t expected_size = 850,
        double max_step = 0,
        bool_cpp_t limit_solution_to_radius = True,
        bool_cpp_t verbose = False,
        bool_cpp_t nondimensionalize = True,
        ):
    # General indexing
    cdef Py_ssize_t i, j, k, m

    # Pull out key information
    cdef double radius_planet
    cdef Py_ssize_t num_layers, num_interfaces, num_radius

    num_radius     = len(radius_array)
    radius_planet  = radius_array[num_radius - 1]
    num_layers     = len(is_solid_by_layer)
    num_interfaces = num_layers - 1
    if num_layers <= 0:
        raise AttributeError('Radial solver requires at least one layer.')

    # Non-dimensionalize inputs
    cdef double G_to_use

    if nondimensionalize:
        radius_array, gravity_array, density_array, complex_shear_array, bulk_array, frequency, G_to_use = \
            non_dimensionalize_physicals(
                    radius_array, gravity_array, density_array, complex_shear_modulus_array, bulk_modulus_array, frequency,
                    mean_radius=radius_planet, bulk_density=planet_bulk_density
                    )
    else:
        G_to_use = G

    # Find boundary condition at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    cdef np.ndarray[np.complex128_t, ndim=1] bc_array
    cdef double complex[::1] bc_array_view
    bc_array = np.empty(3, dtype=np.complex128, order='C')
    bc_array_view = bc_array

    if surface_boundary_conditions is None:
        # Assume tides
        for i in range(3):
            if i == 2:
                if nondimensionalize:
                    bc_array_view[i] == (2. * degree_l + 1.) / 1.
                else:
                    bc_array_view[i] == (2. * degree_l + 1.) / radius_planet
    else:
        # Use user input
        if len(surface_boundary_conditions) != 3:
            raise AttributeError('Unexpected number of user-provided surface boundary conditions.')
        for i in range(3):
            bc_array_view[i] = surface_boundary_conditions[i]

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
    cdef Py_ssize_t num_slices, num_sols, num_ys, max_slices

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
            atols_by_layer_ptr[i * MAX_NUM_YS_REAL + j] = integration_atol
            # TODO test that these scales are the best.
            if scale_rtols_by_layer_ptr_type:
                if layer_is_solid:
                    # Scale (real and imag) y2 (index 2, 3) and y3 (index 4, 5) by 0.1
                    if (j >= 2) and (j <= 5):
                        rtols_by_layer_ptr[i * MAX_NUM_YS_REAL + j] = integration_rtol * 0.1
                    else:
                        rtols_by_layer_ptr[i * MAX_NUM_YS_REAL + j] = integration_rtol
                else:
                    if layer_is_static:
                        # Don't scale
                        rtols_by_layer_ptr[i * MAX_NUM_YS_REAL + j] = integration_rtol
                    else:
                        # Scale (real and imag) y2 by 0.01
                        if (j == 2) or (j == 3):
                            rtols_by_layer_ptr[i * MAX_NUM_YS_REAL + j] = integration_rtol * 0.01
                        else:
                            rtols_by_layer_ptr[i * MAX_NUM_YS_REAL + j] = integration_rtol

        # Determine how many slices are in this layer
        if i == 0:
            # First layer starts at 0.
            start_index_by_layer_ptr[0] = 0
        else:
            # Not first layer. Starting point is based on number of slices in previous layer.
            start_index_by_layer_ptr[i] = start_index_by_layer_ptr[i - 1] + num_slices_by_layer_ptr[i - 1]

        num_slices = 0
        for j in range(start_index_by_layer_ptr[i], num_radius):
            radius_check = radius_array[j]
            if radius_check > layer_upper_radius:
                # We have passed this layer.
                break
            num_slices += 1
        num_slices_by_layer_ptr[i] = num_slices
        if max_slices < num_slices:
            max_slices = num_slices

    # We have all the size information needed to build storage pointers
    # Main storage pointer is setup like [layer_i][solution_i][y_i + r_i]
    cdef double complex*** main_storage = <double complex ***> PyMem_Malloc(num_layers * sizeof(double complex**))
    if not main_storage:
        raise MemoryError()

    cdef double complex** storage_by_solution
    cdef double complex* storage_by_y

    for i in range(num_layers):
        num_sols   = num_solutions_by_layer_ptr[i]
        num_slices = num_slices_by_layer_ptr[i]
        # Number of ys = 2x num sols
        num_ys = 2 * num_sols

        storage_by_solution = <double complex**> PyMem_Malloc(num_sols * sizeof(double complex*))
        if not storage_by_solution:
            raise MemoryError()

        for j in range(num_sols):
            storage_by_y = <double complex*> PyMem_Malloc(num_slices * num_ys * sizeof(double complex))
            if not storage_by_y:
                raise MemoryError()

            storage_by_solution[j] = storage_by_y
        main_storage[i] = storage_by_solution

    # Storage for solutions at the top of a layer.
    cdef double complex[:, :] last_layer_top_solutions_view

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

    # Starting solutions (initial conditions / lower boundary conditions)
    cdef double complex[:, ::1] initial_solutions_view
    cdef double[:, ::1] initial_solutions_real_view

    # Solver class
    cdef RadialSolverBase solver

    # Feedback
    cdef str feedback_str
    feedback_str = ''
    cdef bool_cpp_t error
    error = False

    cdef Py_ssize_t start_index, end_index
    for i in range(num_layers):
        # Get layer's index information
        start_index = start_index_by_layer_ptr[i]
        num_slices  = num_slices_by_layer_ptr[i]
        end_index   = start_index + num_slices

        # Get solution and y information
        num_sols = num_solutions_by_layer_ptr[i]
        num_ys   = 2 * num_sols

        # Setup pointer array slices
        layer_radius_ptr    = &radius_array[start_index]
        layer_density_ptr   = &density_array[start_index]
        layer_gravity_ptr   = &gravity_array[start_index]
        layer_bulk_mod_ptr  = &bulk_modulus_array[start_index]
        layer_shear_mod_ptr = &complex_shear_modulus_array[start_index]

        # Setup other pointers
        layer_rtols_ptr = &rtols_by_layer_ptr[i * MAX_NUM_YS_REAL]
        layer_atols_ptr = &atols_by_layer_ptr[i * MAX_NUM_YS_REAL]

        # Get physical parameters at the top and bottom of the layer
        radius_lower  = layer_radius_ptr[0]
        density_lower = layer_density_ptr[0]
        gravity_lower = layer_gravity_ptr[0]
        bulk_lower    = layer_bulk_mod_ptr[0]
        shear_lower   = layer_shear_mod_ptr[0]
        radius_upper  = layer_radius_ptr[num_slices - 1]
        density_upper = layer_density_ptr[num_slices - 1]
        gravity_upper = layer_gravity_ptr[num_slices - 1]

        # Determine max step size (if not provided by user)
        if max_step_from_arrays:
            # Maximum step size during integration can not exceed the average radial slice size.
            max_step_touse = (radius_upper - radius_lower) / <double>num_slices

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
                frequency,
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
            num_slices,
            layer_radius_ptr,
            layer_density_ptr,
            layer_gravity_ptr,
            layer_bulk_mod_ptr,
            layer_shear_mod_ptr,
            frequency,
            degree_l,
            G_to_use,
            radial_span,
            initial_solutions_real_view[0, :],
            layer_rtols_ptr,
            layer_atols_ptr,
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
            solver.solve(reset=True)

            # Check for problems
            if not solver.success:
                # Problem with integration.
                feedback_str = f'Integration problem at layer {i}; solution {j}:\n\t{solver.message}'
                error = True
                break

            # If no problems, store results.
            # Need to make a copy because the solver pointers will be reallocated during the next solution.
            # Get storage pointer for this solution
            storage_by_y = storage_by_solution[j]
            for m in range(num_slices):
                for k in range(num_ys):

                    # Convert 2x real ys to 1x complex ys
                    storage_by_y[m * num_ys + k] = \
                        (solver.solution_y_ptr[m * num_ys + (2 * k)] +
                         1.0j * solver.solution_y_ptr[m * num_ys + (2 * k) + 1])

                    # Store top most result for initial condition for the next layer
                    if m == (num_slices - 1):
                        last_layer_top_solutions_view[j, k] = storage_by_y[m * num_ys + k]

        if error:
            # Error was encountered during integration
            break

        # Prepare for next layer
        last_layer_upper_gravity = gravity_upper
        last_layer_upper_density = density_upper


    if not error:
        feedback_str = 'Integration completed for all layers.'




    # Release pointers
    PyMem_Free(num_solutions_by_layer_ptr)
    PyMem_Free(num_slices_by_layer_ptr)
    PyMem_Free(layer_radius_ptr)
    PyMem_Free(layer_density_ptr)
    PyMem_Free(layer_gravity_ptr)
    PyMem_Free(layer_bulk_mod_ptr)
    PyMem_Free(layer_shear_mod_ptr)
    PyMem_Free(rtols_by_layer_ptr)
    PyMem_Free(atols_by_layer_ptr)

    # Deconstruct main solution pointer
    # Main storage pointer is setup like [layer_i][solution_i][y_i + r_i]
    for i in range(num_layers):
        num_sols   = num_solutions_by_layer_ptr[i]
        num_slices = num_slices_by_layer_ptr[i]
        # Number of ys = 2x num sols
        num_ys = 2 * num_sols
        for j in range(num_sols):
            PyMem_Free(main_storage[i][j])
        PyMem_Free(main_storage[i])
    PyMem_Free(main_storage)

