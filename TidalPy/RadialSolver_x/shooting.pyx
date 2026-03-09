# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Shooting method solver for the viscoelastic-gravitational radial problem.
# This must remain as a Cython cdef function because it calls CyRK's Cython-only API.

from libc.math cimport fmin, fmax, fabs
from libc.stdio cimport printf
from libc.string cimport memcpy

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.string cimport to_string
from libcpp.string cimport string as cpp_string
from libcpp.memory cimport unique_ptr, make_unique

from CyRK cimport CySolverResult, CySolveOutput, DiffeqFuncType, PreEvalFunc, Event, ODEMethod
from CyRK.cy.cysolver_api cimport baseline_cysolve_ivp_noreturn
from CyRK.utils.utils cimport allocate_mem, free_mem
from CyRK.utils.vector cimport vector

from TidalPy.constants cimport d_PI, d_EPS, d_NAN
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution
from TidalPy.RadialSolver_x.constants cimport C_MAX_NUM_Y, C_MAX_NUM_Y_REAL, C_MAX_NUM_SOL
from TidalPy.RadialSolver_x.rs_solution cimport c_RadialSolutionStorage
from TidalPy.RadialSolver_x.starting.driver cimport c_find_starting_conditions
from TidalPy.RadialSolver_x.interfaces.interfaces cimport c_solve_upper_y_at_interface
from TidalPy.RadialSolver_x.interfaces.reversed cimport c_top_to_bottom_interface_bc
from TidalPy.RadialSolver_x.derivatives.odes cimport c_RadialSolverArgs, c_find_layer_diffeq, c_find_num_shooting_solutions
from TidalPy.RadialSolver_x.boundaries.boundaries cimport c_apply_surface_bc
from TidalPy.RadialSolver_x.boundaries.surface_bc cimport c_get_surface_bc
from TidalPy.RadialSolver_x.collapse.collapse cimport c_collapse_layer_solution


cdef double d_EPS_DBL_10000 = 10000.0 * d_EPS

# Convenience complex constants
cdef double complex cmplx_NAN = d_NAN + 1j * d_NAN
cdef double complex cmplx_zero = 0.0 + 0.0j


def find_num_shooting_solutions(int layer_type, int is_static, int is_incompressible):
    return c_find_num_shooting_solutions(layer_type, is_static, is_incompressible)


cdef int cf_shooting_solver(
        c_RadialSolutionStorage* solution_storage_ptr,
        double frequency,
        double planet_bulk_density,
        int* layer_types_ptr,
        int* is_static_by_layer_ptr,
        int* is_incompressible_by_layer_ptr,
        size_t* first_slice_index_by_layer_ptr,
        size_t* num_slices_by_layer_ptr,
        size_t num_layers,
        size_t num_bc_models,
        int* bc_models_ptr,
        double G_to_use,
        int degree_l,
        cpp_bool use_kamata,
        double starting_radius,
        double start_radius_tolerance,
        ODEMethod integration_method,
        double integration_rtol,
        double integration_atol,
        cpp_bool scale_rtols_by_layer_type,
        size_t max_num_steps,
        size_t expected_size,
        size_t max_ram_MB,
        double max_step,
        cpp_bool verbose
        ) noexcept nogil:
    """ Solves the viscoelastic-gravitational problem for planets using a shooting method. """

    # Get raw pointer of radial solver storage and eos storage
    cdef c_EOSSolution* eos_solution_storage_ptr = solution_storage_ptr.get_eos_solution_ptr()

    solution_storage_ptr.message = cpp_string(b'RadialSolver_x.ShootingMethod:: Starting integration\n')
    if verbose:
        printf(solution_storage_ptr.message.c_str())

    # General indexing
    cdef double last_radius_check = d_NAN
    cdef double radius_check      = d_NAN
    cdef size_t i
    cdef size_t current_layer_i
    cdef size_t layer_i_reversed
    cdef size_t first_slice_index
    cdef size_t slice_i
    cdef size_t solution_i
    cdef size_t y_i
    cdef size_t ytype_i

    # Type conversions
    cdef double degree_l_dbl = <double>degree_l

    # Alias pointers to EOS properties
    cdef double* radius_array_ptr  = eos_solution_storage_ptr.radius_array_vec.data()
    cdef double* gravity_array_ptr = eos_solution_storage_ptr.gravity_array_vec.data()
    cdef double* density_array_ptr = eos_solution_storage_ptr.density_array_vec.data()

    # Recast shear/bulk from double to double complex
    cdef double complex* complex_shear_array_ptr = <double complex*>eos_solution_storage_ptr.complex_shear_array_vec.data()
    cdef double complex* complex_bulk_array_ptr  = <double complex*>eos_solution_storage_ptr.complex_bulk_array_vec.data()

    # Pull out key information
    cdef size_t total_slices = eos_solution_storage_ptr.radius_array_size
    cdef size_t top_slice_i  = total_slices - 1

    # Pull out constants
    cdef double planet_radius   = eos_solution_storage_ptr.radius
    cdef double surface_gravity = eos_solution_storage_ptr.surface_gravity

    # Find boundary condition at the top of the planet
    cdef size_t num_ytypes = num_bc_models

    # 15 = 5 (max_num_solutions) * 3 (number of surface conditions)
    cdef double[15] boundary_conditions
    cdef double* bc_pointer = &boundary_conditions[0]
    cdef int bc_err = c_get_surface_bc(
        bc_pointer,
        bc_models_ptr,
        num_ytypes,
        planet_radius,
        planet_bulk_density,
        degree_l_dbl
        )

    if bc_err != 0:
        solution_storage_ptr.error_code = bc_err
        solution_storage_ptr.success = False
        solution_storage_ptr.message = cpp_string(b'RadialSolver_x.ShootingMethod:: Error computing surface BCs.\n')
        return bc_err

    # Integration information
    cdef size_t num_extra              = 0
    cdef double first_step_size        = 0.0
    cdef cpp_bool capture_dense_output = False

    # Max step size
    cdef double max_step_to_use        = d_NAN
    cdef cpp_bool max_step_from_arrays = False
    if max_step == 0.0:
        max_step_from_arrays = True
    else:
        max_step_to_use = max_step

    # Setup tolerance arrays
    cdef vector[double] rtols_vec = vector[double](1)
    rtols_vec[0] = integration_rtol
    cdef vector[double] atols_vec = vector[double](1)
    atols_vec[0] = integration_atol

    # Create storage for flags and information about each layer.
    cdef vector[size_t] num_solutions_by_layer_vec = vector[size_t]()
    num_solutions_by_layer_vec.resize(num_layers)
    cdef size_t* num_solutions_by_layer_ptr = num_solutions_by_layer_vec.data()

    cdef int layer_type
    cdef int layer_is_static
    cdef int layer_is_incomp
    cdef int layer_below_type
    cdef int layer_below_is_static
    cdef int layer_below_is_incomp
    cdef size_t layer_slices
    cdef size_t num_sols
    cdef size_t num_ys
    cdef size_t num_ys_dbl
    cdef size_t layer_below_num_sols
    cdef size_t layer_below_num_ys
    cdef double layer_upper_radius = d_NAN
    cdef double layer_rtol_real    = d_NAN
    cdef double layer_rtol_imag    = d_NAN
    cdef double layer_atol_real    = d_NAN
    cdef double layer_atol_imag    = d_NAN

    for current_layer_i in range(num_layers):
        layer_type      = layer_types_ptr[current_layer_i]
        layer_is_static = is_static_by_layer_ptr[current_layer_i]
        layer_is_incomp = is_incompressible_by_layer_ptr[current_layer_i]

        num_sols = c_find_num_shooting_solutions(layer_type, layer_is_static, layer_is_incomp)
        num_solutions_by_layer_ptr[current_layer_i] = num_sols

    # Build main storage [layer][solution][y * slice]
    cdef double complex*** main_storage_ptr = <double complex***> allocate_mem(
        num_layers * sizeof(double complex**),
        'main_storage_ptr (shooting; init)'
        )

    cdef double complex** storage_by_solution_ptr = NULL
    cdef double complex* storage_by_y_ptr = NULL

    for current_layer_i in range(num_layers):
        num_sols     = num_solutions_by_layer_ptr[current_layer_i]
        layer_slices = num_slices_by_layer_ptr[current_layer_i]
        num_ys = 2 * num_sols

        storage_by_solution_ptr = <double complex**> allocate_mem(
            num_sols * sizeof(double complex*),
            'storage_by_solution_ptr (shooting; init)'
            )

        for solution_i in range(num_sols):
            storage_by_y_ptr = <double complex*> allocate_mem(
                layer_slices * num_ys * sizeof(double complex),
                'storage_by_y_ptr (shooting; init)'
                )
            storage_by_solution_ptr[solution_i] = storage_by_y_ptr
            storage_by_y_ptr = NULL
        main_storage_ptr[current_layer_i] = storage_by_solution_ptr
        storage_by_solution_ptr = NULL

    # Storage for uppermost ys for each solution (max: 6 ys * 3 sols = 18)
    cdef double complex[18] uppermost_y_per_solution
    cdef double complex* uppermost_y_per_solution_ptr = &uppermost_y_per_solution[0]
    for i in range(18):
        uppermost_y_per_solution_ptr[i] = cmplx_NAN

    # Layer specific pointers
    cdef vector[double] layer_radius_vec = vector[double](eos_solution_storage_ptr.radius_array_vec.size())

    cdef double* layer_radius_ptr
    cdef double* layer_density_ptr
    cdef double* layer_gravity_ptr
    cdef double complex* layer_bulk_mod_ptr
    cdef double complex* layer_shear_mod_ptr

    cdef double radius_lower  = d_NAN
    cdef double radius_upper  = d_NAN
    cdef double density_lower = d_NAN
    cdef double density_upper = d_NAN
    cdef double gravity_lower = d_NAN
    cdef double gravity_upper = d_NAN
    cdef double complex bulk_lower  = cmplx_NAN
    cdef double complex bulk_upper  = cmplx_NAN
    cdef double complex shear_upper = cmplx_NAN
    cdef double complex shear_lower = cmplx_NAN

    # Properties at interfaces
    cdef double static_liquid_density    = d_NAN
    cdef double interface_gravity        = d_NAN
    cdef double last_layer_upper_gravity = d_NAN
    cdef double last_layer_upper_density = d_NAN
    cdef double last_layer_upper_radius  = d_NAN

    # Starting solutions (max: 6 ys * 3 sols = 18 complex, 36 real)
    cdef double complex[18] initial_y
    cdef double complex* initial_y_ptr = &initial_y[0]
    cdef double[36] initial_y_only_real
    cdef double* initial_y_only_real_ptr = &initial_y_only_real[0]
    cdef cpp_bool starting_y_check = False

    cdef double complex dcomplex_tmp = cmplx_NAN

    # Final solution
    cdef size_t num_output_ys = C_MAX_NUM_Y * num_ytypes
    cdef double complex* solution_ptr = <double complex*>solution_storage_ptr.full_solution_vec.data()

    # Collapse variables
    cdef double layer_above_lower_gravity
    cdef double layer_above_lower_density
    cdef double liquid_density_at_interface
    cdef int layer_above_type
    cdef int layer_above_is_static
    cdef int layer_above_is_incomp

    # Integration ODE setup
    cdef double[9] eos_interp_array
    cdef double* eos_interp_array_ptr = &eos_interp_array[0]
    cdef unique_ptr[CySolverResult] integration_solution_uptr = make_unique[CySolverResult](integration_method)
    cdef CySolverResult* integration_solution_ptr = NULL
    cdef double* integrator_data_ptr = NULL

    # Diffeq args
    cdef vector[char] diffeq_args_vec = vector[char](sizeof(c_RadialSolverArgs))
    cdef c_RadialSolverArgs* diffeq_args_ptr = <c_RadialSolverArgs*>diffeq_args_vec.data()
    cdef PreEvalFunc diffeq_preeval_ptr
    cdef DiffeqFuncType layer_diffeq = NULL
    cdef vector[double] y0_vec = vector[double](C_MAX_NUM_Y_REAL)

    # Set diffeq inputs that do not change with layer
    diffeq_args_ptr.degree_l         = degree_l_dbl
    diffeq_args_ptr.lp1              = degree_l_dbl + 1.0
    diffeq_args_ptr.lm1              = degree_l_dbl - 1.0
    diffeq_args_ptr.llp1             = degree_l_dbl * (degree_l_dbl + 1.0)
    diffeq_args_ptr.G                = G_to_use
    diffeq_args_ptr.grav_coeff       = 4.0 * d_PI * G_to_use
    diffeq_args_ptr.frequency        = frequency
    diffeq_args_ptr.layer_index      = 0
    diffeq_args_ptr.eos_solution_ptr = eos_solution_storage_ptr

    # Constant vectors
    cdef double complex[3] constant_vector
    cdef double complex* constant_vector_ptr = &constant_vector[0]
    cdef double complex[3] layer_above_constant_vector
    cdef double complex* layer_above_constant_vector_ptr = &layer_above_constant_vector[0]
    cdef double complex[6] surface_solutions
    cdef double complex* surface_solutions_ptr = &surface_solutions[0]

    for i in range(6):
        if i < 3:
            constant_vector_ptr[i] = cmplx_NAN
            layer_above_constant_vector_ptr[i] = cmplx_NAN
        surface_solutions_ptr[i] = cmplx_NAN

    cdef int bc_solution_info = -999

    # Events (empty)
    cdef vector[Event] events_vec = vector[Event]()

    # Determine starting radius
    if starting_radius == 0.0:
        starting_radius = planet_radius * start_radius_tolerance**(1.0 / degree_l_dbl)
        starting_radius = fmin(starting_radius, 0.95 * planet_radius)

    cdef double starting_gravity       = d_NAN
    cdef double starting_density       = d_NAN
    cdef double complex starting_shear = d_NAN
    cdef double complex starting_bulk  = d_NAN

    # Find starting layer
    cdef size_t start_layer_i           = 0
    cdef size_t last_index_before_start = 0
    cdef size_t start_index_in_layer    = 0

    last_radius_check = 0.0
    for current_layer_i in range(num_layers):
        layer_upper_radius = eos_solution_storage_ptr.upper_radius_bylayer_vec[current_layer_i]
        if current_layer_i == 0:
            last_layer_upper_radius = 0.0
        else:
            last_layer_upper_radius = eos_solution_storage_ptr.upper_radius_bylayer_vec[current_layer_i - 1]

        if last_layer_upper_radius < starting_radius <= layer_upper_radius:
            start_layer_i = current_layer_i
            first_slice_index = first_slice_index_by_layer_ptr[current_layer_i]

            start_index_in_layer = 0
            for slice_i in range(first_slice_index, first_slice_index + num_slices_by_layer_ptr[current_layer_i]):
                radius_check = radius_array_ptr[slice_i]
                if last_radius_check < starting_radius <= radius_check:
                    if slice_i == 0:
                        last_index_before_start = 0
                    else:
                        last_index_before_start = slice_i - 1
                    break
                else:
                    start_index_in_layer += 1
                    last_radius_check = radius_check
            break
        else:
            last_radius_check = last_layer_upper_radius

    # NAN out data below the starting radius
    for slice_i in range(last_index_before_start + 1):
        for ytype_i in range(num_ytypes):
            for y_i in range(C_MAX_NUM_Y):
                solution_ptr[slice_i * C_MAX_NUM_Y * num_ytypes + ytype_i * C_MAX_NUM_Y + y_i] = cmplx_NAN

    # =================== Main Integration Loop ===================
    for current_layer_i in range(start_layer_i, num_layers):
        layer_slices = num_slices_by_layer_ptr[current_layer_i]
        if current_layer_i == start_layer_i:
            first_slice_index = last_index_before_start + 1
            layer_slices -= start_index_in_layer
        else:
            first_slice_index = first_slice_index_by_layer_ptr[current_layer_i]

        num_sols   = num_solutions_by_layer_ptr[current_layer_i]
        num_ys     = 2 * num_sols
        num_ys_dbl = 2 * num_ys

        # Setup pointer slices for this layer
        layer_radius_ptr    = &radius_array_ptr[first_slice_index]
        layer_density_ptr   = &density_array_ptr[first_slice_index]
        layer_gravity_ptr   = &gravity_array_ptr[first_slice_index]
        layer_shear_mod_ptr = &complex_shear_array_ptr[first_slice_index]
        layer_bulk_mod_ptr  = &complex_bulk_array_ptr[first_slice_index]

        layer_radius_vec.resize(layer_slices)
        memcpy(layer_radius_vec.data(), layer_radius_ptr, sizeof(double) * layer_slices)

        # Get physical parameters at top and bottom
        if current_layer_i == start_layer_i:
            eos_solution_storage_ptr.call(current_layer_i, starting_radius, eos_interp_array_ptr)
            starting_gravity = eos_interp_array_ptr[0]
            if starting_gravity < d_EPS:
                starting_gravity = d_EPS
            starting_density = eos_interp_array_ptr[4]
            starting_shear   = <double complex>eos_interp_array_ptr[5] + 1j * <double complex>eos_interp_array_ptr[6]
            starting_bulk    = <double complex>eos_interp_array_ptr[7] + 1j * <double complex>eos_interp_array_ptr[8]

            radius_lower  = starting_radius
            gravity_lower = starting_gravity
            density_lower = starting_density
            shear_lower   = starting_shear
            bulk_lower    = starting_bulk
        else:
            radius_lower  = layer_radius_ptr[0]
            gravity_lower = layer_gravity_ptr[0]
            density_lower = layer_density_ptr[0]
            shear_lower   = layer_shear_mod_ptr[0]
            bulk_lower    = layer_bulk_mod_ptr[0]

        radius_upper  = layer_radius_ptr[layer_slices - 1]
        gravity_upper = layer_gravity_ptr[layer_slices - 1]
        density_upper = layer_density_ptr[layer_slices - 1]
        shear_upper   = layer_shear_mod_ptr[layer_slices - 1]
        bulk_upper    = layer_bulk_mod_ptr[layer_slices - 1]

        # Determine max step size
        if max_step_from_arrays:
            max_step_to_use = fabs(0.33 * (radius_upper - radius_lower))
            max_step_to_use = fmax(max_step_to_use, d_EPS_DBL_10000)

        layer_type      = layer_types_ptr[current_layer_i]
        layer_is_static = is_static_by_layer_ptr[current_layer_i]
        layer_is_incomp = is_incompressible_by_layer_ptr[current_layer_i]

        # Scale rtols by layer type
        if scale_rtols_by_layer_type:
            rtols_vec.resize(num_ys * 2)
            atols_vec.resize(num_ys * 2)

            for y_i in range(num_ys):
                layer_rtol_real = integration_rtol
                layer_rtol_imag = integration_rtol
                layer_atol_real = integration_atol
                layer_atol_imag = integration_atol

                if layer_type == 0:
                    if (y_i == 1) or (y_i == 2):
                        layer_rtol_real *= 0.1
                        layer_rtol_imag *= 0.1
                else:
                    if not layer_is_static:
                        if y_i == 1:
                            layer_rtol_real *= 0.01
                            layer_rtol_imag *= 0.01
                rtols_vec[2 * y_i]     = layer_rtol_real
                rtols_vec[2 * y_i + 1] = layer_rtol_imag
                atols_vec[2 * y_i]     = layer_atol_real
                atols_vec[2 * y_i + 1] = layer_atol_imag

        # Initialize initial_y to NaN
        for y_i in range(36):
            if y_i < 18:
                initial_y_ptr[y_i] = cmplx_NAN
            initial_y_only_real_ptr[y_i] = d_NAN

        if current_layer_i == start_layer_i:
            # First layer: use starting conditions
            c_find_starting_conditions(
                &solution_storage_ptr.success,
                solution_storage_ptr.message,
                layer_type,
                layer_is_static,
                layer_is_incomp,
                use_kamata,
                frequency,
                radius_lower,
                density_lower,
                <cpp_complex[double]>bulk_lower,
                <cpp_complex[double]>shear_lower,
                degree_l,
                G_to_use,
                C_MAX_NUM_Y,
                <cpp_complex[double]*>initial_y_ptr,
                starting_y_check
                )

            if not solution_storage_ptr.success:
                solution_storage_ptr.error_code = -10
                break
        else:
            # Not first layer: use interface function
            layer_below_type      = layer_types_ptr[current_layer_i - 1]
            layer_below_is_static = is_static_by_layer_ptr[current_layer_i - 1]
            layer_below_is_incomp = is_incompressible_by_layer_ptr[current_layer_i - 1]

            interface_gravity = 0.5 * (gravity_lower + last_layer_upper_gravity)

            # Find the density needed for some initial conditions
            if (layer_type == 0) and (layer_below_type == 0):
                static_liquid_density = d_NAN
            elif not (layer_type == 0) and (layer_below_type == 0):
                static_liquid_density = density_lower
            elif (layer_type == 0) and not (layer_below_type == 0):
                static_liquid_density = last_layer_upper_density
            else:
                if layer_is_static and layer_below_is_static:
                    static_liquid_density = density_lower
                elif layer_is_static and not layer_below_is_static:
                    static_liquid_density = density_lower
                elif not layer_is_static and layer_below_is_static:
                    static_liquid_density = last_layer_upper_density
                else:
                    static_liquid_density = d_NAN

            c_solve_upper_y_at_interface(
                <cpp_complex[double]*>uppermost_y_per_solution_ptr,
                <cpp_complex[double]*>initial_y_ptr,
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
                )

        # Reset uppermost y
        for i in range(18):
            uppermost_y_per_solution_ptr[i] = cmplx_NAN

        # Convert initial conditions from complex to 2x real for integration
        for solution_i in range(num_sols):
            for y_i in range(num_ys):
                dcomplex_tmp = initial_y_ptr[solution_i * C_MAX_NUM_Y + y_i]
                initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL + 2 * y_i]     = dcomplex_tmp.real
                initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL + 2 * y_i + 1] = dcomplex_tmp.imag

        # Find correct diffeq
        layer_diffeq = c_find_layer_diffeq(layer_type, layer_is_static, layer_is_incomp)

        # Set the layer index
        diffeq_args_ptr.layer_index = current_layer_i

        # Get storage pointer for this layer
        storage_by_solution_ptr = main_storage_ptr[current_layer_i]

        # Solve for each solution
        for solution_i in range(num_sols):

            # Ensure integration solution is initialized
            if not integration_solution_uptr:
                integration_solution_uptr = make_unique[CySolverResult](integration_method)
            integration_solution_ptr = integration_solution_uptr.get()

            y0_vec.resize(num_ys_dbl)
            memcpy(y0_vec.data(), &initial_y_only_real_ptr[solution_i * C_MAX_NUM_Y_REAL], sizeof(double) * num_ys_dbl)

            ###### Integrate! #######
            baseline_cysolve_ivp_noreturn(
                integration_solution_ptr,
                layer_diffeq,
                radius_lower,
                radius_upper,
                y0_vec,
                expected_size,
                num_extra,
                diffeq_args_vec,
                max_num_steps,
                max_ram_MB,
                capture_dense_output,
                layer_radius_vec,
                diffeq_preeval_ptr,
                events_vec,
                rtols_vec,
                atols_vec,
                max_step_to_use,
                first_step_size,
                True
                )
            #########################

            # Store diagnostic data
            solution_storage_ptr.shooting_method_steps_taken_vec[(3 * current_layer_i) + solution_i] = \
                integration_solution_ptr.steps_taken

            # Check for problems
            if not integration_solution_ptr.success:
                solution_storage_ptr.error_code = -11
                solution_storage_ptr.success = False
                solution_storage_ptr.message = cpp_string(b'RadialSolver_x.ShootingMethod:: Integration problem at layer ') + to_string(current_layer_i) + cpp_string(b'; solution ') + to_string(solution_i) + cpp_string(b':\n\t') + integration_solution_ptr.message + cpp_string(b'\n')
                if verbose:
                    printf(solution_storage_ptr.message.c_str())
                return solution_storage_ptr.error_code

            # Store results
            integrator_data_ptr = &integration_solution_ptr.solution[0]
            storage_by_y_ptr = storage_by_solution_ptr[solution_i]

            for slice_i in range(layer_slices):
                for y_i in range(num_ys):
                    # Convert 2x real ys to 1x complex ys
                    storage_by_y_ptr[num_ys * slice_i + y_i] = (
                        <double complex>integrator_data_ptr[num_ys_dbl * slice_i + (2 * y_i)] +
                        1j * <double complex>integrator_data_ptr[num_ys_dbl * slice_i + (2 * y_i) + 1]
                        )

                    if slice_i == (layer_slices - 1):
                        uppermost_y_per_solution_ptr[solution_i * C_MAX_NUM_Y + y_i] = storage_by_y_ptr[num_ys * slice_i + y_i]

        if solution_storage_ptr.error_code != 0:
            solution_storage_ptr.success = False
            return solution_storage_ptr.error_code

        # Prepare for next layer
        layer_below_num_sols     = num_sols
        layer_below_num_ys       = num_ys
        last_layer_upper_gravity = gravity_upper
        last_layer_upper_density = density_upper

    # =================== Collapse Phase ===================
    if solution_storage_ptr.error_code < 0 or not solution_storage_ptr.success:
        solution_storage_ptr.success = False
        if integration_solution_ptr:
            solution_storage_ptr.message = cpp_string(b'RadialSolver_x.ShootingMethod:: Integration failed:\n\t') + integration_solution_ptr.message + cpp_string(b'\n')
        if verbose:
            printf(solution_storage_ptr.message.c_str())
        return solution_storage_ptr.error_code

    else:
        solution_storage_ptr.message = cpp_string(b'Integration completed for all layers. Beginning solution collapse.\n')

        for ytype_i in range(num_ytypes):
            solution_storage_ptr.message = cpp_string(b'Collapsing radial solutions for "') + to_string(ytype_i) + cpp_string(b'" solver.\n')

            # Reset variables for this solver
            bc_solution_info = -999
            constant_vector_ptr[0] = cmplx_NAN
            constant_vector_ptr[1] = cmplx_NAN
            constant_vector_ptr[2] = cmplx_NAN
            layer_above_lower_gravity   = d_NAN
            layer_above_lower_density   = d_NAN
            liquid_density_at_interface = d_NAN
            layer_above_type            = 9
            layer_above_is_static       = 0
            layer_above_is_incomp       = 0

            if verbose:
                printf(solution_storage_ptr.message.c_str())

            # Work from surface down to starting layer
            for current_layer_i in range(num_layers - start_layer_i):
                layer_i_reversed = num_layers - (current_layer_i + 1)

                layer_slices = num_slices_by_layer_ptr[layer_i_reversed]
                if layer_i_reversed == start_layer_i:
                    first_slice_index = last_index_before_start + 1
                    layer_slices -= start_index_in_layer
                else:
                    first_slice_index = first_slice_index_by_layer_ptr[layer_i_reversed]

                num_sols = num_solutions_by_layer_ptr[layer_i_reversed]
                num_ys   = 2 * num_sols

                layer_radius_ptr    = &radius_array_ptr[first_slice_index]
                layer_density_ptr   = &density_array_ptr[first_slice_index]
                layer_gravity_ptr   = &gravity_array_ptr[first_slice_index]
                layer_bulk_mod_ptr  = &complex_bulk_array_ptr[first_slice_index]
                layer_shear_mod_ptr = &complex_shear_array_ptr[first_slice_index]

                if layer_i_reversed == start_layer_i:
                    radius_lower  = starting_radius
                    gravity_lower = starting_gravity
                    density_lower = starting_density
                    shear_lower   = starting_shear
                    bulk_lower    = starting_bulk
                else:
                    radius_lower  = layer_radius_ptr[0]
                    density_lower = layer_density_ptr[0]
                    gravity_lower = layer_gravity_ptr[0]
                    bulk_lower    = layer_bulk_mod_ptr[0]
                    shear_lower   = layer_shear_mod_ptr[0]

                radius_upper  = layer_radius_ptr[layer_slices - 1]
                density_upper = layer_density_ptr[layer_slices - 1]
                gravity_upper = layer_gravity_ptr[layer_slices - 1]

                layer_type      = layer_types_ptr[layer_i_reversed]
                layer_is_static = is_static_by_layer_ptr[layer_i_reversed]
                layer_is_incomp = is_incompressible_by_layer_ptr[layer_i_reversed]

                storage_by_solution_ptr = main_storage_ptr[layer_i_reversed]

                # Get radial solution values at the top of the layer
                for solution_i in range(num_sols):
                    for y_i in range(num_ys):
                        uppermost_y_per_solution_ptr[solution_i * C_MAX_NUM_Y + y_i] = \
                            storage_by_solution_ptr[solution_i][(layer_slices - 1) * num_ys + y_i]

                if current_layer_i == 0:
                    # Surface layer: Apply surface boundary conditions
                    c_apply_surface_bc(
                        <cpp_complex[double]*>constant_vector_ptr,
                        &bc_solution_info,
                        bc_pointer,
                        <cpp_complex[double]*>uppermost_y_per_solution_ptr,
                        surface_gravity,
                        G_to_use,
                        num_sols,
                        C_MAX_NUM_Y,
                        ytype_i,
                        layer_type,
                        layer_is_static,
                        layer_is_incomp
                        )

                    if bc_solution_info != 0:
                        solution_storage_ptr.error_code = -12
                        solution_storage_ptr.success = False
                        solution_storage_ptr.message = cpp_string(b'RadialSolver_x.ShootingMethod:: Error in surface BC. ZGESV code: ') + to_string(bc_solution_info) + cpp_string(b'\n')
                        if verbose:
                            printf(solution_storage_ptr.message.c_str())
                        return solution_storage_ptr.error_code
                else:
                    # Interior layers: find constants from layer above
                    c_top_to_bottom_interface_bc(
                        <cpp_complex[double]*>constant_vector_ptr,
                        <cpp_complex[double]*>layer_above_constant_vector_ptr,
                        <cpp_complex[double]*>uppermost_y_per_solution_ptr,
                        gravity_upper, layer_above_lower_gravity,
                        density_upper, layer_above_lower_density,
                        layer_type, layer_above_type,
                        layer_is_static, layer_above_is_static,
                        layer_is_incomp, layer_above_is_incomp,
                        num_sols, C_MAX_NUM_Y
                        )

                # Collapse solutions
                c_collapse_layer_solution(
                    <cpp_complex[double]*>solution_ptr,
                    <cpp_complex[double]*>constant_vector_ptr,
                    <cpp_complex[double]**>storage_by_solution_ptr,
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
                    )

                # Setup for next layer
                layer_above_lower_gravity = gravity_lower
                layer_above_lower_density = density_lower
                layer_above_type          = layer_type
                layer_above_is_static     = layer_is_static
                layer_above_is_incomp     = layer_is_incomp

                if num_sols == 1:
                    layer_above_constant_vector_ptr[0] = constant_vector_ptr[0]
                    layer_above_constant_vector_ptr[1] = cmplx_NAN
                    layer_above_constant_vector_ptr[2] = cmplx_NAN
                elif num_sols == 2:
                    layer_above_constant_vector_ptr[0] = constant_vector_ptr[0]
                    layer_above_constant_vector_ptr[1] = constant_vector_ptr[1]
                    layer_above_constant_vector_ptr[2] = cmplx_NAN
                elif num_sols == 3:
                    layer_above_constant_vector_ptr[0] = constant_vector_ptr[0]
                    layer_above_constant_vector_ptr[1] = constant_vector_ptr[1]
                    layer_above_constant_vector_ptr[2] = constant_vector_ptr[2]

    # Free memory
    if not (main_storage_ptr is NULL):
        storage_by_solution_ptr = NULL
        storage_by_y_ptr = NULL
        for current_layer_i in range(num_layers):
            num_sols = num_solutions_by_layer_ptr[current_layer_i]
            if not (main_storage_ptr[current_layer_i] is NULL):
                for solution_i in range(num_sols):
                    if not (main_storage_ptr[current_layer_i][solution_i] is NULL):
                        free_mem(main_storage_ptr[current_layer_i][solution_i])
                        main_storage_ptr[current_layer_i][solution_i] = NULL
                free_mem(main_storage_ptr[current_layer_i])
                main_storage_ptr[current_layer_i] = NULL
        free_mem(main_storage_ptr)
        main_storage_ptr = NULL

    # Update solution status and return
    if solution_storage_ptr.error_code != 0:
        solution_storage_ptr.success = False
    else:
        solution_storage_ptr.success = True
        solution_storage_ptr.message = cpp_string(b'RadialSolver_x.ShootingMethod: Completed without any noted issues.')

    return solution_storage_ptr.error_code
