# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Import cythonized functions
from libc.math cimport fmin
from libc.stdio cimport printf, sprintf
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.string cimport strcpy

from CyRK cimport cysolve_ivp, CySolverResult, CySolveOutput, DiffeqFuncType, PreEvalFunc
from CyRK.utils.utils cimport allocate_mem, free_mem
from CyRK.utils.vector cimport vector

from TidalPy.utilities.math.complex cimport cmplx_NAN, cf_build_dblcmplx
from TidalPy.constants cimport d_PI_DBL, d_EPS_DBL, d_NAN_DBL

from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC
from TidalPy.RadialSolver.constants cimport MAX_NUM_Y, MAX_NUM_Y_REAL, MAX_NUM_SOL
from TidalPy.RadialSolver.starting.driver cimport cf_find_starting_conditions
from TidalPy.RadialSolver.interfaces.interfaces cimport cf_solve_upper_y_at_interface
from TidalPy.RadialSolver.interfaces.reversed cimport cf_top_to_bottom_interface_bc
from TidalPy.RadialSolver.derivatives.odes cimport RadialSolverDiffeqArgStruct, cf_find_layer_diffeq
from TidalPy.RadialSolver.boundaries.boundaries cimport cf_apply_surface_bc
from TidalPy.RadialSolver.boundaries.surface_bc cimport cf_get_surface_bc
from TidalPy.RadialSolver.collapse.collapse cimport cf_collapse_layer_solution


cdef size_t cf_find_num_shooting_solutions(
        int layer_type,
        bint is_static,
        bint is_incompressible
        ) noexcept nogil:
    """ Determine number of solutions required for layer based on assumptions.
    
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

    """

    # Initialize
    cdef Py_ssize_t num_sols
    num_sols = 0

    if (layer_type == 0):
        # Solid
        if is_static:
            if is_incompressible:
                # TODO: Confirm
                num_sols = 3
            else:
                num_sols = 3
        else:
            # Dynamic
            if is_incompressible:
                # TODO: Confirm
                num_sols = 3
            else:
                num_sols = 3
    else:
        # Liquid
        if is_static:
            if is_incompressible:
                # TODO: Confirm
                num_sols = 1
            else:
                num_sols = 1
        else:
            # Dynamic
            if is_incompressible:
                # TODO: Confirm
                num_sols = 2
            else:
                num_sols = 2
    return num_sols


def find_num_shooting_solutions(
        int layer_type,
        bint is_static,
        bint is_incompressible
        ):
    
    return cf_find_num_shooting_solutions(layer_type, is_static, is_incompressible)


cdef int cf_shooting_solver(
        RadialSolutionStorageCC* solution_storage_ptr,
        double frequency,
        double planet_bulk_density,
        int* layer_types_ptr,
        bint* is_static_by_layer_ptr,
        bint* is_incompressible_by_layer_ptr,
        vector[size_t] first_slice_index_by_layer_vec,
        vector[size_t] num_slices_by_layer_vec,
        size_t num_bc_models,
        int* bc_models_ptr,
        double G_to_use,
        int degree_l,
        cpp_bool use_kamata,
        double starting_radius,
        double start_radius_tolerance,
        int integration_method,
        double integration_rtol,
        double integration_atol,
        cpp_bool scale_rtols_by_layer_type,
        size_t max_num_steps,
        size_t expected_size,
        size_t max_ram_MB,
        double max_step,
        cpp_bool verbose
        ) noexcept nogil:
    """ Solves the viscoelastic-gravitational problem for planets using a shooting method.
    """

    # Feedback
    cdef char[256] message
    cdef char* message_ptr = &message[0]

    # Get raw pointer of radial solver storage and eos storage
    cdef EOSSolutionCC* eos_solution_storage_ptr = solution_storage_ptr.get_eos_solution_ptr()

    solution_storage_ptr.set_message('RadialSolver.ShootingMethod:: Starting integration\n')
    if verbose:
        printf(solution_storage_ptr.message_ptr)

    # General indexing
    cdef double last_radius_check = d_NAN_DBL
    cdef double radius_check      = d_NAN_DBL
    cdef size_t i
    cdef size_t current_layer_i
    cdef size_t layer_i_reversed
    cdef size_t first_slice_index
    cdef size_t slice_i
    # Indexing for: Solution | ys | ytypes
    cdef size_t solution_i
    cdef size_t y_i
    cdef size_t ytype_i

    # Type conversions
    cdef double degree_l_dbl = <double>degree_l

    # Alias pointers to EOS properties
    cdef double* radius_array_ptr  = eos_solution_storage_ptr.radius_array_vec.data()
    cdef double* gravity_array_ptr = eos_solution_storage_ptr.gravity_array_vec.data()
    # cdef double* pressure_array_ptr = eos_solution_storage_ptr.pressure_array_vec   # Unused
    cdef double* density_array_ptr = eos_solution_storage_ptr.density_array_vec.data()

    # Need to recast the storage's shear/bulk double arrays to double complex for local use
    cdef double complex* complex_shear_array_ptr = <double complex*>eos_solution_storage_ptr.complex_shear_array_vec.data()
    cdef double complex* complex_bulk_array_ptr  = <double complex*>eos_solution_storage_ptr.complex_bulk_array_vec.data()

    # Pull out key information
    cdef size_t num_layers   = eos_solution_storage_ptr.num_layers
    cdef size_t total_slices = eos_solution_storage_ptr.radius_array_size
    cdef size_t top_slice_i  = total_slices - 1
    
    # Pull out any constants now that arrays have had dimensional protocol applied to them.
    cdef double planet_radius   = eos_solution_storage_ptr.radius
    cdef double surface_gravity = eos_solution_storage_ptr.surface_gravity

    # Find boundary condition at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    # The [5] represents the maximum number of solvers that can be invoked with a single call to radial_solver
    cdef size_t num_ytypes = num_bc_models

    # 15 = 5 (max_num_solutions) * 3 (number of surface conditions)
    cdef double[15] boundary_conditions
    cdef double* bc_pointer = &boundary_conditions[0]
    cf_get_surface_bc(
        bc_pointer,  # Changed parameter
        bc_models_ptr,
        num_ytypes,
        planet_radius,
        planet_bulk_density,
        degree_l_dbl
        )

    # Integration information
    # Max step size
    cdef double max_step_to_use        = d_NAN_DBL
    cdef cpp_bool max_step_from_arrays = False
    if max_step == 0.0:
        # If max_step is zero use the array information to determine max_step_size
        max_step_from_arrays = True
    else:
        # Otherwise use user input.
        max_step_to_use = max_step

    # Setup tolerance arrays
    # For simplicity just make these all as large as the maximum number of ys.
    # Maximum number of ys = 6. Then 2x for conversion from complex to real
    cdef double[12] rtols_array
    cdef double[12] atols_array
    cdef double* rtols_ptr = &rtols_array[0]
    cdef double* atols_ptr = &atols_array[0]

    for i in range(12):
        rtols_ptr[i] = d_NAN_DBL
        atols_ptr[i] = d_NAN_DBL

    # Create storage for flags and information about each layer.
    cdef vector[size_t] num_solutions_by_layer_vec = vector[size_t]()
    num_solutions_by_layer_vec.resize(num_layers)
    cdef size_t* num_solutions_by_layer_ptr = num_solutions_by_layer_vec.data()

    # Opt: The bools above could be stores in a single char variable (per layer).
    #  Eg., 0x00 All false, 0x01 is solid, 0x10 is static and liquid, 0x11 is static and solid, etc.

    cdef int layer_type
    cdef bint layer_is_static
    cdef bint layer_is_incomp
    cdef int layer_below_type
    cdef bint layer_below_is_static
    cdef bint layer_below_is_incomp
    cdef size_t layer_slices
    cdef size_t num_sols
    cdef size_t num_ys
    cdef size_t num_ys_dbl
    cdef size_t layer_below_num_sols
    cdef size_t layer_below_num_ys
    cdef double layer_upper_radius = d_NAN_DBL
    cdef double layer_rtol_real    = d_NAN_DBL
    cdef double layer_rtol_imag    = d_NAN_DBL
    cdef double layer_atol_real    = d_NAN_DBL
    cdef double layer_atol_imag    = d_NAN_DBL

    for current_layer_i in range(num_layers):
        # Pull out information on this layer
        layer_type         = layer_types_ptr[current_layer_i]
        layer_is_static    = is_static_by_layer_ptr[current_layer_i]
        layer_is_incomp    = is_incompressible_by_layer_ptr[current_layer_i]
        layer_upper_radius = eos_solution_storage_ptr.upper_radius_bylayer_vec[current_layer_i]

        # Find number of solutions based on this layer's assumptions
        num_sols = cf_find_num_shooting_solutions(
            layer_type,
            layer_is_static,
            layer_is_incomp
            )
        num_ys = 2 * num_sols
        num_solutions_by_layer_ptr[current_layer_i] = num_sols

    # We have all the size information needed to build storage pointers
    # Main storage pointer is setup like [current_layer_i][solution_i][y_i + r_i]
    cdef double complex*** main_storage_ptr = <double complex ***> allocate_mem(
        num_layers * sizeof(double complex**),
        'main_storage_ptr (radial_solver; init)'
        )

    cdef double complex** storage_by_solution_ptr = NULL
    cdef double complex* storage_by_y_ptr = NULL

    for current_layer_i in range(num_layers):
        num_sols     = num_solutions_by_layer_ptr[current_layer_i]
        layer_slices = num_slices_by_layer_vec[current_layer_i]
        # Number of ys = 2x num sols
        num_ys = 2 * num_sols

        storage_by_solution_ptr = <double complex **> allocate_mem(
            num_sols * sizeof(double complex*), 
            'storage_by_solution_ptr (radial_solver; init)'
            ) 

        for solution_i in range(num_sols):
            storage_by_y_ptr = <double complex *> allocate_mem(
                layer_slices * num_ys * sizeof(double complex),
                'storage_by_y_ptr (radial_solver; init)'
                )

            storage_by_solution_ptr[solution_i] = storage_by_y_ptr
            storage_by_y_ptr = NULL
        main_storage_ptr[current_layer_i] = storage_by_solution_ptr
        storage_by_solution_ptr = NULL

    # Create storage for uppermost ys for each solution. We don't know how many solutions or ys per layer so assume the
    #  worst.
    cdef double complex[18] uppermost_y_per_solution
    cdef double complex* uppermost_y_per_solution_ptr = &uppermost_y_per_solution[0]

    for i in range(18):
        uppermost_y_per_solution_ptr[i] = cmplx_NAN

    # Layer specific pointers; set the size based on the layer with the most slices.
    cdef double* layer_radius_ptr
    cdef double* layer_density_ptr
    cdef double* layer_gravity_ptr
    cdef double complex* layer_bulk_mod_ptr
    cdef double complex* layer_shear_mod_ptr

    # Properties at top and bottom of layer
    cdef double[2] radial_span
    cdef double* radial_span_ptr = &radial_span[0]
    radial_span_ptr[0] = d_NAN_DBL
    radial_span_ptr[1] = d_NAN_DBL

    cdef double radius_lower  = d_NAN_DBL
    cdef double radius_upper  = d_NAN_DBL
    cdef double density_lower = d_NAN_DBL
    cdef double density_upper = d_NAN_DBL
    cdef double gravity_lower = d_NAN_DBL
    cdef double gravity_upper = d_NAN_DBL
    cdef double complex bulk_lower  = cmplx_NAN
    cdef double complex bulk_upper  = cmplx_NAN
    cdef double complex shear_upper = cmplx_NAN
    cdef double complex shear_lower = cmplx_NAN

    # Properties at interfaces between layers
    cdef double static_liquid_density    = d_NAN_DBL
    cdef double interface_gravity        = d_NAN_DBL
    cdef double last_layer_upper_gravity = d_NAN_DBL
    cdef double last_layer_upper_density = d_NAN_DBL
    cdef double last_layer_upper_radius  = d_NAN_DBL

    # Starting solutions (initial conditions / lower boundary conditions)
    # Allocate memory for the initial value arrays now. We don't know the number of solutions or ys. But the max number
    #  is not all that different from the min. So it is more efficient to assume the largest size.
    # Largest size = 6 (ys) x 3 (sols) = 18
    # For the "only_real" this is further multiplied by 2 since we convert the number of _complex_ ys to 2x _real_ ys
    cdef double complex[18] initial_y
    cdef double complex* initial_y_ptr = &initial_y[0]
    cdef double[36] initial_y_only_real
    cdef double* initial_y_only_real_ptr = &initial_y_only_real[0]
    cdef cpp_bool starting_y_check = False

    # Intermediate values
    cdef double complex dcomplex_tmp = cmplx_NAN

    # No matter the results of the integration, we know the shape and size of the final solution.
    # The number of rows will depend on if the user wants to simultaneously calculate loading Love numbers.
    cdef size_t num_output_ys = MAX_NUM_Y * num_ytypes

    # Get a reference pointer to solution array
    # Cast the solution pointer from double to double complex
    cdef double complex* solution_ptr = <double complex*>solution_storage_ptr.full_solution_vec.data()

    # During collapse, variables for the layer above the target one are used. Declare these and preset them.
    cdef double layer_above_lower_gravity
    cdef double layer_above_lower_density
    cdef double liquid_density_at_interface
    cdef int layer_above_type
    cdef bint layer_above_is_static
    cdef bint layer_above_is_incomp

    # Layer's differential equation will vary by layer type
    cdef double[9] eos_interp_array  # The "9" here is the number of variables (dependent + extra) found in the EOS solver. See TidalPy.Material.eos.solver.pyx for details.
    cdef double* eos_interp_array_ptr = &eos_interp_array[0]
    cdef CySolveOutput integration_solution
    cdef CySolverResult* integration_solution_ptr = NULL
    cdef double* integrator_data_ptr = NULL

    # The radial solver diffeq's require additional inputs other than "y" and "r". Build a structure that stores these
    # extra arguments, a point (cast to char*) will be passed to the diffeq's while solving the ODE.
    cdef RadialSolverDiffeqArgStruct diffeq_args
    cdef size_t sizeof_args                           = sizeof(RadialSolverDiffeqArgStruct)
    cdef RadialSolverDiffeqArgStruct* diffeq_args_ptr = &diffeq_args
    cdef PreEvalFunc diffeq_preeval_ptr               = NULL
    cdef DiffeqFuncType layer_diffeq                  = NULL
    cdef double* y0_ptr                               = NULL

    # The diffeq needs to be able to call the EOS solution which is different for each layer.
    # We will load the EOS solution into this argument structure for each layer. For now, set it to null.

    # Set diffeq inputs that do not change with layer
    diffeq_args_ptr.degree_l         = degree_l_dbl
    diffeq_args_ptr.lp1              = degree_l_dbl + 1.0
    diffeq_args_ptr.lm1              = degree_l_dbl - 1.0
    diffeq_args_ptr.llp1             = degree_l_dbl * (degree_l_dbl + 1.0)
    diffeq_args_ptr.G                = G_to_use
    diffeq_args_ptr.grav_coeff       = 4.0 * d_PI_DBL * G_to_use
    diffeq_args_ptr.frequency        = frequency
    diffeq_args_ptr.layer_index      = 0
    diffeq_args_ptr.eos_solution_ptr = eos_solution_storage_ptr

    # The constant vectors are the same size as the number of solutions in the layer. But since the largest they can
    #  ever be is 3, it is more efficient to just preallocate them on the stack rather than dynamically allocate them
    #  on the heap.
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

    # Variables used to solve the linear equation at the planet's surface.
    # Info = flag set by the solver. Set equal to -999. This will indicate that the solver has not been called yet.
    cdef int bc_solution_info = -999

    # The radial solver can not start at r=0 (singularity) and the higher in the planet the better the stability of
    # the solution. The flip side is that the more of the planet you skip the less accurate your results will be. 
    # Below we determine a reasonable starting radius for the radial solver solution. 
    if starting_radius == 0.0:
        # User did not provide a starting radius.
        # Use a model involving planet radius and degree l to determine a reasonable choice. 
        # This formulism is based on Hilary Martens thesis and LoadDef manual.
        starting_radius = planet_radius * start_radius_tolerance**(1. / degree_l_dbl)

        # Ensure the starting radius is not too close to the surface of the planet.
        starting_radius = fmin(starting_radius, 0.95 * planet_radius)
    
    # Other physical properties at this starting radius (these will be set later)
    cdef double starting_gravity       = d_NAN_DBL
    cdef double starting_density       = d_NAN_DBL
    cdef double complex starting_shear = d_NAN_DBL
    cdef double complex starting_bulk  = d_NAN_DBL

    # Determine which layer this starting radius resides in. We will skip the lower layers
    cdef size_t start_layer_i           = 0
    cdef size_t last_index_before_start = 0
    cdef size_t start_index_in_layer    = 0
    last_radius_check = 0.0
    for current_layer_i in range(num_layers):
        # Check if the radius is in this layer
        layer_upper_radius = eos_solution_storage_ptr.upper_radius_bylayer_vec[current_layer_i]
        if current_layer_i == 0:
            last_layer_upper_radius = 0.0
        else:
            last_layer_upper_radius = eos_solution_storage_ptr.upper_radius_bylayer_vec[current_layer_i - 1]

        if last_layer_upper_radius < starting_radius <= layer_upper_radius:
            # It is!
            start_layer_i = current_layer_i
            first_slice_index = first_slice_index_by_layer_vec[current_layer_i]
            
            # Now find the last radial slice before the starting radius
            start_index_in_layer = 0
            for slice_i in range(first_slice_index, first_slice_index + num_slices_by_layer_vec[current_layer_i]):
                radius_check = radius_array_ptr[slice_i]
                if last_radius_check < starting_radius <= radius_check:
                    if slice_i == 0:
                        last_index_before_start = 0
                    else:
                        # We found that the starting radius is in-between this slice and the last slice.
                        # We want to set the last index to the slice before this one.
                        last_index_before_start = slice_i - 1
                    break
                else:
                    start_index_in_layer += 1
                    last_radius_check = radius_check
            break
        else:
            last_radius_check = last_layer_upper_radius

    # Step through the solution vector and NAN out data below the starting radius
    for slice_i in range(last_index_before_start + 1):
        for ytype_i in range(num_ytypes):
            for y_i in range(MAX_NUM_Y):
                solution_ptr[slice_i * MAX_NUM_Y * num_ytypes + ytype_i * MAX_NUM_Y + y_i] = cmplx_NAN

    for current_layer_i in range(start_layer_i, num_layers):
        # Get layer's index information
        layer_slices = num_slices_by_layer_vec[current_layer_i]
        if current_layer_i == start_layer_i:
            # The first slice index is not going to actually be the bottom-most index
            # Instead it is the one right at or above our starting radius.
            first_slice_index = last_index_before_start + 1

            # When we loop through slices we only want to loop between the starting slice and the top of the layer
            layer_slices -= start_index_in_layer
        else:
            first_slice_index = first_slice_index_by_layer_vec[current_layer_i]

        # Get solution and y information
        num_sols   = num_solutions_by_layer_ptr[current_layer_i]
        num_ys     = 2 * num_sols
        num_ys_dbl = 2 * num_ys

        # Setup pointer array slices starting at the start of this layer (either base or at starting index)
        layer_radius_ptr    = &radius_array_ptr[first_slice_index]
        layer_density_ptr   = &density_array_ptr[first_slice_index]
        layer_gravity_ptr   = &gravity_array_ptr[first_slice_index]
        layer_shear_mod_ptr = &complex_shear_array_ptr[first_slice_index]
        layer_bulk_mod_ptr  = &complex_bulk_array_ptr[first_slice_index]

        # Get physical parameters at the top and bottom of the layer
        if current_layer_i == start_layer_i:
            # In the first layer we can not use the physical properties at the bottom of the arrays
            # because the starting radius may not be at the bottom.
            # Even worse, it may not be at any of the slice indices that are stored.
            # To get the most accurate result we need to perform an interpolation to find various properties at
            # this starting radius. 
            eos_solution_storage_ptr.call(current_layer_i, starting_radius, eos_interp_array_ptr)

            # Save the values, look at the "TidalPy.Material.eos.eos_solution_.hpp" to see how these are saved. 
            # We are storing these in function-global variables because they will be used again during collapse
            starting_gravity = eos_interp_array_ptr[0]
            # TODO: Sometimes at very small r the g can be negative. Probably an issue with the EOS but for now just put a force check in?
            if starting_gravity < d_EPS_DBL:
                starting_gravity = d_EPS_DBL
            starting_density = eos_interp_array_ptr[4]
            starting_shear   = cf_build_dblcmplx(eos_interp_array_ptr[5], eos_interp_array_ptr[6])
            starting_bulk    = cf_build_dblcmplx(eos_interp_array_ptr[7], eos_interp_array_ptr[8])
            
            # Now set local variables used in this loop
            radius_lower  = starting_radius
            gravity_lower = starting_gravity
            density_lower = starting_density
            shear_lower   = starting_shear
            bulk_lower    = starting_bulk
        else:
            # Otherwise we can just use the values at the base of the layer
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

        radial_span_ptr[0] = radius_lower
        radial_span_ptr[1] = radius_upper

        # Determine max step size (if not provided by user)
        if max_step_from_arrays:
            # Maximum step size during integration can not exceed ~1/3 of the layer size.
            max_step_to_use = 0.33 * (radial_span_ptr[1] - radial_span_ptr[0])

        # Get assumptions for layer
        layer_type      = layer_types_ptr[current_layer_i]
        layer_is_static = is_static_by_layer_ptr[current_layer_i]
        layer_is_incomp = is_incompressible_by_layer_ptr[current_layer_i]

        # Determine rtols and atols for this layer.
        # Scale rtols by layer type
        for y_i in range(num_ys):
            # Default is that each layer's rtol and atol equal user input.
            # TODO: Change up the tolerance scaling between real and imaginary?
            layer_rtol_real = integration_rtol
            layer_rtol_imag = integration_rtol
            layer_atol_real = integration_atol
            layer_atol_imag = integration_atol

            if scale_rtols_by_layer_type:
                # Certain layer assumptions can affect solution stability so use additional scales on the relevant rtols
                # TODO test that these scales are the best. See Issue #44
                if layer_type == 0:
                    # Solid layer
                    # Scale y2 and y3 by 0.1
                    if (y_i == 1) or (y_i == 2):
                        # Scale both the real and imaginary portions by the same amount.
                        layer_rtol_real *= 0.1
                        layer_rtol_imag *= 0.1
                else:
                    # Liquid layer
                    if not layer_is_static:
                        # Scale dynamic liquid layer's y2 by additional 0.01.
                        if (y_i == 1):
                            # Scale both the real and imaginary portions by the same amount.
                            layer_rtol_real *= 0.01
                            layer_rtol_imag *= 0.01
            # Populate rtol and atol pointers.
            rtols_ptr[2 * y_i]     = layer_rtol_real
            rtols_ptr[2 * y_i + 1] = layer_rtol_imag
            atols_ptr[2 * y_i]     = layer_atol_real
            atols_ptr[2 * y_i + 1] = layer_atol_imag

        # Set initial y to nan for now.
        # OPT: this is for debugging purposes. could likely be commented out in the future for small performance gain.
        for y_i in range(36):
            if y_i < 18:
                initial_y_ptr[y_i] = cmplx_NAN
            initial_y_only_real_ptr[y_i] = d_NAN_DBL
        
        if current_layer_i == start_layer_i:
            # In the first layer. Use initial condition function to find initial conditions
            cf_find_starting_conditions(
                &solution_storage_ptr.success,
                solution_storage_ptr.message_ptr,
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
                MAX_NUM_Y, 
                initial_y_ptr,  # Modified Variable
                starting_y_check
                )

            if not solution_storage_ptr.success:
                solution_storage_ptr.error_code = -10
                break
                
        else:
            # Not in the first layer. Use interface function to find initial conditions
            layer_below_type      = layer_types_ptr[current_layer_i - 1]
            layer_below_is_static = is_static_by_layer_ptr[current_layer_i - 1]
            layer_below_is_incomp = is_incompressible_by_layer_ptr[current_layer_i - 1]

            # Find gravity at the base interface using bottom of this layer and top of previous.
            interface_gravity = 0.5 * (gravity_lower + last_layer_upper_gravity)

            # Find the density needed for some initial conditions.
            if (layer_type == 0) and (layer_below_type == 0):
                # Both layers are solid. A liquid interface density is not needed.
                static_liquid_density = d_NAN_DBL
            elif not (layer_type == 0) and (layer_below_type == 0):
                # Layer below is solid, this layer is liquid. Use its density.
                static_liquid_density = density_lower
            elif (layer_type == 0) and not (layer_below_type == 0):
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
                    static_liquid_density = d_NAN_DBL

            # Find the starting values for this layer using the results a the top of the previous layer + an interface
            #  function.
            cf_solve_upper_y_at_interface(
                uppermost_y_per_solution_ptr,
                initial_y_ptr,
                layer_below_num_sols,
                num_sols,
                MAX_NUM_Y,
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

        # Reset the uppermost y value array
        for i in range(18):
            uppermost_y_per_solution_ptr[i] = cmplx_NAN

        # Change initial conditions into 2x real values instead of complex for integration
        for solution_i in range(num_sols):
            for y_i in range(num_ys):
                dcomplex_tmp = initial_y_ptr[solution_i * MAX_NUM_Y + y_i]
                initial_y_only_real_ptr[solution_i * MAX_NUM_Y_REAL + 2 * y_i]     = dcomplex_tmp.real
                initial_y_only_real_ptr[solution_i * MAX_NUM_Y_REAL + 2 * y_i + 1] = dcomplex_tmp.imag

        # Find correct diffeq
        layer_diffeq = cf_find_layer_diffeq(layer_type, layer_is_static, layer_is_incomp)

        # Set the layer index to current layer in the additional diffeq arg struct
        diffeq_args_ptr.layer_index = current_layer_i

        # Get storage pointer for this layer
        storage_by_solution_ptr = main_storage_ptr[current_layer_i]

        # Solve for each solution
        for solution_i in range(num_sols):
            y0_ptr = &initial_y_only_real_ptr[solution_i * MAX_NUM_Y_REAL]

            ###### Integrate! #######
            integration_solution = cysolve_ivp(
                layer_diffeq,            # Differential equation [DiffeqFuncType]
                radial_span_ptr,         # Radial span [const double*]
                y0_ptr,                  # y0 array [const double*]
                num_ys_dbl,              # Number of ys [size_t]
                integration_method,      # Integration method [int]
                d_NAN_DBL,                     # Relative Tolerance (as scalar) [double]
                d_NAN_DBL,                     # Absolute Tolerance (as scalar) [double]
                <char*>diffeq_args_ptr,  # Extra input args to diffeq [char*]
                sizeof_args,             # Size of additional argument structure [size_t]
                0,                       # Number of extra outputs tracked [size_t]
                max_num_steps,           # Max number of steps (0 = find good value) [size_t]
                max_ram_MB,              # Max amount of RAM allowed [size_t]
                False,                   # Use dense output [bint]
                layer_radius_ptr,        # Interpolate at radius array [double*]
                layer_slices,            # Size of interpolation array [size_t]
                diffeq_preeval_ptr,      # Pre-eval function used in diffeq [PreEvalFunc]
                rtols_ptr,               # Relative Tolerance (as array) [double*]
                atols_ptr,               # Absolute Tolerance (as array) [double*]
                max_step_to_use,         # Maximum step size [double]
                0.0,                     # Initial step size (0 = find good value) [doub;e]
                expected_size            # Expected final integration size (0 = find good value) [size_t]
                )
            integration_solution_ptr = integration_solution.get()
            #########################

            # Store diagnostic data
            solution_storage_ptr.shooting_method_steps_taken_vec[(3 * current_layer_i) + solution_i] = \
                integration_solution_ptr.steps_taken

            # Check for problems
            if not integration_solution_ptr.success:
                # Problem with integration.
                solution_storage_ptr.error_code = -11
                solution_storage_ptr.success = False
                sprintf(message_ptr, 'RadialSolver.ShootingMethod:: Integration problem at layer %d; solution %d:\n\t%s.\n', current_layer_i, solution_i, integration_solution_ptr.message_ptr)
                solution_storage_ptr.set_message(message_ptr)
                if verbose :
                    printf(message_ptr)
                return solution_storage_ptr.error_code

            # If no problems, store results.
            integrator_data_ptr = &integration_solution_ptr.solution[0]
            # Need to make a copy because the solver pointers will be reallocated during the next solution.
            # Get storage pointer for this solution
            storage_by_y_ptr = storage_by_solution_ptr[solution_i]

            for y_i in range(num_ys):
                for slice_i in range(layer_slices):
                    # Convert 2x real ys to 1x complex ys
                    storage_by_y_ptr[num_ys * slice_i + y_i] = cf_build_dblcmplx(
                        integrator_data_ptr[num_ys_dbl * slice_i + (2 * y_i)],
                        integrator_data_ptr[num_ys_dbl * slice_i + (2 * y_i) + 1]
                        )

                # Store top most result for initial condition for the next layer
                # slice_i should already be set to the top of this layer after the end of the previous loop.
                uppermost_y_per_solution_ptr[solution_i * MAX_NUM_Y + y_i] = storage_by_y_ptr[num_ys * slice_i + y_i]
            
            # Decrement the shared pointer for the solution
            integration_solution.reset()

        if solution_storage_ptr.error_code != 0:
            # Error was encountered during integration
            solution_storage_ptr.success = False
            return solution_storage_ptr.error_code

        # Prepare for next layer
        layer_below_num_sols     = num_sols
        layer_below_num_ys       = num_ys
        last_layer_upper_gravity = gravity_upper
        last_layer_upper_density = density_upper

    if solution_storage_ptr.error_code < 0 or not solution_storage_ptr.success:
        solution_storage_ptr.success = False
        if integration_solution_ptr:
            sprintf(message_ptr, 'RadialSolver.ShootingMethod:: Integration failed:\n\t%s.\n', integration_solution_ptr.message_ptr)
            solution_storage_ptr.set_message(message_ptr)

        if verbose:
            printf(message_ptr)
        return solution_storage_ptr.error_code

    else:
        # No errors. Proceed with collapsing all sub-solutions into final full solution.
        strcpy(message_ptr, 'Integration completed for all layers. Beginning solution collapse.\n')

        for ytype_i in range(num_ytypes):
            sprintf(message_ptr, 'Collapsing radial solutions for "%d" solver.\n', ytype_i)

            # Reset variables for this solver
            bc_solution_info = -999
            constant_vector_ptr[0] = cmplx_NAN
            constant_vector_ptr[1] = cmplx_NAN
            constant_vector_ptr[2] = cmplx_NAN
            layer_above_lower_gravity   = d_NAN_DBL
            layer_above_lower_density   = d_NAN_DBL
            liquid_density_at_interface = d_NAN_DBL
            layer_above_type            = 9
            layer_above_is_static       = False
            layer_above_is_incomp       = False

            if verbose:
                printf(message_ptr)
            # Collapse the multiple solutions for each layer into one final combined solution.

            # Work from the surface down to the starting layer.
            for current_layer_i in range(num_layers - start_layer_i):
                layer_i_reversed = num_layers - (current_layer_i + 1)

                # Pull out layer information.
                layer_slices       = num_slices_by_layer_vec[layer_i_reversed]
                if layer_i_reversed == start_layer_i:
                    # The first slice index is not going to actually be the bottom-most index
                    # Instead it is the one right at or above our starting radius.
                    first_slice_index = last_index_before_start + 1

                    # When we loop through slices we only want to loop between the starting slice and the top of the layer
                    layer_slices -= (last_index_before_start + 1)
                else:
                    first_slice_index = first_slice_index_by_layer_vec[layer_i_reversed]

                # Get solution and y information
                num_sols     = num_solutions_by_layer_ptr[layer_i_reversed]
                num_ys       = 2 * num_sols

                # Setup pointer array slices starting at this layer's beginning
                layer_radius_ptr    = &radius_array_ptr[first_slice_index]
                layer_density_ptr   = &density_array_ptr[first_slice_index]
                layer_gravity_ptr   = &gravity_array_ptr[first_slice_index]
                layer_bulk_mod_ptr  = &complex_bulk_array_ptr[first_slice_index]
                layer_shear_mod_ptr = &complex_shear_array_ptr[first_slice_index]

                # Get physical parameters at the top and bottom of the layer
                if layer_i_reversed == start_layer_i:
                    # In the starting layer.
                    # Found earlier via interpolation at the starting radius
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

                # Get assumptions for layer
                layer_type      = layer_types_ptr[layer_i_reversed]
                layer_is_static = is_static_by_layer_ptr[layer_i_reversed]
                layer_is_incomp = is_incompressible_by_layer_ptr[layer_i_reversed]

                # Get full solutions for this layer
                storage_by_solution_ptr = main_storage_ptr[layer_i_reversed]

                # Get radial solution values at the top of the layer
                for solution_i in range(num_sols):
                    for y_i in range(num_ys):
                        uppermost_y_per_solution_ptr[solution_i * MAX_NUM_Y + y_i] = \
                            storage_by_solution_ptr[solution_i][(layer_slices - 1) * num_ys + y_i]

                if current_layer_i == 0:
                    # Working on surface (uppermost) layer -- Apply surface boundary conditions.
                    cf_apply_surface_bc(
                        constant_vector_ptr,  # Modified Variable
                        &bc_solution_info,  # Modified Variable
                        bc_pointer,
                        uppermost_y_per_solution_ptr,
                        surface_gravity,
                        G_to_use,
                        num_sols,
                        MAX_NUM_Y,
                        ytype_i,
                        layer_type,
                        layer_is_static,
                        layer_is_incomp
                        )

                    # Check that the boundary condition was successfully applied.
                    if bc_solution_info != 0:
                        solution_storage_ptr.error_code = -12
                        sprintf(message_ptr, 'RadialSolver.ShootingMethod:: Error encountered while applying surface boundary condition. ZGESV code: %d.\nThe solutions may not be valid at the surface.\n', bc_solution_info)
                        solution_storage_ptr.set_message(message_ptr)
                        if verbose:
                            printf(message_ptr)
                        return solution_storage_ptr.error_code
                else:
                    # Working on interior layers. Will need to find the constants of integration based on the layer above.
                    cf_top_to_bottom_interface_bc(
                        constant_vector_ptr,  # Modified Variable
                        layer_above_constant_vector_ptr,
                        uppermost_y_per_solution_ptr,
                        gravity_upper, layer_above_lower_gravity,
                        density_upper, layer_above_lower_density,
                        layer_type, layer_above_type,
                        layer_is_static, layer_above_is_static,
                        layer_is_incomp, layer_above_is_incomp,
                        num_sols, MAX_NUM_Y
                        )

                # Use constant vectors to find the full y from all of the solutions in this layer
                cf_collapse_layer_solution(
                    solution_ptr,  # Modified Variable
                    constant_vector_ptr,
                    storage_by_solution_ptr,
                    layer_radius_ptr,
                    layer_density_ptr,
                    layer_gravity_ptr,
                    frequency,
                    first_slice_index,
                    layer_slices,
                    num_sols,
                    MAX_NUM_Y,
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

            # Ready for next y-type

    # Free memory
    # Deconstruct main solution pointer
    # Main storage pointers are structured like [current_layer_i][solution_i][y_i + slice_i]
    # Then main storage
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
        solution_storage_ptr.set_message('RadialSolver.ShootingMethod: Completed without any noted issues.')

    # Done!
    return solution_storage_ptr.error_code
