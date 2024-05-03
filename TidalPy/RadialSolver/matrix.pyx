# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport NAN, isnan
from libc.stdio cimport printf
from libc.stdlib cimport exit

from scipy.linalg.cython_lapack cimport zgesv
from CyRK.utils.utils cimport allocate_mem, free_mem

from TidalPy.utilities.constants_x cimport G
from TidalPy.utilities.math.complex cimport cf_build_dblcmplx, cmplx_NAN
from TidalPy.utilities.dimensions.nondimensional cimport cf_non_dimensionalize_physicals, cf_redimensionalize_physicals
from TidalPy.RadialSolver.boundaries.surface_bc cimport cf_get_surface_bc
from TidalPy.RadialSolver.matrix_types.solid_matrix cimport cf_fundamental_matrix
from TidalPy.RadialSolver.solver cimport RadialSolverSolution, MAX_NUM_Y



cdef RadialSolverSolution cf_matrix_propagate(
            size_t total_slices,
            double* radius_array_ptr,
            double* density_array_ptr,
            double* gravity_array_ptr,
            double* bulk_modulus_ptr,
            double complex* complex_shear_modulus_array_ptr,
            double frequency,
            double planet_bulk_density,
            # TODO: In the future the propagation matrix should take in layer types and multiple layers
            # size_t num_layers,
            # int* layer_types_ptr,
            # int* is_static_by_layer_ptr,
            # int* is_incompressible_by_layer_ptr,
            # double* upper_radius_by_layer_ptr,
            size_t num_bc_models,
            size_t* bc_models_ptr,
            unsigned int degree_l = 2,
            unsigned char core_condition = 0,
            bint nondimensionalize = True,
            bint verbose = False,
            bint raise_on_fail = False
            ) noexcept nogil:

        # Setup
        cdef size_t r_i, j, k, jj, ytype_i
        cdef size_t last_index_shift_36, index_shift_36, last_index_shift_18, index_shift_18, index_shift_max_y
        cdef bint error = False

        # Nondimensional variables
        cdef double mean_radius = radius_array_ptr[total_slices - 1]
        cdef double planet_radius_to_use = NAN
        cdef double bulk_density_to_use = NAN
        cdef double frequency_to_use = NAN
        cdef double G_to_use = NAN

        if nondimensionalize:
            cf_non_dimensionalize_physicals(
                total_slices, frequency, mean_radius, planet_bulk_density, radius_array_ptr, density_array_ptr,
                gravity_array_ptr, bulk_modulus_ptr, complex_shear_modulus_array_ptr,
                &planet_radius_to_use, &bulk_density_to_use, &frequency_to_use, &G_to_use
                )
            
            # Ensure that no errors occured during the non-dim process
            if isnan(planet_radius_to_use) or isnan(bulk_density_to_use) or isnan(frequency_to_use) or isnan(G_to_use):
                raise ValueError('NaNs encountered after non-dimensionalize call.')
        else:
            planet_radius_to_use = mean_radius
            bulk_density_to_use = planet_bulk_density
            frequency_to_use = frequency
            G_to_use = G

        # Find boundary condition at the top of the planet -- this is dependent on the forcing type.
        #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
        # The [5] represents the maximum number of solvers that can be invoked with a single call to radial_solver
        cdef double degree_l_dbl = <double>degree_l
        cdef size_t max_num_solutions = 5
        cdef size_t num_ytypes = num_bc_models
        cdef str solver_name

        # 15 = 5 (max_num_solutions) * 3 (number of surface conditions)
        cdef double[15] boundary_conditions
        cdef double* bc_pointer = &boundary_conditions[0]
        cf_get_surface_bc(
            bc_pointer,  # Changed parameter
            bc_models_ptr,
            num_ytypes,
            planet_radius_to_use,
            bulk_density_to_use,
            degree_l_dbl
            )

        # Define memory for our fundamental matricies. Each matrix is 6 x 6 x radial slices
        # These have to be heap allocated because we do not know the number of radial slices at compile time (and it could be large)
        cdef size_t matrix_size = 6 * 6 * total_slices
        cdef double complex* fundamental_mtx_ptr
        cdef double complex* inverse_fundamental_mtx_ptr
        cdef double complex* derivative_mtx_ptr
        cdef double complex* propagation_mtx_ptr

        fundamental_mtx_ptr = <double complex *>allocate_mem(
            sizeof(double complex) * matrix_size,
            "`fundamental_mtx_ptr` (cf_matrix_propagate)"
            )
        inverse_fundamental_mtx_ptr = <double complex *>allocate_mem(
            sizeof(double complex) * matrix_size,
            "`inverse_fundamental_mtx_ptr` (cf_matrix_propagate)"
            )
        derivative_mtx_ptr = <double complex *>allocate_mem(
            sizeof(double complex) * matrix_size,
            "`derivative_mtx_ptr` (cf_matrix_propagate)"
            )
        
        # Propagation matrix has 6 rows but only 3 columns.
        cdef size_t prop_mat_size = 6 * 3 * total_slices
        propagation_mtx_ptr = <double complex *>allocate_mem(
            sizeof(double complex) * prop_mat_size,
            "`propagation_mtx_ptr` (cf_matrix_propagate)"
            )
        
        # Popualte matricies with the correct layer type. 
        # TODO: Currently only solid, static, incompressible layers are supported for matrix propagation.
        cf_fundamental_matrix(
            total_slices,
            radius_array_ptr,
            density_array_ptr,
            gravity_array_ptr,
            complex_shear_modulus_array_ptr,
            fundamental_mtx_ptr,  # Changed variable
            inverse_fundamental_mtx_ptr,  # Changed variable
            derivative_mtx_ptr,  # Changed variable
            degree_l,
            G_to_use
            ) 

        # Initialize the base of the propagation matrix to the initial conditions
        ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."

        # TODO Add more of these.
        if core_condition == 0:
            # Henning & Hurford (2014): "At the core, a special seed matrix Bcore is created with only three columns,
            # equal to the first, second, and third columns of Y for the properties at the base layer."
            for j in range(6):
                for k in range(3):
                    propagation_mtx_ptr[j * 6 + k] = fundamental_mtx_ptr[j * 6 + k]
        elif core_condition == 1:
            # Roberts & Nimmo (2008): liquid innermost zone.
            for j in range(6):
                for k in range(3):
                    if (j == 2) and (k == 0):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(1., 0.)
                    elif (j == 3) and (k == 1):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(1., 0.)
                    elif (j == 5) and (k == 2):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(1., 0.)
                    else:
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(0., 0.)
        elif core_condition == 2:
            # Solid innermost zone
            for j in range(6):
                for k in range(3):
                    if (j == 0) and (k == 0):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(1., 0.)
                    elif (j == 1) and (k == 1):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(1., 0.)
                    elif (j == 2) and (k == 2):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(1., 0.)
                    else:
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(0., 0.)
        elif core_condition == 3:
            # Solid innermost zone
            for j in range(6):
                for k in range(3):
                    if (j == 0) and (k == 0):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(0.05, 0.)
                    elif (j == 1) and (k == 1):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(0.01, 0.)
                    elif (j == 5) and (k == 2):
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(1., 0.)
                    else:
                        propagation_mtx_ptr[j * 6 + k] = cf_build_dblcmplx(0., 0.)
        else:
            printf("Unknown starting core conditions encountered in `cf_matrix_propagate`: %d (acceptible values: 0, 1, 2, 3)", core_condition)
            exit(-1)

        # Step through the planet's shells and build the propagation matrix
        cdef double complex temp_cmplx
        cdef double complex[18] temp_matrix
        cdef double complex* temp_matrix_ptr = &temp_matrix[0]
        cdef double complex[9] surface_matrix
        cdef double complex* surface_matrix_ptr = &surface_matrix[0]

        # Need to find the inverse of the surface matrix (3 x 3)
        # We will use a linear equation solver to solve: surface_solution = surface_matrix^-1 @ surface_bc
        # Re-arranged: surface_matrix @ surface_solution = surface_bc
        # ZGESV computes the solution to system of linear equations A * X = B for GE matrices
        # See https://www.netlib.org/lapack/explore-html/d8/da6/group__gesv_ga0850dc117a6c7ec3cb64905d5de1cd23.html#ga0850dc117a6c7ec3cb64905d5de1cd23
        
        # Create the ZGESV "A" variable
        cdef double complex[9] surface_matrix
        cdef double complex* surface_matrix_ptr = &surface_matrix[0]
        cdef double complex[9] surface_matrix_copy
        cdef double complex* surface_matrix_copy_ptr = &surface_matrix_copy[0]

        # Create the ZGESV "X" variable
        cdef double complex[3] surface_solution
        cdef double complex* surface_solution_ptr = &surface_solution[0]

        # Create a copy of ZGESV "B" surface boundary condition because it is overwritten with each call of zgesv
        cdef double complex[3] bc_copy
        cdef double complex* bc_copy_ptr = &bc_copy[0]

        # Initialize surface matrix and temp_matrix to nan to help with debugging
        for j in range(18):
            if j < 3:
                surface_solution_ptr[j] = cmplx_NAN
                bc_copy_ptr[j] = cmplx_NAN
            if j < 9:
                surface_matrix_ptr[j] = cmplx_NAN
                surface_matrix_copy_ptr[j] = cmplx_NAN
            temp_matrix_ptr[j] = cmplx_NAN

        for r_i in range(1, total_slices):
            
            # Need to start the index for this radial slice. The shift is based on the size of the respective matrix.
            # Fundamental matrix is 6x6
            index_shift_36 = r_i * 36
            last_index_shift_36 = (r_i - 1) * 36
            # Propagation matrix is 3x3
            index_shift_18 = r_i * 18
            last_index_shift_18 = (r_i - 1) * 18

            # The function we are performing here is:
            # P_{i} = Y_{i} @ ( Y_{i-1}^{-1} @ P_{i-1} )

            # Perform the first matrix multiplication
            # A = Y_{i-1}^{-1} @ P_{i-1}
            for j in range(6):
                for k in range(3):
                    temp_cmplx = cf_build_dblcmplx(0., 0.)
                    for jj in range(6):
                        temp_cmplx += (
                            inverse_fundamental_mtx_ptr[last_index_shift_36 + j * 6 + jj] * 
                            propagation_mtx_ptr[last_index_shift_18 + jj * 6 + k]
                            )
                    temp_matrix_ptr[j * 6 + k] = temp_cmplx
            
            # Now perform the outer matrix multiplication
            # P_{i} = Y_{i} @ A
            for j in range(6):
                for k in range(3):
                    temp_cmplx = cf_build_dblcmplx(0., 0.)
                    for jj in range(6):
                        temp_cmplx += (
                            fundamental_mtx_ptr[index_shift_36 + j * 6 + jj] * 
                            temp_matrix_ptr[jj * 6 + k]
                            )
                    propagation_mtx_ptr[index_shift_18 + j * 6 + k] = temp_cmplx
        
            # We need to define a matrix that equals the propagation matrix at the surface value.
            # Surface condition matrix is a 3x3 matrix of the top-most shell of the propagation matrix's rows [3, 4, 6]
            if r_i == (total_slices - 1):
                index_shift_18 = r_i * 18
                for i in range(3):
                    surface_matrix_ptr[0 + i] = propagation_mtx_ptr[index_shift_18 + (2 * 6) + i]
                    surface_matrix_ptr[3 + i] = propagation_mtx_ptr[index_shift_18 + (3 * 6) + i]
                    surface_matrix_ptr[6 + i] = propagation_mtx_ptr[index_shift_18 + (5 * 6) + i]
        
        # Other information required by zgesv
        # Size of matrix; 3x3
        cdef int mat_size = 3
        cdef int* mat_size_ptr = &mat_size
        # Size of bc 3x1
        cdef int bc_columns = 1
        cdef int* bc_columns_ptr = &bc_columns

        # IPIV = Integer pivot array that is an additional output provided by ZGESV. It is not used but must be provided.
        #  It must be at least as large as the largest dimension of the input matrix, for this work that is 3.
        cdef int[10] lapack_ipiv
        cdef int* lapack_ipiv_ptr = &lapack_ipiv[0]

        # Info = flag set by the solver. 
        cdef int bc_solution_info = -999
        cdef int* bc_solution_info_ptr = &bc_solution_info

        # Build radial solution storage
        cdef double complex* y_sv_ptr
        y_sv_ptr = <double complex *>allocate_mem(
            sizeof(double complex) * MAX_NUM_Y * total_slices,
            "`y_sv_ptr` (cf_matrix_propagate)"
            )
        
        # Build solution
        cdef RadialSolverSolution solution
        solution = RadialSolverSolution(total_slices, bc_models_ptr, num_bc_models)
        cdef double complex* solution_ptr = solution.full_solution_ptr

# slice_i_shifted * num_output_ys + (lhs_y_index + 4)

        for ytype_i in range(num_bc_models):
            
            # New y-type being solved (tidal, loading, free)
            # Set linear solver flag to -999. This will indicate that the solver has not been called yet.
            bc_solution_info = -999

            # Set / reset the values of the RHS and LHS of the equation
            for j in range(9):
                if j < 3:
                    bc_copy_ptr[j] = bc_pointer[ytype_i * 3 + j]
                surface_matrix_copy_ptr[j] = surface_matrix_ptr[j]

            zgesv(
                mat_size_ptr,
                bc_columns_ptr,
                surface_matrix_ptr,  # Modified Variable
                mat_size_ptr,
                lapack_ipiv_ptr,  # Modified Variable
                bc_copy_ptr,  # Modified Variable - Both the RHS of the equation and where the result will be stored.
                mat_size_ptr,
                bc_solution_info_ptr  # Modified Variable
                )
            # bc_copy_ptr
        
            # Check for errors
            # TODO: Convert to char pointer?
            if bc_solution_info != 0:
                feedback_str = \
                    (f'Error encountered while applying surface boundary condition. ZGESV code: {bc_solution_info}'
                    f'\nThe solutions may not be valid at the surface.')
                if verbose:
                    print(feedback_str)
                if raise_on_fail:
                    raise RuntimeError(feedback_str)
                error = True
            
            for i in range(total_slices):
                
                index_shift_18 = r_i * 18
                index_shift_max_y = r_i * MAX_NUM_Y

                # Perform matrix multiplication: prop_matrix @ surface_solution 
                for j in range(6):
                    temp_cmplx = cf_build_dblcmplx(0., 0.)
                    for jj in range(3):
                        temp_cmplx += (
                            propagation_mtx_ptr[index_shift_18 + j * 6 + jj] * 
                            bc_copy_ptr[jj]
                            )
                    y_sv_ptr[index_shift_max_y + j] = temp_cmplx

                    slice_i_shifted * num_output_ys + (lhs_y_index + 4)
                    solution_ptr
                

        # Invert the surface matrix and solve using the surface boundary condition
        surface_matrix_inv = np.linalg.inv(surface_matrix)
        surface_solution = surface_matrix_inv @ surface_bc
            
        S = S

        # Redimensionalize
        if nondimensionalize:
            cf_redimensionalize_physicals(
                total_slices, frequency, mean_radius, planet_bulk_density, radius_array_ptr, density_array_ptr,
                gravity_array_ptr, bulk_modulus_ptr, complex_shear_modulus_array_ptr,
                &planet_radius_to_use, &bulk_density_to_use, &frequency_to_use, &G_to_use
                )
            


        # Free memory
        free_mem(fundamental_mtx_ptr)
        free_mem(inverse_fundamental_mtx_ptr)
        free_mem(derivative_mtx_ptr)
        free_mem(propagation_mtx_ptr)
        free_mem(y_sv_ptr)



    fundamental_matrix: np.ndarray, fundamental_matrix_inverse: np.ndarray, derivative_matrix: np.ndarray,
    inner_boundary_condition: np.ndarray, world_radius: float,
    order_l: int = 2
    ) -> np.ndarray:
    """ This function will propagate the incompressible tidal equations, via the fundamental matrix, through a world or
    layers sub-shells.

    Assumptions:
    - Each shell is homogeneous (However, it can have different properties from the shells above/below it).

    References
    ----------
    See appendix of HH14 and lines ~2355 of thermal.h in ID

    Parameters
    ----------
    fundamental_matrix : np.ndarray
        Fundamental (incompressible) Matrix (6 x 6 x N); See fundamental.py
    fundamental_matrix_inverse : np.ndarray
        Inverse of the fundamental (incompressible) Matrix (6 x 6 x N); See fundamental.py
    derivative_matrix : np.ndarray
        Derivative matrix, A, defined by the function dy/dr = A dot y
    inner_boundary_condition : np.ndarray
        Boundary condition Matrix (6 x 3) of the tidal problem at the inner surface.
    world_radius : float
        Radius of the world [m]
    order_l : int = 2
        Tidal harmonic order

    Returns
    -------
    tidal_y : np.ndarray
        Matrix [6 x N] of tidal solutions. See decompression.py on how useful information is extracted.

    """
    num_shells = fundamental_matrix.shape[2]

    # Find propagation matrix
    #    Prop matrix has 6 rows but only 3 columns.
    propagation_mtx = np.zeros((6, 3, num_shells), dtype=np.complex128)
    #    It is initialized with the inner boundary condition matrix.
    propagation_mtx[:, :, 0] = inner_boundary_condition

    # Step through the planet's shells
    for i in range(num_shells):
        if i == 0:
            # Skip central shell - boundary conditions provided.
            continue

        # Dot product of between the fundamental matrix and (the dot product between the inverse fundamental matrix one
        #    layer below and the propagation matrix one layer below).
        propagation_mtx[:, :, i] = \
            fundamental_matrix[:, :, i] @ (fundamental_matrix_inverse[:, :, i - 1] @ propagation_mtx[:, :, i - 1])

    # Surface condition matrix is a 3x3 matrix of the top-most shell of the aggregate matrix's rows [3, 4, 6]
    surface_matrix = np.vstack((propagation_mtx[2, :, -1], propagation_mtx[3, :, -1], propagation_mtx[5, :, -1]))

    # The surface boundary conditions are (always?) static and only based on the radius and order-l
    surface_bc = np.zeros((3,), dtype=np.complex128)
    surface_bc[2] = -1. * (2. * order_l + 1.) / world_radius

    # Invert the surface matrix and solve using the surface boundary condition
    surface_matrix_inv = np.linalg.inv(surface_matrix)
    surface_solution = surface_matrix_inv @ surface_bc

    # Using the aggregate matrix, solve for the tidal "y"s
    y_sv = np.empty((6, num_shells), dtype=np.complex128)
    y_derivative_sv = np.empty((6, num_shells), dtype=np.complex128)
    for i in range(num_shells):
        y_sv[:, i] = propagation_mtx[:, :, i] @ surface_solution

        # Calculate the derivatives of the tidal solution with radius
        y_derivative_sv[:, i] = derivative_matrix[:, :, i] @ y_sv[:, i]

    # As discussed in B13 (discussed near their equation 7), SVC16 (and the earlier 2004 book) use a different
    #    convention for tidal_y than is used by Takeuchi and Saito (1972). Since a good chunk of the field follows the
    #    latter, we will do the same. Below are the conversions from SVC16 to TS72
    y = np.empty_like(y_sv)
    y[0, :] = y_sv[0, :]
    y[1, :] = y_sv[2, :]  # Flip y3 and y2
    y[2, :] = y_sv[1, :]  # Flip y3 and y2
    y[3, :] = y_sv[3, :]
    y[4, :] = y_sv[4, :] * -1.
    y[5, :] = y_sv[5, :] * -1.

    return y