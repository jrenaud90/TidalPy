# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

""" Propagation of tidal solution using the fundamental matrix

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
B13   : Beuthe (2013, DOI: 10.1016/j.icarus.2012.11.020)
"""

from libc.math cimport NAN, isnan
from libc.stdio cimport printf, sprintf
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.string cimport strcpy

from scipy.linalg.cython_lapack cimport zgesv
from CyRK.utils.utils cimport allocate_mem, free_mem

from TidalPy.constants cimport d_G
from TidalPy.utilities.math.complex cimport cmplx_zero, cmplx_NAN, cf_build_dblcmplx
from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC
from TidalPy.RadialSolver.boundaries.surface_bc cimport cf_get_surface_bc
from TidalPy.RadialSolver.matrix_types.solid_matrix cimport cf_fundamental_matrix
from TidalPy.RadialSolver.constants cimport MAX_NUM_Y


cdef void cf_matrix_propagate(
        shared_ptr[RadialSolutionStorageCC] solution_storage_sptr,
        double frequency,
        double planet_bulk_density,
        # TODO: In the future the propagation matrix should take in layer types and multiple layers
        # int* layer_types_ptr,
        # int* is_static_by_layer_ptr,
        # int* is_incompressible_by_layer_ptr,
        size_t num_bc_models,
        int* bc_models_ptr,
        double G_to_use = d_G,
        unsigned int degree_l = 2,
        unsigned char core_condition = 0,
        cpp_bool verbose = False,
        cpp_bool raise_on_fail = False
        ) noexcept nogil:

    # Setup
    cdef size_t r_i, i, j, k, jj, ytype_i, slice_i
    cdef size_t last_index_shift_36, index_shift_36, last_index_shift_18, index_shift_18, index_shift_max_y, full_shift

    # Get raw pointer of radial solver storage and eos storage
    cdef RadialSolutionStorageCC* solution_storage_ptr = solution_storage_sptr.get()
    cdef EOSSolutionCC* eos_solution_storage_ptr       = solution_storage_sptr.get().eos_solution_sptr.get()

    strcpy(solution_storage_ptr.message_ptr, "RadialSolver.PropMatrixMethod:: Propagator Matrix Method Called.\n")

    # Pull out key information
    cdef size_t num_layers     = eos_solution_storage_ptr.num_layers
    cdef size_t total_slices   = eos_solution_storage_ptr.radius_array_size
    cdef size_t top_slice_i    = total_slices - 1
    cdef size_t num_interfaces = num_layers - 1

    # Alias pointers to EOS properties
    cdef double* radius_array_ptr  = &eos_solution_storage_ptr.radius_array_vec[0]
    cdef double* gravity_array_ptr = &eos_solution_storage_ptr.gravity_array_vec[0]
    # cdef double* pressure_array_ptr = solution_storage_ptr.pressure_ptr   # Unused
    cdef double* density_array_ptr = &eos_solution_storage_ptr.density_array_vec[0]
    
    # Need to recast the storage's shear/bulk double arrays to double complex for local use
    cdef double complex* complex_shear_array_ptr = <double complex*>&eos_solution_storage_ptr.complex_shear_array_vec[0]
    cdef double complex* complex_bulk_array_ptr  = <double complex*>&eos_solution_storage_ptr.complex_bulk_array_vec[0]

    # Pull out constants
    cdef double planet_radius = radius_array_ptr[top_slice_i]

    # Find boundary condition at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    # The [5] represents the maximum number of solvers that can be invoked with a single call to radial_solver
    cdef double degree_l_dbl      = <double>degree_l
    cdef size_t max_num_solutions = 5
    cdef size_t num_ytypes        = num_bc_models

    # Boundary condition size: 15 = 5 (max_num_solutions) * 3 (number of surface conditions)
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

    # Define memory for our fundamental matricies.
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

    # Populate matricies with the correct layer type. 
    # TODO: Currently only solid, static, incompressible layers are supported for matrix propagation.
    cf_fundamental_matrix(
        total_slices,
        radius_array_ptr,
        density_array_ptr,
        gravity_array_ptr,
        complex_shear_array_ptr,
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
                    propagation_mtx_ptr[j * 6 + k] = cmplx_zero
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
                    propagation_mtx_ptr[j * 6 + k] = cmplx_zero
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
                    propagation_mtx_ptr[j * 6 + k] = cmplx_zero
    else:
        sprintf(solution_storage_ptr.message_ptr, "RadialSolver.PropMatrixMethod:: Unknown starting core conditions encountered in `cf_matrix_propagate`: %d (acceptible values: 0, 1, 2, 3)\n", core_condition)
        solution_storage_ptr.error_code = -20
        if verbose or raise_on_fail:
            printf(solution_storage_ptr.message_ptr)
        if raise_on_fail:
            exit(EXIT_FAILURE)

    # Step through the planet's shells and build the propagation matrix
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

    cdef double complex temp_cmplx
    cdef double complex[18] temp_matrix
    cdef double complex* temp_matrix_ptr = &temp_matrix[0]

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
        # Propagation matrix is 6x3
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
            for i in range(3):
                surface_matrix_ptr[0 + i] = propagation_mtx_ptr[index_shift_18 + (2 * 6) + i]
                surface_matrix_ptr[3 + i] = propagation_mtx_ptr[index_shift_18 + (3 * 6) + i]
                surface_matrix_ptr[6 + i] = propagation_mtx_ptr[index_shift_18 + (5 * 6) + i]
    
    # Next we need to solve the linear equation U = S^-1 @ B
    # Where U is the solution at the surface, S is the surface matrix constructed in the previous step, and B
    # is the surface boundary condition. 
    # To do this we will use the linear equation solver for complex numbers, ZGESV. In order to use this we need to
    # change the equation to A * X = B; or for our example: S * U = B
    # More info on ZGESV: https://www.math.utah.edu/software/lapack/lapack-z/zgesv.html

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

    # Info = flag set by the solver. -999 indicates ZGESV has not been called yet.
    cdef int bc_solution_info = -999
    cdef int* bc_solution_info_ptr = &bc_solution_info
    
    # Used to convert from SVC radial solutions to T&S format
    cdef double complex[6] ts_conversion
    cdef double complex* ts_conversion_ptr = &ts_conversion[0]
    for i in range(6):
        ts_conversion_ptr[i] = cmplx_NAN

    # Various index shifts
    cdef size_t ytype_shift
    cdef size_t solution_slice_ishift

    # Build solution
    cdef double* solution_dbl_ptr = solution_storage_ptr.full_solution_ptr
    # Cast the solution pointer from double to double complex
    cdef double complex* solution_ptr = <double complex*>solution_dbl_ptr

    ytype_i = 0
    while solution_storage_ptr.error_code == 0:

        # New y-type being solved (tidal, loading, free)
        if ytype_i == num_bc_models:
            break

        # Set linear solver flag to -999. This will indicate that the solver has not been called yet.
        bc_solution_info_ptr[0] = -999

        # Set / reset the values of the RHS and LHS of the equation
        # We use copies of these pointers since ZGESV overwrites the values on exit.
        for j in range(9):
            if j < 3:
                bc_copy_ptr[j] = bc_pointer[ytype_i * 3 + j]
            surface_matrix_copy_ptr[j] = surface_matrix_ptr[j]

        # Solve the linear equation
        zgesv(
            mat_size_ptr,              # (Input)
            bc_columns_ptr,            # (Input)
            surface_matrix_copy_ptr,   # A; (Input & Output)
            mat_size_ptr,              # (Input)
            lapack_ipiv_ptr,           # (Output)
            bc_copy_ptr,               # B -> X (Input & Output)
            mat_size_ptr,              # (Input)
            bc_solution_info_ptr       # (Ouput)
            )

        # Check for errors
        if bc_solution_info_ptr[0] != 0:
            sprintf(solution_storage_ptr.message_ptr, "RadialSolver.PropMatrixMethod:: Error encountered while applying surface boundary condition. ZGESV code: %d \nThe solutions may not be valid at the surface.\n", bc_solution_info)
            solution_storage_ptr.error_code = -21
            if verbose or raise_on_fail:
                printf(solution_storage_ptr.message_ptr)
            if raise_on_fail:
                exit(EXIT_FAILURE)
        
        # Step through each radial step and apply the propagation matrix to the surface solution
        for slice_i in range(total_slices):
            index_shift_18        = slice_i * 18
            ytype_shift           = ytype_i * MAX_NUM_Y
            solution_slice_ishift = num_bc_models * slice_i
            full_shift            = ytype_shift + solution_slice_ishift

            # Perform matrix multiplication: prop_matrix @ surface_solution 
            for j in range(6):
                temp_cmplx = cf_build_dblcmplx(0., 0.)
                for jj in range(3):
                    temp_cmplx += (
                        propagation_mtx_ptr[index_shift_18 + j * 6 + jj] * 
                        bc_copy_ptr[jj]  # Recall that bc_copy_ptr now contains the solution to the A * X = B linear system
                        )
                solution_ptr[full_shift + j] = temp_cmplx
            
            # As discussed in B13 (discussed near their equation 7), SVC16 (and the earlier 2004 book) use a different
            #    convention for tidal_y than is used by Takeuchi and Saito (1972). Since a good chunk of the field follows the
            #    latter, we will do the same. Below are the conversions from SVC16 to TS72
            ts_conversion_ptr[0] = solution_ptr[full_shift + 0]        # No Change
            ts_conversion_ptr[1] = solution_ptr[full_shift + 2]        # Flip y3 for y2
            ts_conversion_ptr[2] = solution_ptr[full_shift + 1]        # Flip y2 for y3
            ts_conversion_ptr[3] = solution_ptr[full_shift + 3]        # No Change
            ts_conversion_ptr[4] = -1. * solution_ptr[full_shift + 4]  # Change sign
            ts_conversion_ptr[5] = -1. * solution_ptr[full_shift + 5]  # Change sign

            # Store converted values back into solution pointer
            for i in range(6):
                solution_ptr[full_shift + i] = ts_conversion_ptr[i]
        
        # Get ready for next y-type solution
        ytype_i += 1

    if solution_storage_ptr.error_code == 0:
        # Update status message
        solution_storage_ptr.success = True
        strcpy(solution_storage_ptr.message_ptr, 'RadialSolver (propagation matrix method) completed without any noted issues.\n')

    # Free memory
    free_mem(fundamental_mtx_ptr)
    free_mem(inverse_fundamental_mtx_ptr)
    free_mem(derivative_mtx_ptr)
    free_mem(propagation_mtx_ptr)
