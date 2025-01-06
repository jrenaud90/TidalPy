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

from libc.stdio cimport printf, sprintf
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.string cimport strcpy
from libc.math cimport fmin
from libcpp cimport bool as cpp_bool

import numpy as np
cimport numpy as cnp
cnp.import_array()
from scipy.linalg.cython_lapack cimport zgesv

from TidalPy.constants cimport d_PI_DBL
from TidalPy.utilities.math.complex cimport cmplx_one, cmplx_zero, cmplx_NAN, cf_build_dblcmplx
from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC
from TidalPy.RadialSolver.rs_solution cimport RadialSolutionStorageCC
from TidalPy.RadialSolver.boundaries.surface_bc cimport cf_get_surface_bc
from TidalPy.RadialSolver.matrix_types.solid_matrix cimport cf_fundamental_matrix
from TidalPy.RadialSolver.constants cimport MAX_NUM_Y

ctypedef double complex double_complex_t


cdef int cf_matrix_propagate(
        RadialSolutionStorageCC* solution_storage_ptr,
        double frequency,
        double planet_bulk_density,
        # TODO: In the future the propagation matrix should take in layer types and multiple layers
        # int* layer_types_ptr,
        # int* is_static_by_layer_ptr,
        # int* is_incompressible_by_layer_ptr,
        vector[size_t] first_slice_index_by_layer_vec,
        vector[size_t] num_slices_by_layer_vec,
        size_t num_bc_models,
        int* bc_models_ptr,
        double G_to_use,
        int degree_l,
        double starting_radius,
        double start_radius_tolerance,
        int core_model,
        cpp_bool verbose
        ) noexcept nogil:

    # Setup
    cdef size_t i, j, k, jj, ytype_i, slice_i
    cdef size_t last_index_shift_36, index_shift_36, row_shift_index, last_index_shift_18, index_shift_18, index_shift_max_y, full_shift

    # Get raw pointer of radial solver storage and eos storage
    cdef EOSSolutionCC* eos_solution_storage_ptr = solution_storage_ptr.get_eos_solution_ptr()

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
    cdef size_t num_output_ys     = MAX_NUM_Y * num_ytypes

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

    # Make conversions from TS72 to SV04
    # The propagation matrix used by Sabadini & Veermerson (2004) uses different sign conventions so we need to convert
    for ytype_i in range(num_ytypes):
        full_shift = 3 * ytype_i
        if bc_models_ptr[ytype_i] == 0:
            # Free surface, as far as I know no conversion needed.
            pass
        elif bc_models_ptr[ytype_i] == 1:
            # Tidal boundary condition - the last component is negative of the TS74 (see Eq. 1.127 of S&V2004)
            bc_pointer[full_shift + 2] *= -1.0
        elif bc_models_ptr[ytype_i] == 2:
            # Loading boundary condition.
            # There is some confusion here:
            #   Sabadini+ (1982) has this as (y3, y4, y6) = [-g * (2l + 1) / (4 * pi * R^2), 0, G * (2l + 1) / R^2]  ** I think there is a typo here or some other sign convention was used in this paper for y6, it is negative in all of Sabadini's other work, even that which cites this paper.
            #   S&V (2004) has this as (y3, y4, y6) = [-g * (2l + 1) / (4 * pi * R^2), 0, -G * (2l + 1) / R^2]
            #   SVC (2016) has this as (y3, y4, y6) = [-g * (2l + 1) / (4 * pi * R^2), 0, -G * (2l + 1) / R^2]
            #  ^ These all match Farrell (1972; right after Eq. 31)
            # This is compared to H. Martens thesis which is:
            #   [-g^2 * (2l + 1) / (4 * pi * G), 0, (2l + 1) * g]
            # ^ These have different dimensions "Love numbers of [martens] vary by factor of R * g, so the surface conditions do as well" (paraphrased from Martens thesis)
            # These differ from what is used by the shooting method solver which is based on Saito (1974)
            #   [-(1 / 3) (2l + 1) rho, 0, (2l + 1) / R]
            # If we divide Martens conditions by R g we get (After some Algebra)
            #   [-rho * (2l + 1) / 3, 0, (2l + 1) / R]  -- which matches Saito.
            # As for Sabadini BCs:
            #   We know to convert from TS to SV we have to negative y6 and flip y2 for y3. So the SV BC's listed above
            #   Are really (y2, y4, -y6). Meaning that they should match the output used for the shooting method BC
            #   Except we need to take the negative of y6 as we did for the tidal case.
            # ALMA uses [-g * (2l + 1) / (4 * pi * R^2), 0, -G * (2l + 1) / R^2]
            # If I perform a 1:1 match then the SV method is a factor of G/R different for y2 condition and -G/R for y6
            # TODO: Check if this is the correct conversion. the y6 one does not make sense to me; why is the scale different than the tidal case?
            # bc_pointer[full_shift + 0] *= (G_to_use / planet_radius)
            # bc_pointer[full_shift + 1] *= 1.0  # No scale for y4
            # bc_pointer[full_shift + 2] *= (-G_to_use / planet_radius)
            
            # Okay that did not work. Simply switching the y6 condition to negative reproduces the shooting method results. _shrug_
            bc_pointer[full_shift + 2] *= -1.0

    # Determine starting radius slice
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

    # Determine which layer this starting radius resides in. We will skip the lower layers
    cdef double layer_upper_radius, last_layer_upper_radius, radius_check
    cdef size_t current_layer_i
    cdef size_t start_layer_i           = 0
    cdef size_t last_index_before_start = 0
    cdef size_t start_index_in_layer    = 0
    cdef size_t first_slice_index       = 0
    cdef double last_radius_check = 0.0
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
    
    # For the propagation matrix we have to start at index 2. 
    # The reason is that if we start at index 0 then the starting condition would be at index -1 which wont work with our setup.
    # If we start at index 1 then the starting condition is okay but we need the fundamental matrix at index 0 (start - 1)
    # But index 0 corresponds to r = 0 where many elements of the fundamental matrix are nan or inf. 
    # We sacrifice a bit of accuracy by starting slightly higher. If there are many slices in the planet then this error is small.
    first_slice_index = last_index_before_start + 1
    if first_slice_index == 0 or first_slice_index == 1:
        first_slice_index = 2

    # Define memory for our fundamental matrices.
    # These have to be heap allocated because we do not know the number of radial slices at compile time (and it could be large)
    cdef size_t matrix_size = 6 * 6 * total_slices
    # Propagation matrix has 6 rows but only 3 columns.
    cdef size_t prop_mat_size = 6 * 3 * total_slices

    cdef vector[double_complex_t] fundamental_mtx_vec = vector[double_complex_t]()
    fundamental_mtx_vec.resize(matrix_size)
    cdef vector[double_complex_t] inverse_fundamental_mtx_vec = vector[double_complex_t]()
    inverse_fundamental_mtx_vec.resize(matrix_size)
    cdef vector[double_complex_t] derivative_mtx_vec = vector[double_complex_t]()
    derivative_mtx_vec.resize(matrix_size)
    cdef vector[double_complex_t] propagation_mtx_vec = vector[double_complex_t]()
    propagation_mtx_vec.resize(prop_mat_size)

    cdef double complex* fundamental_mtx_ptr         = &fundamental_mtx_vec[0]
    cdef double complex* inverse_fundamental_mtx_ptr = &inverse_fundamental_mtx_vec[0]
    cdef double complex* derivative_mtx_ptr          = &derivative_mtx_vec[0]
    cdef double complex* propagation_mtx_ptr         = &propagation_mtx_vec[0]

    # Only populate matrix values that are used, starting at first slice index - 1 (because we need the matrix just before the solution starts)
    # Populate matrices with the correct layer type. 
    # TODO: Currently only solid, static, incompressible layers are supported for matrix propagation.
    cf_fundamental_matrix(
        <Py_ssize_t>first_slice_index - 1,
        <Py_ssize_t>total_slices,
        radius_array_ptr,
        density_array_ptr,
        gravity_array_ptr,
        complex_shear_array_ptr,
        fundamental_mtx_ptr,          # Changed variable
        inverse_fundamental_mtx_ptr,  # Changed variable
        derivative_mtx_ptr,           # Changed variable
        degree_l,
        G_to_use
        )

    # Initialize the base of the propagation matrix to the initial conditions
    ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
    # ^JPR: Confirmed via testing that the solutions quickly converge away from these starting values if everything
    #    else is well behaved. There are differences near the core which could matter for heat generation in that region.
    #    Found that the core_model=0 approach matches the shooting method the best.
    index_shift_18 = (first_slice_index - 1) * 18
    index_shift_36 = (first_slice_index - 1) * 36
    cdef double grav_constant
    if core_model == 0:
        # Henning & Hurford (2014): "At the core, a special seed matrix Bcore is created with only three columns,
        # equal to the first, second, and third columns of Y for the properties at the base layer."
        for j in range(6):
            row_shift_index = index_shift_18 + (j * 3)
            for k in range(3):
                propagation_mtx_ptr[row_shift_index + k] = fundamental_mtx_ptr[index_shift_36 + j * 6 + k]
    elif core_model == 1:
        # Roberts & Nimmo (2008): liquid innermost zone.
        for j in range(6):
            row_shift_index = index_shift_18 + (j * 3)
            for k in range(3):
                if (j == 2) and (k == 0):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                elif (j == 3) and (k == 1):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                elif (j == 5) and (k == 2):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                else:
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_zero
    elif core_model == 2:
        # Solid Inner Core (Based on Henning & Hurford 2014)
        for j in range(6):
            row_shift_index = index_shift_18 + (j * 3)
            for k in range(3):
                if (j == 0) and (k == 0):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                elif (j == 1) and (k == 1):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                elif (j == 2) and (k == 2):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                else:
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_zero
    elif core_model == 3:
        # Liquid Inner Core (based on Tobie+2005; As determined by Marc Neveu for IcyDwarf).
        for j in range(6):
            row_shift_index = index_shift_18 + (j * 3)
            for k in range(3):
                if (j == 0) and (k == 0):
                    propagation_mtx_ptr[row_shift_index + k] = cf_build_dblcmplx(0.05, 0.)
                elif (j == 1) and (k == 1):
                    propagation_mtx_ptr[row_shift_index + k] = cf_build_dblcmplx(0.01, 0.)
                elif (j == 5) and (k == 2):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                else:
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_zero
    elif core_model == 4:
        # Interface matrix from SVC Eq. 1.150
        grav_constant = (4. / 3.) * d_PI_DBL * G_to_use * density_array_ptr[first_slice_index - 1]
        for j in range(6):
            row_shift_index = index_shift_18 + (j * 3)
            for k in range(3):
                if (j == 0) and (k == 0):
                    propagation_mtx_ptr[row_shift_index + k] = -radius_array_ptr[first_slice_index - 1]**(degree_l - 1) / grav_constant
                elif (j == 0) and (k == 2):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                elif (j == 1) and (k == 1):
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_one
                elif (j == 2) and (k == 2):
                    propagation_mtx_ptr[row_shift_index + k] = density_array_ptr[first_slice_index - 1] * grav_constant * radius_array_ptr[first_slice_index - 1]
                elif (j == 4) and (k == 0):
                    propagation_mtx_ptr[row_shift_index + k] = radius_array_ptr[first_slice_index - 1]**degree_l
                elif (j == 5) and (k == 0):
                    propagation_mtx_ptr[row_shift_index + k] = 2.0 * (degree_l - 1) * radius_array_ptr[first_slice_index - 1]**(degree_l - 1)
                elif (j == 5) and (k == 2):
                    propagation_mtx_ptr[row_shift_index + k] = 3.0 * grav_constant
                else:
                    propagation_mtx_ptr[row_shift_index + k] = cmplx_zero
    else:
        sprintf(solution_storage_ptr.message_ptr, "RadialSolver.PropMatrixMethod:: Unknown starting core conditions encountered in `cf_matrix_propagate`: %d (acceptable values: 0, 1, 2, 3)\n", core_model)
        solution_storage_ptr.error_code = -20
        solution_storage_ptr.success = False
        if verbose:
            printf(solution_storage_ptr.message_ptr)
        return solution_storage_ptr.error_code

    # Step through the planet's shells and build the propagation matrix
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

    for slice_i in range(first_slice_index, total_slices):
        
        # Need to start the index for this radial slice. The shift is based on the size of the respective matrix.
        # Fundamental matrix is 6x6
        index_shift_36 = slice_i * 36
        last_index_shift_36 = (slice_i - 1) * 36
        # Propagation matrix is 6x3
        index_shift_18 = slice_i * 18
        last_index_shift_18 = (slice_i - 1) * 18

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
                        propagation_mtx_ptr[last_index_shift_18 + jj * 3 + k]
                        )
                temp_matrix_ptr[j * 3 + k] = temp_cmplx

        # Now perform the outer matrix multiplication
        # P_{i} = Y_{i} @ A
        for j in range(6):
            for k in range(3):
                temp_cmplx = cf_build_dblcmplx(0., 0.)
                for jj in range(6):
                    temp_cmplx += (
                        fundamental_mtx_ptr[index_shift_36 + j * 6 + jj] * 
                        temp_matrix_ptr[jj * 3 + k]
                        )
                propagation_mtx_ptr[index_shift_18 + j * 3 + k] = temp_cmplx

        # We need to define a matrix that equals the propagation matrix at the surface value.
        # Surface condition matrix is a 3x3 matrix of the top-most shell of the propagation matrix's rows [3, 4, 6]
        if slice_i == (total_slices - 1):
            for i in range(3):
                surface_matrix_ptr[0 + i] = propagation_mtx_ptr[index_shift_18 + (2 * 3) + i]
                surface_matrix_ptr[3 + i] = propagation_mtx_ptr[index_shift_18 + (3 * 3) + i]
                surface_matrix_ptr[6 + i] = propagation_mtx_ptr[index_shift_18 + (5 * 3) + i]

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
    cdef double* solution_dbl_ptr = solution_storage_ptr.full_solution_vec.data()
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
        # While we make a copy, use this chance to take the transpose of the surface matrix. This is required because
        # the zgesv solver requires FORTRAN-ordered arrays, whereas we have been using C-ordered arrays
        for i in range(3):
            bc_copy_ptr[i] = bc_pointer[ytype_i * 3 + i]
            for j in range(3):
                # Make a copy of transpose
                surface_matrix_copy_ptr[3 * j + i] = surface_matrix_ptr[3 * i + j]

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
            solution_storage_ptr.success = False
            if verbose:
                printf(solution_storage_ptr.message_ptr)
            return solution_storage_ptr.error_code
        
        # Step through each radial step and apply the propagation matrix to the surface solution
        for slice_i in range(total_slices):
            index_shift_18        = slice_i * 18
            ytype_shift           = ytype_i * MAX_NUM_Y
            full_shift            = num_output_ys * slice_i + ytype_shift

            if slice_i < first_slice_index:
                # We are not in the solution region yet. Just nan the results
                for i in range(6):
                    solution_ptr[full_shift + i] = cmplx_NAN
            else:
                # Perform matrix multiplication: prop_matrix @ surface_solution 
                for j in range(6):
                    temp_cmplx = cf_build_dblcmplx(0., 0.)
                    for jj in range(3):
                        temp_cmplx += (
                            propagation_mtx_ptr[index_shift_18 + j * 3 + jj] * 
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

    # Update solution status and return
    if solution_storage_ptr.error_code != 0:
        solution_storage_ptr.success = False
    else:
        solution_storage_ptr.success = True
        solution_storage_ptr.set_message('RadialSolver.MatrixPropagation: Completed without any noted issues.')
    
    return solution_storage_ptr.error_code
