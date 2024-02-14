# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from scipy.linalg.cython_lapack cimport zgesv
from libc.math cimport NAN, pi

cdef void cf_apply_surface_bc(
        double complex* constant_vector_ptr,
        int* bc_solution_info_ptr,
        double* bc_pointer,
        double complex* uppermost_y_per_solution_ptr,
        double surface_gravity,
        double G_to_use,
        unsigned char num_sols,
        unsigned char max_num_y,
        unsigned char ytype_i,
        int layer_type,
        bint layer_is_static,
        bint layer_is_incomp
        ) noexcept nogil:
    
    # Variables used to solve the linear equation at the planet's surface.
    # These are passed to the LAPACK solver.
    # NRHS = number of solutions that will be solved at the same time. Only one will be solved per radial_solver call.
    cdef int lapack_nrhs = 1
    cdef int* lapack_nrhs_ptr = &lapack_nrhs

    # IPIV = Integer pivot array that is an additional output provided by ZGESV. It is not used but must be provided.
    #  It must be at least as large as the largest dimension of the input matrix, for this work that is 3.
    cdef int[10] lapack_ipiv
    cdef int* lapack_ipiv_ptr = &lapack_ipiv[0]

    cdef int num_sols_int = <int>num_sols
    cdef int* num_sols_int_ptr = &num_sols_int

    # Allocate surface matrices. We only need one of these (which one depends on the uppermost layer type).
    # But, since there are only 3 and all of them are small, we will just allocate all of them separately on the stack.
    cdef double complex[3][3] surface_matrix_solid
    cdef double complex[2][2] surface_matrix_liquid_dynamic
    cdef double complex[1][1] surface_matrix_liquid_static
    cdef double complex* surface_matrix_ptr

    # Create coefficient matrix based on surface layer type.
    if (layer_type == 0):
        # Solid layer
        # Set pointer to correct matrix
        surface_matrix_ptr = &surface_matrix_solid[0][0]

        # At the surface: y_2 = S_1; y_4 = S_4; y_6 = S_6 [See: B.37 in KTC21; 16 in KMN15]
        # We will set the constant vector equal to the surface boundary condition.
        #  It will be overwritten with the solution to the linear solution.
        constant_vector_ptr[0] = bc_pointer[ytype_i * 3 + 0]
        constant_vector_ptr[1] = bc_pointer[ytype_i * 3 + 1]
        constant_vector_ptr[2] = bc_pointer[ytype_i * 3 + 2]

        # The definitions above need to be transposed as the LAPACK solver TidalPy uses requires
        #  FORTRAN-ordered arrays.
        surface_matrix_ptr[0] = uppermost_y_per_solution_ptr[0 * max_num_y + 1]
        surface_matrix_ptr[1] = uppermost_y_per_solution_ptr[0 * max_num_y + 3]
        surface_matrix_ptr[2] = uppermost_y_per_solution_ptr[0 * max_num_y + 5]
        surface_matrix_ptr[3] = uppermost_y_per_solution_ptr[1 * max_num_y + 1]
        surface_matrix_ptr[4] = uppermost_y_per_solution_ptr[1 * max_num_y + 3]
        surface_matrix_ptr[5] = uppermost_y_per_solution_ptr[1 * max_num_y + 5]
        surface_matrix_ptr[6] = uppermost_y_per_solution_ptr[2 * max_num_y + 1]
        surface_matrix_ptr[7] = uppermost_y_per_solution_ptr[2 * max_num_y + 3]
        surface_matrix_ptr[8] = uppermost_y_per_solution_ptr[2 * max_num_y + 5]

    else:
        if layer_is_static:
            # Set pointer to correct matrix
            surface_matrix_ptr = &surface_matrix_liquid_dynamic[0][0]

            # Unlike the dynamic liquid layer, a static liquid layer's y_2 is undefined. That leads to one less boundary condition
            #  and one less solution (1 total).
            #  At the surface, y_7 = S_7 [See: Eq. 17, 10 in S74]

            # We will set the constant vector equal to the surface boundary condition.
            #  It will be overwritten with the solution to the linear solution.
            # y_7 = y_6 + (4 pi G / g) y_2
            constant_vector_ptr[0] = \
                bc_pointer[ytype_i * 3 + 2] + \
                bc_pointer[ytype_i * 3 + 0] * (4. * pi * G_to_use / surface_gravity)

            # These are unused. Set to NAN so if they do get used we might be able to catch it.
            constant_vector_ptr[1] = NAN
            constant_vector_ptr[2] = NAN

            # Note: for a static liquid layer, y_7 held in index 1 (index 0 is y_5).
            surface_matrix_ptr[0] = uppermost_y_per_solution_ptr[0 * max_num_y + 1]
        else:
            # Set pointer to correct matrix
            surface_matrix_ptr = &surface_matrix_liquid_static[0][0]

            # Unlike the solid layer, a liquid layer's y_4 is undefined. That leads to one less boundary condition and one
            #  less solution (2 total).
            #  At the surface, y_2 = S_1; y_6 = S_6 [See: Eq. B.38 in KTC21; Eq. 17 in KMN15

            # We will set the constant vector equal to the surface boundary condition.
            #  It will be overwritten with the solution to the linear solution.
            # The surface boundary condition will still have 3 members. Drop the one related to y_4
            constant_vector_ptr[0] = bc_pointer[ytype_i * 3 + 0]
            constant_vector_ptr[1] = bc_pointer[ytype_i * 3 + 2]

            # The last constant is unused. Set to NAN so if they do get used we might be able to catch it.
            constant_vector_ptr[2] = NAN

            # Build y-solution matrix to be applied to the surface.
            # Note: for a dynamic liquid, y_2 and y_6 are held at indices 1 and 3 respectively
            surface_matrix_ptr[0] = uppermost_y_per_solution_ptr[0 * max_num_y + 1]
            surface_matrix_ptr[1] = uppermost_y_per_solution_ptr[0 * max_num_y + 3]
            surface_matrix_ptr[2] = uppermost_y_per_solution_ptr[1 * max_num_y + 1]
            surface_matrix_ptr[3] = uppermost_y_per_solution_ptr[1 * max_num_y + 3]

    # Find the solution to the linear equation at the surface surface_matrix @ constant_vector = bc
    # ZGESV computes the solution to system of linear equations A * X = B for GE matrices
    # See https://www.netlib.org/lapack/explore-html/d6/d10/group__complex16_g_esolve_ga531713dfc62bc5df387b7bb486a9deeb.html#ga531713dfc62bc5df387b7bb486a9deeb
    zgesv(
        num_sols_int_ptr,
        lapack_nrhs_ptr,
        surface_matrix_ptr,   # Modified Variable
        num_sols_int_ptr,
        lapack_ipiv_ptr,      # Modified Variable
        constant_vector_ptr,  # Modified Variable
        num_sols_int_ptr,
        bc_solution_info_ptr  # Modified Variable
        )
