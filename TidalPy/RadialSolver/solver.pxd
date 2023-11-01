from libcpp cimport bool as bool_cpp_t

cdef class RadialSolverSolution():

    cdef public str message
    cdef public bool_cpp_t success

    cdef double complex* full_solution_ptr
    cdef double complex[::1] full_solution_view

    cdef size_t num_ys
    cdef size_t num_slices
    cdef size_t total_size
    cdef size_t num_ytypes

    cdef tuple ytypes