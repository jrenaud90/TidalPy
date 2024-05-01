cdef class RadialSolverSolution():

    cdef public bint success
    cdef char* _message

    # Result structure information
    cdef size_t num_ys
    cdef size_t num_slices
    cdef size_t total_size
    cdef size_t num_ytypes
    cdef char* ytypes[5]

    # Result pointers and data
    cdef double complex* full_solution_ptr
    cdef double complex[::1] full_solution_view

    # Love number information
    cdef double complex* complex_love_ptr
    cdef double complex[::1] complex_love_view
