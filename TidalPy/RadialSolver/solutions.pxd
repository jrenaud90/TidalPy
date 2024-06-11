cdef class RadialSolverSolution():

    cdef public bint success
    cdef char* _message
    cdef bytes _message_bytes

    # Result structure information
    cdef size_t num_ys
    cdef size_t num_slices
    cdef size_t total_size
    cdef size_t num_ytypes
    cdef char* ytypes[5]
    cdef bint ytypes_set

    # Result pointers and data
    cdef double complex* full_solution_ptr
    cdef double complex[::1] full_solution_view

    # Love number information
    cdef double complex* complex_love_ptr
    cdef double complex[::1] complex_love_view

    cdef void set_models(self, int* bc_models_ptr) noexcept nogil
    cdef void set_message(self, str new_message) noexcept

cdef size_t cf_find_num_shooting_solutions(
    int layer_type,
    bint is_static,
    bint is_incompressible
    ) noexcept nogil
