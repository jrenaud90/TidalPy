cdef extern from "interp_.hpp" nogil:
    size_t c_binary_search_with_guess(
        double key,
        double* array,
        size_t length,
        size_t guess
        )

    void c_interp(
        double* desired_x_ptr,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        size_t len_x,
        size_t* provided_j_ptr,
        double* result_ptr
        )

    void c_interp_complex(
        double desired_x,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        size_t len_x,
        size_t* provided_j_ptr,
        double* result_ptr
        )
