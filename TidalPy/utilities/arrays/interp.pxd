cdef extern from "interp_.cpp" nogil:
    size_t cf_binary_search_with_guess(
        double key,
        double* array,
        size_t length,
        size_t guess
        )

    void cf_interp(
        double* desired_x_ptr,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        size_t len_x,
        size_t* provided_j_ptr,
        double* result_ptr
        )

    void cf_interp_complex(
        double desired_x,
        double* x_domain_ptr,
        double* dependent_values_ptr,
        size_t len_x,
        size_t* provided_j_ptr,
        double* result_ptr
        )
