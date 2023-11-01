from libcpp cimport bool as bool_cpp_t

cdef void cf_apply_surface_bc(
    double complex* constant_vector_ptr,
    int* bc_solution_info,
    double* bc_pointer,
    double complex* uppermost_y_per_solution_ptr,
    double surface_gravity,
    double G_to_use,
    unsigned char num_sols,
    unsigned char max_num_y,
    unsigned char ytype_i,
    bool_cpp_t layer_is_solid,
    bool_cpp_t layer_is_static,
    bool_cpp_t layer_is_incomp
    ) noexcept nogil