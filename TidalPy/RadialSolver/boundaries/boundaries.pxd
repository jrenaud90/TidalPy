from libcpp cimport bool as cpp_bool


cdef void cf_apply_surface_bc(
    double complex* constant_vector_ptr,
    int* bc_solution_info,
    double* bc_pointer,
    double complex* uppermost_y_per_solution_ptr,
    double surface_gravity,
    double G_to_use,
    size_t num_sols,
    size_t max_num_y,
    size_t ytype_i,
    int layer_type,
    cpp_bool layer_is_static,
    cpp_bool layer_is_incomp
    ) noexcept nogil
