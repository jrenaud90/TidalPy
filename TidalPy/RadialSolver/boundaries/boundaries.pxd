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
    int layer_type,
    bint layer_is_static,
    bint layer_is_incomp
    ) noexcept nogil
