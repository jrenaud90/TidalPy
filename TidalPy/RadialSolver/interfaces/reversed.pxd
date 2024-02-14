cdef void cf_top_to_bottom_interface_bc(
    double complex* constant_vector_ptr,
    double complex* layer_above_constant_vector_ptr,
    double complex* uppermost_y_per_solution_ptr,
    double gravity_upper,
    double layer_above_lower_gravity,
    double density_upper,
    double layer_above_lower_density,
    int layer_type,
    int layer_above_type,
    bint layer_is_static,
    bint layer_above_is_static,
    bint layer_is_incomp,
    bint layer_above_is_incomp,
    unsigned char num_sols,
    unsigned char max_num_y
    ) noexcept nogil
