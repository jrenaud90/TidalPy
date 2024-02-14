cdef void cf_solve_upper_y_at_interface(
        double complex* lower_layer_y_ptr,
        double complex* upper_layer_y_ptr,
        size_t num_sols_lower,
        size_t num_sols_upper,
        size_t num_y_lower,
        size_t num_y_upper,
        int lower_layer_type,
        bint lower_is_static,
        bint lower_is_incompressible,
        int upper_layer_type,
        bint upper_is_static,
        bint upper_is_incompressible,
        double interface_gravity,
        double liquid_density,
        double G_to_use = *
        ) noexcept nogil
