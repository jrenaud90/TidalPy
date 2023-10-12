from libcpp cimport bool as bool_cpp_t


cdef void cf_solve_upper_y_at_interface(
        double complex* lower_layer_y_ptr,
        double complex* upper_layer_y_ptr,
        size_t num_sols_lower,
        size_t num_sols_upper,
        size_t num_y_lower,
        size_t num_y_upper,
        bool_cpp_t lower_is_solid,
        bool_cpp_t lower_is_static,
        bool_cpp_t lower_is_incompressible,
        bool_cpp_t upper_is_solid,
        bool_cpp_t upper_is_static,
        bool_cpp_t upper_is_incompressible,
        double interface_gravity,
        double liquid_density,
        double G_to_use = *
        ) noexcept nogil
