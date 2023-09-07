from libcpp cimport bool as bool_cpp_t

cpdef Py_ssize_t find_solution_num(
        bool_cpp_t is_solid,
        bool_cpp_t is_static,
        bool_cpp_t is_incompressible) noexcept nogil

cpdef void interface_x(
        double complex[:, :] lower_layer_y,
        double complex[:, :] upper_layer_y,
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