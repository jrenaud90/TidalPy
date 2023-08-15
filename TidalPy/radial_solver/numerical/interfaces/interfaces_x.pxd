from libcpp cimport bool as bool_cpp_t

cdef unsigned int solution_num(bool_cpp_t is_solid, bool_cpp_t is_static, bool_cpp_t is_compressible) nogil

cdef void interface_x(
        double complex[:, :] lower_layer_y, double complex[:, :] upper_layer_y,
        bool_cpp_t lower_is_solid, bool_cpp_t lower_is_static, bool_cpp_t lower_is_compressible,
        bool_cpp_t upper_is_solid, bool_cpp_t upper_is_static, bool_cpp_t upper_is_compressible,
        double interface_gravity, double liquid_density, double G_to_use = *
        ) nogil