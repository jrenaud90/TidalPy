from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex


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
    cpp_bool layer_is_static,
    cpp_bool layer_above_is_static,
    cpp_bool layer_is_incomp,
    cpp_bool layer_above_is_incomp,
    size_t num_sols,
    size_t max_num_y
    ) noexcept nogil
