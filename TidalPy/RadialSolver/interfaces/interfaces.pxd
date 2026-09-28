from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex


cdef extern from "interfaces_.hpp" nogil:

    cdef void c_solve_upper_y_at_interface(
        cpp_complex[double]* lower_layer_y_ptr,
        cpp_complex[double]* upper_layer_y_ptr,
        size_t num_sols_lower,
        size_t num_sols_upper,
        size_t max_num_y,
        int lower_layer_type,
        cpp_bool lower_is_static,
        cpp_bool lower_is_incompressible,
        int upper_layer_type,
        cpp_bool upper_is_static,
        cpp_bool upper_is_incompressible,
        double interface_gravity,
        double liquid_density,
        double G_to_use) noexcept nogil
