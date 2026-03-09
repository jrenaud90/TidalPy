from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.string cimport string as cpp_string


cdef extern from "driver_.hpp" nogil:

    cdef void c_find_starting_conditions(
        cpp_bool* success_ptr,
        cpp_string& message,
        int layer_type,
        int is_static,
        int is_incompressible,
        cpp_bool use_kamata,
        double frequency,
        double radius,
        double density,
        cpp_complex[double] bulk_modulus,
        cpp_complex[double] shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr,
        cpp_bool run_y_checks) noexcept nogil
