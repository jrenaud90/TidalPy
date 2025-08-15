from libcpp cimport bool as cpp_bool
from libcpp.string cimport string as cpp_string


cdef void cf_find_starting_conditions(
    cpp_bool* success_ptr,
    cpp_string& message,
    int layer_type,
    bint is_static,
    bint is_incompressible,
    cpp_bool use_kamata,
    double frequency,
    double radius,
    double density,
    double complex bulk_modulus,
    double complex shear_modulus,
    int degree_l,
    double G_to_use,
    size_t num_ys, 
    double complex* starting_conditions_ptr,
    cpp_bool run_y_checks = *
    ) noexcept nogil
