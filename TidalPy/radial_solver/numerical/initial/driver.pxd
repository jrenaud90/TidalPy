from libcpp cimport bool as bool_cpp_t


cdef void cf_find_initial_conditions(
    bool_cpp_t is_solid,
    bool_cpp_t is_static,
    bool_cpp_t is_incompressible,
    bool_cpp_t use_kamata,
    double frequency,
    double radius,
    double density,
    double bulk_modulus,
    double complex shear_modulus,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* initial_conditions_ptr
    )