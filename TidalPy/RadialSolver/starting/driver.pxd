from libcpp cimport bool as bool_cpp_t


cdef void cf_find_starting_conditions(
    int layer_type,
    bint is_static,
    bint is_incompressible,
    bool_cpp_t use_kamata,
    double frequency,
    double radius,
    double density,
    double bulk_modulus,
    double complex shear_modulus,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* starting_conditions_ptr,
    bool_cpp_t run_y_checks = *
    )