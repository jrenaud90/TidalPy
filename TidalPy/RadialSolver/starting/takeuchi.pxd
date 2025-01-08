cdef void cf_takeuchi_solid_dynamic_compressible(
    double frequency,
    double radius,
    double density,
    double complex bulk_modulus,
    double complex shear_modulus,
    int degree_l,
    double G_to_use,
    size_t num_ys, 
    double complex* starting_conditions_ptr
    ) noexcept nogil

cdef void cf_takeuchi_solid_static_compressible(
    double radius,
    double density,
    double complex bulk_modulus,
    double complex shear_modulus,
    int degree_l,
    double G_to_use,
    size_t num_ys, 
    double complex* starting_conditions_ptr
    ) noexcept nogil

########################################################################################################################
#### Liquid Layers
########################################################################################################################

cdef void cf_takeuchi_liquid_dynamic_compressible(
    double frequency,
    double radius,
    double density,
    double complex bulk_modulus,
    int degree_l,
    double G_to_use,
    size_t num_ys, 
    double complex* starting_conditions_ptr
    ) noexcept nogil