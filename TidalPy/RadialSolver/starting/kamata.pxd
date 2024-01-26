

cdef void cf_kamata_solid_dynamic_compressible(
    double frequency,
    double radius,
    double density,
    double bulk_modulus,
    double complex shear_modulus,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* starting_conditions_ptr
    ) noexcept nogil

cdef void cf_kamata_solid_static_compressible(
    double radius,
    double density,
    double bulk_modulus,
    double complex shear_modulus,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* starting_conditions_ptr
    ) noexcept nogil

cdef void cf_kamata_solid_dynamic_incompressible(
    double frequency,
    double radius,
    double density,
    double complex shear_modulus,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* starting_conditions_ptr
    ) noexcept nogil

########################################################################################################################
#### Liquid Layers
########################################################################################################################

cdef void cf_kamata_liquid_dynamic_compressible(
    double frequency,
    double radius,
    double density,
    double bulk_modulus,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* starting_conditions_ptr
    ) noexcept nogil

cdef void cf_kamata_liquid_dynamic_incompressible(
    double frequency,
    double radius,
    double density,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* starting_conditions_ptr
    ) noexcept nogil