

cdef void kamata_solid_dynamic_compressible(
    double frequency,
    double radius,
    double density,
    double bulk_modulus,
    double complex shear_modulus,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* initial_conditions
    ) noexcept nogil

cdef void kamata_solid_static_compressible(
        double radius,
        double density,
        double bulk_modulus,
        double complex shear_modulus,
        unsigned int degree_l,
        double G_to_use,
        ssize_t num_ys, 
        double complex* initial_conditions
        ) noexcept nogil

########################################################################################################################
#### Liquid Layers
########################################################################################################################

cdef void kamata_liquid_dynamic_compressible(
    double frequency,
    double radius,
    double density,
    double bulk_modulus,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* initial_conditions
    ) noexcept nogil

cdef void kamata_liquid_dynamic_incompressible(
    double frequency,
    double radius,
    double density,
    unsigned int degree_l,
    double G_to_use,
    ssize_t num_ys, 
    double complex* initial_conditions
    ) noexcept nogil