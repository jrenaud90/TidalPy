cdef void dy_solid_static_compressible_x(double radius, double[:] radial_functions, double[:] dy, double complex shear_modulus, double bulk_modulus, double density, double gravity, unsigned int degree_l = *, double G_to_use = *) nogil

cdef void dy_liquid_static_compressible_x(double radius, double[:] radial_functions, double[:] dy, double density, double gravity, unsigned int degree_l = *, double G_to_use = *) nogil
