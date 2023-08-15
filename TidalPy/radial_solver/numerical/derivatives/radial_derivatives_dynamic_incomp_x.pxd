cdef void dy_solid_dynamic_incompressible_x(double radius, double[:] radial_functions, double[:] dy, double complex shear_modulus, double density, double gravity, double frequency, unsigned int degree_l = *, double G_to_use = *) nogil

cdef void dy_liquid_dynamic_incompressible_x(double radius, double[:] radial_functions, double[:] dy, double density, double gravity, double frequency, unsigned int degree_l = *, double G_to_use = *) nogil
