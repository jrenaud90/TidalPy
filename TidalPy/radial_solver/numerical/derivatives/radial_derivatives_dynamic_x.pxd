cdef void radial_derivatives_solid_general_x(double radius, double[:] radial_functions, double[:] dy, double complex shear_modulus, double bulk_modulus, double density, double gravity, double frequency, unsigned int degree_l = *, double G_to_use = *) nogil

cdef void radial_derivatives_liquid_general_x(double radius, double[:] radial_functions, double[:] dy, double bulk_modulus, double density, double gravity, double frequency, unsigned int degree_l = *, double G_to_use = *) nogil
