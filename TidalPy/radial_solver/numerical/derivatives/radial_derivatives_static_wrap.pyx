""" Cython wrappers for `radial_derivatives_static.pyx` """

from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_static_x cimport dy_solid_static_compressible_x, dy_liquid_static_compressible_x

cdef double G
G = 6.67430e-11
cpdef void dy_solid_static_compressible_wrap(
        double radius, double[:] radial_functions, double[:] dy,
        double complex shear_modulus, double bulk_modulus, double density, double gravity,
        unsigned int degree_l=2, double G_to_use=G):
    """ Test wrapper for the radial derivatives (for solid layers) function. """

    dy_solid_static_compressible_x(
        radius, radial_functions, dy,
        shear_modulus, bulk_modulus, density, gravity,
        degree_l=degree_l, G_to_use=G_to_use)


cpdef void dy_liquid_static_compressible_wrap(
        double radius, double[:] radial_functions, double[:] dy,
        double density, double gravity,
        unsigned int degree_l=2, double G_to_use=G):
    """ Test wrapper for the radial derivatives (for liquid layers) function. """

    dy_liquid_static_compressible_x(
        radius, radial_functions, dy,
        density, gravity,
        degree_l=degree_l, G_to_use=G_to_use)
