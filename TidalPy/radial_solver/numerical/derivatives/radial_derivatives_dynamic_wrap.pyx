""" Cython tests for `radial_derivatives_dynamic.pyx` """

from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_dynamic_x cimport radial_derivatives_solid_general_x, radial_derivatives_liquid_general_x

cdef double G
G = 6.67430e-11
cpdef void radial_derivatives_solid_general_wrap(
        double radius, double[:] radial_functions, double[:] dy,
        double complex shear_modulus, double bulk_modulus, double density, double gravity, double frequency,
        unsigned int degree_l=2, double G_to_use=G):
    """ Test wrapper for the radial derivatives (for solid layers) function. """

    radial_derivatives_solid_general_x(
        radius, radial_functions, dy,
        shear_modulus, bulk_modulus, density, gravity, frequency,
        degree_l=degree_l, G_to_use=G_to_use)


cpdef void radial_derivatives_liquid_general_wrap(
        double radius, double[:] radial_functions, double[:] dy,
        double bulk_modulus, double density, double gravity, double frequency,
        unsigned int degree_l=2, double G_to_use=G):
    """ Test wrapper for the radial derivatives (for liquid layers) function. """

    radial_derivatives_liquid_general_x(
        radius, radial_functions, dy,
        bulk_modulus, density, gravity, frequency,
        degree_l=degree_l, G_to_use=G_to_use)
