""" Cython wrappers for `radial_derivatives_dynamic.pyx` """

from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_dynamic_incomp_x cimport dy_solid_dynamic_incompressible_x, dy_liquid_dynamic_incompressible_x

cdef double G
G = 6.67430e-11
cpdef void dy_solid_dynamic_incompressible_wrap(
        double radius, double[:] radial_functions, double[:] dy,
        double complex shear_modulus, double density, double gravity, double frequency,
        unsigned int degree_l=2, double G_to_use=G):
    """ Test wrapper for the radial derivatives (for solid layers) function. """

    dy_solid_dynamic_incompressible_x(
        radius, radial_functions, dy,
        shear_modulus, density, gravity, frequency,
        degree_l=degree_l, G_to_use=G_to_use)


cpdef void dy_liquid_dynamic_incompressible_wrap(
        double radius, double[:] radial_functions, double[:] dy,
        double density, double gravity, double frequency,
        unsigned int degree_l=2, double G_to_use=G):
    """ Test wrapper for the radial derivatives (for liquid layers) function. """

    dy_liquid_dynamic_incompressible_x(
        radius, radial_functions, dy,
        density, gravity, frequency,
        degree_l=degree_l, G_to_use=G_to_use)
