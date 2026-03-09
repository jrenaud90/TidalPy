# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp.complex cimport complex as cpp_complex


def takeuchi_solid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Takeuchi starting conditions for a solid dynamic compressible layer.

    TS72 Eqs. 95-102. Three independent solutions.

    Parameters
    ----------
    frequency : float
        Forcing frequency [rad s-1].
    radius : float
        Radius [m].
    density : float
        Density [kg m-3].
    bulk_modulus : complex
        Bulk modulus [Pa].
    shear_modulus : complex
        Shear modulus [Pa].
    degree_l : int
        Tidal harmonic order.
    G_to_use : float
        Gravitational constant.
    starting_conditions_view : complex[:, ::1]
        Output array of shape [num_solutions, num_ys].
    """
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    cdef cpp_complex[double] K = cpp_complex[double](bulk_modulus.real, bulk_modulus.imag)
    cdef cpp_complex[double] mu = cpp_complex[double](shear_modulus.real, shear_modulus.imag)
    c_takeuchi_solid_dynamic_compressible(frequency, radius, density, K, mu, degree_l, G_to_use, num_ys, ptr)


def takeuchi_solid_static_compressible(
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Takeuchi starting conditions for a solid static compressible layer.

    TS72 Eqs. 95-102 (w=0). Three independent solutions.
    """
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    cdef cpp_complex[double] K = cpp_complex[double](bulk_modulus.real, bulk_modulus.imag)
    cdef cpp_complex[double] mu = cpp_complex[double](shear_modulus.real, shear_modulus.imag)
    c_takeuchi_solid_static_compressible(radius, density, K, mu, degree_l, G_to_use, num_ys, ptr)


def takeuchi_liquid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Takeuchi starting conditions for a liquid dynamic compressible layer.

    TS72 Eqs. 95-102 (mu=0). Two independent solutions.
    """
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    cdef cpp_complex[double] K = cpp_complex[double](bulk_modulus.real, bulk_modulus.imag)
    c_takeuchi_liquid_dynamic_compressible(frequency, radius, density, K, degree_l, G_to_use, num_ys, ptr)
