# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp.complex cimport complex as cpp_complex


def kamata_solid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Kamata starting conditions for a solid dynamic compressible layer.

    KMN15 Eqs. B1-B16. Three independent solutions.

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
    c_kamata_solid_dynamic_compressible(frequency, radius, density, K, mu, degree_l, G_to_use, num_ys, ptr)


def kamata_solid_static_compressible(
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Kamata starting conditions for a solid static compressible layer.

    KMN15 Eqs. B1-B16 (w=0). Three independent solutions.
    """
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    cdef cpp_complex[double] K = cpp_complex[double](bulk_modulus.real, bulk_modulus.imag)
    cdef cpp_complex[double] mu = cpp_complex[double](shear_modulus.real, shear_modulus.imag)
    c_kamata_solid_static_compressible(radius, density, K, mu, degree_l, G_to_use, num_ys, ptr)


def kamata_solid_dynamic_incompressible(
        double frequency,
        double radius,
        double density,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Kamata starting conditions for a solid dynamic incompressible layer.

    KMN15 Eqs. B17-B28. Three independent solutions.
    """
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    cdef cpp_complex[double] mu = cpp_complex[double](shear_modulus.real, shear_modulus.imag)
    c_kamata_solid_dynamic_incompressible(frequency, radius, density, mu, degree_l, G_to_use, num_ys, ptr)


def kamata_liquid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Kamata starting conditions for a liquid dynamic compressible layer.

    KMN15 Eqs. B29-B37. Two independent solutions.
    """
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    cdef cpp_complex[double] K = cpp_complex[double](bulk_modulus.real, bulk_modulus.imag)
    c_kamata_liquid_dynamic_compressible(frequency, radius, density, K, degree_l, G_to_use, num_ys, ptr)


def kamata_liquid_dynamic_incompressible(
        double frequency,
        double radius,
        double density,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Kamata starting conditions for a liquid dynamic incompressible layer.

    KMN15 Eqs. B29-B37 (incompressible limit). Two independent solutions.
    """
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    c_kamata_liquid_dynamic_incompressible(frequency, radius, density, degree_l, G_to_use, num_ys, ptr)
