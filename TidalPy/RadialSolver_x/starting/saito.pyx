# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp.complex cimport complex as cpp_complex


def saito_liquid_static_incompressible(
        double radius,
        int degree_l,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate Saito starting conditions for a liquid static incompressible layer.

    S74 Eq. 19. One independent solution.

    Parameters
    ----------
    radius : float
        Radius [m].
    degree_l : int
        Tidal harmonic order.
    starting_conditions_view : complex[:, ::1]
        Output array of shape [1, num_ys].
    """
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    c_saito_liquid_static_incompressible(radius, degree_l, num_ys, ptr)
