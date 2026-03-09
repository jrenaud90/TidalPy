# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex


def solve_upper_y_at_interface(
        double complex[:, ::1] lower_layer_y_view,
        double complex[:, ::1] upper_layer_y_view,
        int lower_layer_type,
        cpp_bool lower_is_static,
        cpp_bool lower_is_incompressible,
        int upper_layer_type,
        cpp_bool upper_is_static,
        cpp_bool upper_is_incompressible,
        double interface_gravity,
        double liquid_density,
        double G_to_use,
        int max_num_y = 6):
    """
    Calculate the initial conditions for an overlying layer given the lower layer's y values.

    Parameters
    ----------
    lower_layer_y_view : complex[:, ::1]
        Lower layer y values, shape [num_sols_lower, max_num_y].
    upper_layer_y_view : complex[:, ::1]
        Output upper layer y values, shape [num_sols_upper, max_num_y].
    lower_layer_type : int
        0 = solid, 1 = liquid.
    lower_is_static : bool
    lower_is_incompressible : bool
    upper_layer_type : int
        0 = solid, 1 = liquid.
    upper_is_static : bool
    upper_is_incompressible : bool
    interface_gravity : float
        Gravity at the interface [m s-2].
    liquid_density : float
        Density of the liquid at the interface [kg m-3].
    G_to_use : float
        Gravitational constant.
    max_num_y : int, optional
        Maximum number of y values per solution. Default 6.
    """
    cdef size_t num_sols_lower = lower_layer_y_view.shape[0]
    cdef size_t num_sols_upper = upper_layer_y_view.shape[0]

    c_solve_upper_y_at_interface(
        <cpp_complex[double]*>&lower_layer_y_view[0, 0],
        <cpp_complex[double]*>&upper_layer_y_view[0, 0],
        num_sols_lower,
        num_sols_upper,
        max_num_y,
        lower_layer_type,
        lower_is_static,
        lower_is_incompressible,
        upper_layer_type,
        upper_is_static,
        upper_is_incompressible,
        interface_gravity,
        liquid_density,
        G_to_use
        )
