# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex


def top_to_bottom_interface_bc(
        double complex[::1] constant_vector_view,
        double complex[::1] layer_above_constant_vector_view,
        double complex[:, ::1] uppermost_y_per_solution_view,
        double gravity_upper,
        double layer_above_lower_gravity,
        double density_upper,
        double layer_above_lower_density,
        int layer_type,
        int layer_above_type,
        cpp_bool layer_is_static,
        cpp_bool layer_above_is_static,
        cpp_bool layer_is_incomp,
        cpp_bool layer_above_is_incomp,
        int max_num_y = 6):
    """
    Calculate the constant vector for a layer given the layer above's constants (top-to-bottom).

    Used during the collapse phase of the shooting method.

    Parameters
    ----------
    constant_vector_view : complex[::1]
        Output constant vector for this layer.
    layer_above_constant_vector_view : complex[::1]
        Constant vector from the layer above.
    uppermost_y_per_solution_view : complex[:, ::1]
        Y values at the top of this layer, shape [num_sols, max_num_y].
    gravity_upper : float
        Gravity at the top of this layer [m s-2].
    layer_above_lower_gravity : float
        Gravity at the bottom of the layer above [m s-2].
    density_upper : float
        Density at the top of this layer [kg m-3].
    layer_above_lower_density : float
        Density at the bottom of the layer above [kg m-3].
    layer_type : int
        0 = solid, 1 = liquid.
    layer_above_type : int
        0 = solid, 1 = liquid.
    layer_is_static : bool
    layer_above_is_static : bool
    layer_is_incomp : bool
    layer_above_is_incomp : bool
    max_num_y : int, optional
        Maximum number of y values per solution. Default 6.
    """
    cdef size_t num_sols = uppermost_y_per_solution_view.shape[0]

    c_top_to_bottom_interface_bc(
        <cpp_complex[double]*>&constant_vector_view[0],
        <cpp_complex[double]*>&layer_above_constant_vector_view[0],
        <cpp_complex[double]*>&uppermost_y_per_solution_view[0, 0],
        gravity_upper,
        layer_above_lower_gravity,
        density_upper,
        layer_above_lower_density,
        layer_type,
        layer_above_type,
        layer_is_static,
        layer_above_is_static,
        layer_is_incomp,
        layer_above_is_incomp,
        num_sols,
        max_num_y
        )
