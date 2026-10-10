# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex

from TidalPy.RadialSolver.buffer_checks cimport (
    c_layer_num_solutions, cy_check_at_least, cy_check_solution_rows, cy_resolve_num_ys)


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
        object max_num_y = None):
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
    max_num_y : int, optional
        The y values per solution, which must equal the array's column count; None (default) takes it from it.

    Raises
    ------
    ValueError
        If constant_vector_view is shorter than 3, layer_above_constant_vector_view shorter than the layer above's
        solutions, or uppermost_y_per_solution_view too small for this layer's solutions and ys.
    """
    cy_check_at_least("top_to_bottom_interface_bc's constant_vector_view", constant_vector_view.shape[0], 3, "entries")
    cy_check_at_least(
        "top_to_bottom_interface_bc's layer_above_constant_vector_view", layer_above_constant_vector_view.shape[0],
        <Py_ssize_t>c_layer_num_solutions(layer_above_type, layer_above_is_static), "entries")
    cy_check_solution_rows(
        "top_to_bottom_interface_bc's uppermost_y_per_solution_view", uppermost_y_per_solution_view.shape[0],
        layer_type, layer_is_static)
    cdef size_t num_ys = cy_resolve_num_ys(
        "top_to_bottom_interface_bc", uppermost_y_per_solution_view.shape[1], max_num_y, layer_type,
        layer_is_static)
    cdef size_t num_sols = c_layer_num_solutions(layer_type, layer_is_static)

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
        num_sols,
        num_ys
        )
