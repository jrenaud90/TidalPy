# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex

from TidalPy.RadialSolver.buffer_checks cimport (
    c_layer_num_solutions, cy_check_solution_rows, cy_resolve_num_ys)


def solve_upper_y_at_interface(
        double complex[:, ::1] lower_layer_y_view,
        double complex[:, ::1] upper_layer_y_view,
        int lower_layer_type,
        cpp_bool lower_is_static,
        int upper_layer_type,
        cpp_bool upper_is_static,
        double interface_gravity,
        double liquid_density,
        double G_to_use,
        object max_num_y = None):
    """
    Calculate the initial conditions for an overlying layer given the lower layer's y values.

    Parameters
    ----------
    lower_layer_y_view : complex[:, ::1]
        Lower layer y values, shape [num_sols_lower, num_y]: a row for each of the lower layer's solutions (3 for a
        solid, 2 for a dynamic liquid, 1 for a static liquid).
    upper_layer_y_view : complex[:, ::1]
        Output upper layer y values, shape [num_sols_upper, num_y], a row for each of the upper layer's solutions
        and the same column count.
    lower_layer_type : int
        0 = solid, 1 = liquid.
    lower_is_static : bool
    upper_layer_type : int
        0 = solid, 1 = liquid.
    upper_is_static : bool
    interface_gravity : float
        Gravity at the interface [m s-2].
    liquid_density : float
        Density of the liquid at the interface [kg m-3].
    G_to_use : float
        Gravitational constant.
    max_num_y : int, optional
        The y values per solution, which must equal the arrays' column count; None (default) takes it from them.

    Raises
    ------
    ValueError
        If either array has too few rows for its layer's solutions, too few columns for its ys, or the two column
        counts differ.
    """
    cy_check_solution_rows(
        "solve_upper_y_at_interface's lower_layer_y_view", lower_layer_y_view.shape[0], lower_layer_type,
        lower_is_static)
    cy_check_solution_rows(
        "solve_upper_y_at_interface's upper_layer_y_view", upper_layer_y_view.shape[0], upper_layer_type,
        upper_is_static)
    if lower_layer_y_view.shape[1] != upper_layer_y_view.shape[1]:
        raise ValueError(
            f"TidalPy: solve_upper_y_at_interface needs the two arrays to have the same column count; got "
            f"{lower_layer_y_view.shape[1]} and {upper_layer_y_view.shape[1]}.")
    cy_resolve_num_ys(
        "solve_upper_y_at_interface", lower_layer_y_view.shape[1], max_num_y, lower_layer_type, lower_is_static)
    cdef size_t num_ys = cy_resolve_num_ys(
        "solve_upper_y_at_interface", upper_layer_y_view.shape[1], max_num_y, upper_layer_type, upper_is_static)
    cdef size_t num_sols_lower = c_layer_num_solutions(lower_layer_type, lower_is_static)
    cdef size_t num_sols_upper = c_layer_num_solutions(upper_layer_type, upper_is_static)

    c_solve_upper_y_at_interface(
        <cpp_complex[double]*>&lower_layer_y_view[0, 0],
        <cpp_complex[double]*>&upper_layer_y_view[0, 0],
        num_sols_lower,
        num_sols_upper,
        num_ys,
        lower_layer_type,
        lower_is_static,
        upper_layer_type,
        upper_is_static,
        interface_gravity,
        liquid_density,
        G_to_use
        )
