# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex

from TidalPy.RadialSolver.buffer_checks cimport (
    cy_check_at_least, cy_check_solution_rows, cy_resolve_num_ys)


def apply_surface_bc(
        double complex[::1] constant_vector_view,
        double[::1] bc_view,
        double complex[:, ::1] uppermost_y_per_solution_view,
        double surface_gravity,
        double G_to_use,
        size_t ytype_i,
        int layer_type,
        cpp_bool layer_is_static,
        cpp_bool layer_is_incomp,
        object max_num_y = None):
    """
    Apply surface boundary conditions by solving a linear system.

    Parameters
    ----------
    constant_vector_view : complex[::1]
        Output constant vector (overwritten with solution), at least 3 long; entries past the layer's solutions are
        NaN.
    bc_view : float[::1]
        Boundary condition values, 3 per y-type, at least 3 * (ytype_i + 1) long.
    uppermost_y_per_solution_view : complex[:, ::1]
        Y values at surface for each solution, a row for each of the layer's solutions.
    surface_gravity : float
        Gravitational acceleration at surface [m s-2].
    G_to_use : float
        Gravitational constant.
    ytype_i : int
        Y-type index (tidal=0, loading=1, etc.).
    layer_type : int
        0=solid, 1=liquid.
    layer_is_static : bool
    layer_is_incomp : bool
    max_num_y : int, optional
        The y values per solution, which must equal the array's column count; None (default) takes it from it.

    Returns
    -------
    info : int
        LAPACK solver info (0=success).

    Raises
    ------
    ValueError
        If an array is too small for the layer's solutions, ys, or the y-type's boundary conditions.
    """
    cy_check_at_least("apply_surface_bc's constant_vector_view", constant_vector_view.shape[0], 3, "entries")
    cy_check_at_least("apply_surface_bc's bc_view", bc_view.shape[0], 3 * (ytype_i + 1), "entries")
    cy_check_solution_rows(
        "apply_surface_bc's uppermost_y_per_solution_view", uppermost_y_per_solution_view.shape[0], layer_type,
        layer_is_static)
    cdef size_t num_ys = cy_resolve_num_ys(
        "apply_surface_bc", uppermost_y_per_solution_view.shape[1], max_num_y, layer_type, layer_is_static)
    cdef size_t num_sols = uppermost_y_per_solution_view.shape[0]
    cdef int bc_solution_info = 0

    c_apply_surface_bc(
        <cpp_complex[double]*>&constant_vector_view[0],
        &bc_solution_info,
        &bc_view[0],
        <cpp_complex[double]*>&uppermost_y_per_solution_view[0, 0],
        surface_gravity,
        G_to_use,
        num_sols,
        num_ys,
        ytype_i,
        layer_type,
        layer_is_static,
        layer_is_incomp
        )

    return bc_solution_info
