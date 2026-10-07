# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex

from TidalPy.RadialSolver.buffer_checks cimport cy_check_solution_buffer

from TidalPy.constants cimport cy_resolve_G, get_shared_config_address, set_tidalpy_config_ptr

# Wire this DLL's shared pointer to the process-wide TidalPy config singleton, whose G cy_resolve_G reads.
set_tidalpy_config_ptr(get_shared_config_address())

NOT_CONVERGED_MESSAGE = (
    "The power series starting conditions refused to start: the series did not converge at this radius, or the "
    "solutions grow too steeply there (a weak solid, or a dynamic liquid at long periods); use the Takeuchi or Kamata "
    "starting conditions, or a smaller radius.")


cdef void cy_power_series_solid(
        str name,
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        cpp_bool is_incompressible,
        int degree_l,
        object G_to_use,
        double complex[:, ::1] starting_conditions_view) except *:
    """Shared body of the solid wrappers: check the buffer, sum the series, raise when it refused to start."""
    cy_check_solution_buffer(
        name + "'s starting_conditions_view", starting_conditions_view.shape[0], starting_conditions_view.shape[1], 0,
        False)
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    cdef cpp_complex[double] K = cpp_complex[double](bulk_modulus.real, bulk_modulus.imag)
    cdef cpp_complex[double] mu = cpp_complex[double](shear_modulus.real, shear_modulus.imag)
    cdef cpp_bool converged = c_power_series_solid(
        frequency,
        radius,
        density,
        K,
        mu,
        is_incompressible,
        degree_l,
        cy_resolve_G(G_to_use),
        num_ys,
        ptr)
    if not converged:
        raise RuntimeError(NOT_CONVERGED_MESSAGE)


cdef void cy_power_series_liquid(
        str name,
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        cpp_bool is_incompressible,
        int degree_l,
        object G_to_use,
        double complex[:, ::1] starting_conditions_view) except *:
    """Shared body of the dynamic-liquid wrappers."""
    cy_check_solution_buffer(
        name + "'s starting_conditions_view", starting_conditions_view.shape[0], starting_conditions_view.shape[1], 1,
        False)
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]
    cdef cpp_complex[double] K = cpp_complex[double](bulk_modulus.real, bulk_modulus.imag)
    cdef cpp_bool converged = c_power_series_liquid_dynamic(
        frequency,
        radius,
        density,
        K,
        is_incompressible,
        degree_l,
        cy_resolve_G(G_to_use),
        num_ys,
        ptr)
    if not converged:
        raise RuntimeError(NOT_CONVERGED_MESSAGE)


def power_series_solid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        object G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate power series starting conditions for a solid dynamic compressible layer.

    Martens (2016) Sec. 4.2.8 (after Smylie 2013), in TidalPy's y convention and summed to convergence. Three
    independent solutions.

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
    G_to_use : float or None
        Gravitational constant [m3 kg-1 s-2]; None takes the TidalPy configuration's value (SciPy's G).
    starting_conditions_view : complex[:, ::1]
        Output array of shape [num_solutions, num_ys].

    Raises
    ------
    RuntimeError
        When the series refuses to start at this radius (see power_series_.hpp).
    """
    cy_power_series_solid(
        "power_series_solid_dynamic_compressible",
        frequency,
        radius,
        density,
        bulk_modulus,
        shear_modulus,
        False,
        degree_l,
        G_to_use,
        starting_conditions_view)


def power_series_solid_static_compressible(
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        object G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate power series starting conditions for a solid static compressible layer (w = 0).

    Three independent solutions; arguments as ``power_series_solid_dynamic_compressible`` without the frequency.
    """
    cy_power_series_solid(
        "power_series_solid_static_compressible",
        0.0,
        radius,
        density,
        bulk_modulus,
        shear_modulus,
        False,
        degree_l,
        G_to_use,
        starting_conditions_view)


def power_series_solid_dynamic_incompressible(
        double frequency,
        double radius,
        double density,
        double complex shear_modulus,
        int degree_l,
        object G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate power series starting conditions for a solid dynamic incompressible layer.

    Three independent solutions; arguments as ``power_series_solid_dynamic_compressible`` without the bulk modulus.
    """
    cy_power_series_solid(
        "power_series_solid_dynamic_incompressible",
        frequency,
        radius,
        density,
        0.0,
        shear_modulus,
        True,
        degree_l,
        G_to_use,
        starting_conditions_view)


def power_series_solid_static_incompressible(
        double radius,
        double density,
        double complex shear_modulus,
        int degree_l,
        object G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate power series starting conditions for a solid static incompressible layer (w = 0).

    Three independent solutions; arguments as ``power_series_solid_dynamic_incompressible`` without the frequency.
    """
    cy_power_series_solid(
        "power_series_solid_static_incompressible",
        0.0,
        radius,
        density,
        0.0,
        shear_modulus,
        True,
        degree_l,
        G_to_use,
        starting_conditions_view)


def power_series_liquid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        int degree_l,
        object G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate power series starting conditions for a liquid dynamic compressible layer.

    Two independent solutions of (y1, y2, y5, y6); arguments as ``power_series_solid_dynamic_compressible`` without
    the shear modulus.
    """
    cy_power_series_liquid(
        "power_series_liquid_dynamic_compressible",
        frequency,
        radius,
        density,
        bulk_modulus,
        False,
        degree_l,
        G_to_use,
        starting_conditions_view)


def power_series_liquid_dynamic_incompressible(
        double frequency,
        double radius,
        double density,
        int degree_l,
        object G_to_use,
        double complex[:, ::1] starting_conditions_view):
    """
    Calculate power series starting conditions for a liquid dynamic incompressible layer.

    Two independent solutions of (y1, y2, y5, y6); the series ends after its first term.
    """
    cy_power_series_liquid(
        "power_series_liquid_dynamic_incompressible",
        frequency,
        radius,
        density,
        0.0,
        True,
        degree_l,
        G_to_use,
        starting_conditions_view)
