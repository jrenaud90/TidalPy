# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.string cimport string as cpp_string, npos as cpp_npos


def find_starting_conditions(
        int layer_type,
        cpp_bool is_static,
        cpp_bool is_incompressible,
        cpp_bool use_kamata,
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        double complex[:, ::1] starting_conditions_view,
        cpp_bool run_y_checks = True):
    """
    Dispatch to the correct starting condition function.

    Parameters
    ----------
    layer_type : int
        0 = solid, 1 = liquid.
    is_static : bool
        True for static tides.
    is_incompressible : bool
        True for incompressible.
    use_kamata : bool
        True = Kamata (2015), False = Takeuchi & Saito (1972).
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
    run_y_checks : bool, optional
        If True, validate num_ys. Default True.
    """
    # Feedback
    cdef cpp_string message = cpp_string(b"No message set.")
    cdef cpp_bool success = False

    # starting conditions are passed as an array with shape [num_solutions, num_ys]
    cdef size_t num_ys = starting_conditions_view.shape[1]
    cdef cpp_complex[double]* ptr = <cpp_complex[double]*>&starting_conditions_view[0, 0]

    cdef cpp_complex[double] K = cpp_complex[double](bulk_modulus.real, bulk_modulus.imag)
    cdef cpp_complex[double] mu = cpp_complex[double](shear_modulus.real, shear_modulus.imag)

    c_find_starting_conditions(
        &success,
        message,
        layer_type,
        is_static,
        is_incompressible,
        use_kamata,
        frequency,
        radius,
        density,
        K,
        mu,
        degree_l,
        G_to_use,
        num_ys,
        ptr,
        run_y_checks
        )

    if not success:
        if message.find(cpp_string(b'not implemented')) != cpp_npos:
            raise NotImplementedError(message.decode('utf-8'))
        else:
            raise Exception(message.decode('utf-8'))
