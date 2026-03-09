# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libc.math cimport M_PI


def find_num_shooting_solutions(int layer_type, int is_static, int is_incompressible):
    """
    Return the number of independent shooting solutions for a layer.

    Parameters
    ----------
    layer_type : int
        0 = solid, 1 = liquid.
    is_static : int
        1 = static, 0 = dynamic.
    is_incompressible : int
        1 = incompressible, 0 = compressible.

    Returns
    -------
    num_solutions : int
    """
    return c_find_num_shooting_solutions(layer_type, is_static, is_incompressible)


def find_layer_diffeq_name(int layer_type, int is_static, int is_incompressible):
    """
    Return a string name identifying the ODE function for the given layer configuration.

    Parameters
    ----------
    layer_type : int
        0 = solid, 1 = liquid.
    is_static : int
        1 = static, 0 = dynamic.
    is_incompressible : int
        1 = incompressible, 0 = compressible.

    Returns
    -------
    name : str
    """
    cdef str layer_str   = 'solid' if layer_type == 0 else 'liquid'
    cdef str static_str  = 'static' if is_static == 1 else 'dynamic'
    cdef str incomp_str  = 'incompressible' if is_incompressible == 1 else 'compressible'
    return f'{layer_str}_{static_str}_{incomp_str}'
