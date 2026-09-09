# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
interp.pyx
Python/Cython wrapper for TidalPy's 1-D linear interpolation utility.

Exposes ``interp(x, xp, fp)``, a ``numpy.interp``-style linear interpolation
backed by the header-only C++ ``c_interp`` (Utilities_x/arrays/interp_.hpp).
"""

import numpy as np

from TidalPy.Utilities_x.arrays.interp cimport c_interp


def interp(x, xp, fp):
    """Linear interpolation of ``fp`` (sampled on ``xp``) at ``x`` (numpy.interp-style).

    Parameters
    ----------
    x : float or array_like
        Query coordinate(s) at which to interpolate.
    xp : array_like of float
        Sample coordinates, sorted in ascending order (length >= 1).
    fp : array_like of float
        Sample values, the same length as ``xp``.

    Returns
    -------
    float or numpy.ndarray
        The interpolated value(s). A Python ``float`` for scalar ``x``; otherwise
        a ``float64`` array with the shape of ``x``. Queries outside
        ``[xp[0], xp[-1]]`` clamp to the corresponding endpoint value.

    Assumptions
    -----------
    - ``xp`` is sorted ascending; results are undefined otherwise.
    - ``xp`` and ``fp`` have the same length.
    """
    cdef double[::1] xp_v = np.ascontiguousarray(xp, dtype=np.float64)
    cdef double[::1] fp_v = np.ascontiguousarray(fp, dtype=np.float64)
    cdef size_t n = xp_v.shape[0]
    if n == 0:
        raise ValueError("xp must have at least one element.")
    if <size_t>fp_v.shape[0] != n:
        raise ValueError("xp and fp must have the same length.")

    if np.ndim(x) == 0:
        return c_interp(<double>x, &xp_v[0], &fp_v[0], n, 0)

    x_in = np.ascontiguousarray(x, dtype=np.float64)
    cdef double[::1] x_v = x_in.ravel()
    cdef size_t m = x_v.shape[0]
    out = np.empty(m, dtype=np.float64)
    cdef double[::1] out_v = out
    cdef size_t i
    for i in range(m):
        out_v[i] = c_interp(x_v[i], &xp_v[0], &fp_v[0], n, 0)
    return out.reshape(np.shape(x))
