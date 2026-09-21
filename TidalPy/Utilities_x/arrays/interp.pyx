# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python wrapper for TidalPy's 1-D linear interpolation, backed by ``c_interp`` in ``interp_.hpp``."""

import numpy as np

from libcpp.vector cimport vector

from TidalPy.Utilities_x.arrays.interp cimport c_interp, c_partition_radius_by_layer


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
        return c_interp(
                <double>x,
                &xp_v[0],
                &fp_v[0],
                n,
                0)

    x_in = np.ascontiguousarray(x, dtype=np.float64)
    cdef double[::1] x_v = x_in.ravel()
    cdef size_t m = x_v.shape[0]
    out = np.empty(m, dtype=np.float64)
    cdef double[::1] out_v = out
    cdef size_t i
    for i in range(m):
        out_v[i] = c_interp(x_v[i], &xp_v[0], &fp_v[0], n, 0)
    return out.reshape(np.shape(x))


def partition_radius_by_layer(double[::1] radius not None, double[::1] upper_radius_bylayer not None):
    """Split an ascending radius array into one run of slices per layer.

    A layered profile repeats each interface radius, once for the layer below and once for the layer above.
    This applies the rule that decides which copy belongs to which layer, and it is the same C++ routine the
    equation-of-state solution and the world radial solver partition with, so the three cannot disagree.

    Parameters
    ----------
    radius : np.ndarray[dtype=np.float64]
        Slice radii [m], ascending, with interface radii appearing twice.
    upper_radius_bylayer : np.ndarray[dtype=np.float64]
        Upper radius of each layer [m], inner to outer.

    Returns
    -------
    first_slice : np.ndarray[dtype=np.uint64]
        Index of each layer's first slice.
    num_slices : np.ndarray[dtype=np.uint64]
        Number of slices in each layer; zero when a layer caught none.
    """
    cdef size_t num_slices_in = radius.shape[0]
    cdef size_t num_layers    = upper_radius_bylayer.shape[0]
    cdef vector[size_t] first_out
    cdef vector[size_t] count_out
    if num_slices_in == 0 or num_layers == 0:
        return np.zeros(num_layers, dtype=np.uint64), np.zeros(num_layers, dtype=np.uint64)
    with nogil:
        c_partition_radius_by_layer(
            &radius[0], num_slices_in, &upper_radius_bylayer[0], num_layers, first_out, count_out)
    cdef size_t layer_i
    first_arr = np.empty(num_layers, dtype=np.uint64)
    count_arr = np.empty(num_layers, dtype=np.uint64)
    for layer_i in range(num_layers):
        first_arr[layer_i] = first_out[layer_i]
        count_arr[layer_i] = count_out[layer_i]
    return first_arr, count_arr
