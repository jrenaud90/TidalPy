# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport tgamma

import numpy as np
cimport numpy as cnp
cnp.import_array()

from TidalPy.constants cimport d_NAN


cdef double[51] pre_calculated_doubles
cdef double* pre_calculated_doubles_ptr = &pre_calculated_doubles[0]
pre_calculated_doubles_ptr[  0] = 1.00
pre_calculated_doubles_ptr[  1] = 1.00
pre_calculated_doubles_ptr[  2] = 2.00
pre_calculated_doubles_ptr[  3] = 3.00
pre_calculated_doubles_ptr[  4] = 8.00
pre_calculated_doubles_ptr[  5] = 15.00
pre_calculated_doubles_ptr[  6] = 48.00
pre_calculated_doubles_ptr[  7] = 105.00
pre_calculated_doubles_ptr[  8] = 384.00
pre_calculated_doubles_ptr[  9] = 945.00
pre_calculated_doubles_ptr[ 10] = 3840.00
pre_calculated_doubles_ptr[ 11] = 10395.00
pre_calculated_doubles_ptr[ 12] = 46080.00
pre_calculated_doubles_ptr[ 13] = 135135.00
pre_calculated_doubles_ptr[ 14] = 645120.00
pre_calculated_doubles_ptr[ 15] = 2027025.00
pre_calculated_doubles_ptr[ 16] = 10321920.00
pre_calculated_doubles_ptr[ 17] = 34459425.00
pre_calculated_doubles_ptr[ 18] = 185794560.00
pre_calculated_doubles_ptr[ 19] = 654729075.00
pre_calculated_doubles_ptr[ 20] = 3715891200.00
pre_calculated_doubles_ptr[ 21] = 13749310575.00
pre_calculated_doubles_ptr[ 22] = 81749606400.00
pre_calculated_doubles_ptr[ 23] = 316234143225.00
pre_calculated_doubles_ptr[ 24] = 1961990553600.00
pre_calculated_doubles_ptr[ 25] = 7905853580625.01
pre_calculated_doubles_ptr[ 26] = 51011754393599.96
pre_calculated_doubles_ptr[ 27] = 213458046676875.16
pre_calculated_doubles_ptr[ 28] = 1428329123020799.25
pre_calculated_doubles_ptr[ 29] = 6190283353629379.00
pre_calculated_doubles_ptr[ 30] = 42849873690623960.00
pre_calculated_doubles_ptr[ 31] = 191898783962510816.00
pre_calculated_doubles_ptr[ 32] = 1371195958099966720.00
pre_calculated_doubles_ptr[ 33] = 6332659870762856448.00
pre_calculated_doubles_ptr[ 34] = 46620662575398879232.00
pre_calculated_doubles_ptr[ 35] = 221643095476699824128.00
pre_calculated_doubles_ptr[ 36] = 1678343852714360832000.00
pre_calculated_doubles_ptr[ 37] = 8200794532637892935680.00
pre_calculated_doubles_ptr[ 38] = 63777066403145720004608.00
pre_calculated_doubles_ptr[ 39] = 319830986772877752139776.00
pre_calculated_doubles_ptr[ 40] = 2551082656125828464640000.00
pre_calculated_doubles_ptr[ 41] = 13113070457687983475654656.00
pre_calculated_doubles_ptr[ 42] = 107145471557284812694749184.00
pre_calculated_doubles_ptr[ 43] = 563862029680583787669356544.00
pre_calculated_doubles_ptr[ 44] = 4714400748520528253875650560.00
pre_calculated_doubles_ptr[ 45] = 25373791335626273400058544128.00
pre_calculated_doubles_ptr[ 46] = 216862434431944368947512475648.00
pre_calculated_doubles_ptr[ 47] = 1192568192774434172503588864000.00
pre_calculated_doubles_ptr[ 48] = 10409396852733329709480598831104.00
pre_calculated_doubles_ptr[ 49] = 58435841445947288807899666579456.00
pre_calculated_doubles_ptr[ 50] = 520469842636666553028024352112640.00


cdef double cf_double_factorial(int n) noexcept nogil:

    if n < 51:
        # Use precalculated doubles
        return pre_calculated_doubles_ptr[n]
    elif n >= 51 and n < 171:
        # Calculate using recursion
        return tgamma(<double>n + 1.) / cf_double_factorial(n - 1)
    else:
        # ValueError('C function `tgamma` experiences overflow for l > 170.')
        return d_NAN


def double_factorial(int n):
    
    if n >= 171:
        raise ValueError('C function `tgamma` experiences overflow for l > 170.')

    return cf_double_factorial(n)

# ======================================================================================================================
# XSF - Spherical Bessel Functions
# ======================================================================================================================
# Python-visible scalar wrappers
def _spherical_jn_scalar(long n, double x):
    """Compute j_n(x) for scalar real x."""
    return sph_bessel_j(n, x)

def _spherical_jn_scalar_complex(long n, double complex z):
    """Compute j_n(z) for scalar complex z."""
    return sph_bessel_j(n, z)

def _spherical_jn_d_scalar(long n, double x):
    """Compute j_n'(x) for scalar real x."""
    return sph_bessel_j_jac(n, x)


def _spherical_yn_scalar(long n, double x):
    """Compute y_n(x) for scalar real x."""
    return sph_bessel_y(n, x)

def _spherical_yn_scalar_complex(long n, double complex z):
    """Compute y_n(z) for scalar complex z."""
    return sph_bessel_y(n, z)

def _spherical_yn_d_scalar(long n, double x):
    """Compute y_n'(x) for scalar real x."""
    return sph_bessel_y_jac(n, x)


def _spherical_in_scalar(long n, double x):
    """Compute i_n(x) for scalar real x."""
    return sph_bessel_i(n, x)

def _spherical_in_scalar_complex(long n, double complex z):
    """Compute i_n(z) for scalar complex z."""
    return sph_bessel_i(n, z)

def _spherical_in_d_scalar(long n, double x):
    """Compute i_n'(x) for scalar real x."""
    return sph_bessel_i_jac(n, x)


def _spherical_kn_scalar(long n, double x):
    """Compute k_n(x) for scalar real x."""
    return sph_bessel_k(n, x)

def _spherical_kn_scalar_complex(long n, double complex z):
    """Compute k_n(z) for scalar complex z."""
    return sph_bessel_k(n, z)

def _spherical_kn_d_scalar(long n, double x):
    """Compute k_n'(x) for scalar real x."""
    return sph_bessel_k_jac(n, x)


# Vectorized (numpy array) wrappers
def spherical_jn_array(long n, double[::1] x not None, bint derivative=False):
    """
    Compute spherical Bessel function j_n(x) over an array.

    Parameters
    ----------
    n : int
        Order (non-negative).
    x : 1D array of float64
        Arguments.
    derivative : bool
        If True, compute the derivative j_n'(x).

    Returns
    -------
    out : ndarray, shape (len(x),)
    """
    cdef Py_ssize_t N = x.shape[0]
    cdef cnp.ndarray[double, ndim=1] out = np.empty(N, dtype=np.float64)
    cdef double[::1] out_view = out
    cdef Py_ssize_t i

    if derivative:
        for i in range(N):
            out_view[i] = sph_bessel_j_jac(n, x[i])
    else:
        for i in range(N):
            out_view[i] = sph_bessel_j(n, x[i])

    return out


def spherical_yn_array(long n, double[::1] x not None, bint derivative=False):
    """
    Compute spherical Bessel function y_n(x) over an array.

    Parameters
    ----------
    n : int
        Order (non-negative).
    x : 1D array of float64
        Arguments.
    derivative : bool
        If True, compute the derivative y_n'(x).

    Returns
    -------
    out : ndarray, shape (len(x),)
    """
    cdef Py_ssize_t N = x.shape[0]
    cdef cnp.ndarray[double, ndim=1] out = np.empty(N, dtype=np.float64)
    cdef double[::1] out_view = out
    cdef Py_ssize_t i

    if derivative:
        for i in range(N):
            out_view[i] = sph_bessel_y_jac(n, x[i])
    else:
        for i in range(N):
            out_view[i] = sph_bessel_y(n, x[i])

    return out


def spherical_in_array(long n, double[::1] x not None, bint derivative=False):
    """
    Compute modified spherical Bessel function i_n(x) over an array.

    Parameters
    ----------
    n : int
        Order (non-negative).
    x : 1D array of float64
        Arguments.
    derivative : bool
        If True, compute the derivative i_n'(x).

    Returns
    -------
    out : ndarray, shape (len(x),)
    """
    cdef Py_ssize_t N = x.shape[0]
    cdef cnp.ndarray[double, ndim=1] out = np.empty(N, dtype=np.float64)
    cdef double[::1] out_view = out
    cdef Py_ssize_t i

    if derivative:
        for i in range(N):
            out_view[i] = sph_bessel_i_jac(n, x[i])
    else:
        for i in range(N):
            out_view[i] = sph_bessel_i(n, x[i])

    return out


def spherical_kn_array(long n, double[::1] x not None, bint derivative=False):
    """
    Compute modified spherical Bessel function k_n(x) over an array.

    Parameters
    ----------
    n : int
        Order (non-negative).
    x : 1D array of float64
        Arguments.
    derivative : bool
        If True, compute the derivative k_n'(x).

    Returns
    -------
    out : ndarray, shape (len(x),)
    """
    cdef Py_ssize_t N = x.shape[0]
    cdef cnp.ndarray[double, ndim=1] out = np.empty(N, dtype=np.float64)
    cdef double[::1] out_view = out
    cdef Py_ssize_t i

    if derivative:
        for i in range(N):
            out_view[i] = sph_bessel_k_jac(n, x[i])
    else:
        for i in range(N):
            out_view[i] = sph_bessel_k(n, x[i])

    return out


def spherical_jn_complex_array(long n, double complex[::1] z not None):
    """
    Compute spherical Bessel function j_n(z) over a complex array.

    Parameters
    ----------
    n : int
        Order (non-negative).
    z : 1D array of complex128
        Arguments.

    Returns
    -------
    out : ndarray of complex128, shape (len(z),)
    """
    cdef Py_ssize_t N = z.shape[0]
    cdef cnp.ndarray[double complex, ndim=1] out = np.empty(N, dtype=np.complex128)
    cdef double complex[::1] out_view = out
    cdef Py_ssize_t i

    for i in range(N):
        out_view[i] = sph_bessel_j(n, z[i])

    return out


def spherical_yn_complex_array(long n, double complex[::1] z not None):
    """
    Compute spherical Bessel function y_n(z) over a complex array.

    Parameters
    ----------
    n : int
        Order (non-negative).
    z : 1D array of complex128
        Arguments.

    Returns
    -------
    out : ndarray of complex128, shape (len(z),)
    """
    cdef Py_ssize_t N = z.shape[0]
    cdef cnp.ndarray[double complex, ndim=1] out = np.empty(N, dtype=np.complex128)
    cdef double complex[::1] out_view = out
    cdef Py_ssize_t i

    for i in range(N):
        out_view[i] = sph_bessel_y(n, z[i])

    return out


def spherical_in_complex_array(long n, double complex[::1] z not None):
    """
    Compute modified spherical Bessel function i_n(z) over a complex array.

    Parameters
    ----------
    n : int
        Order (non-negative).
    z : 1D array of complex128
        Arguments.

    Returns
    -------
    out : ndarray of complex128, shape (len(z),)
    """
    cdef Py_ssize_t N = z.shape[0]
    cdef cnp.ndarray[double complex, ndim=1] out = np.empty(N, dtype=np.complex128)
    cdef double complex[::1] out_view = out
    cdef Py_ssize_t i

    for i in range(N):
        out_view[i] = sph_bessel_i(n, z[i])

    return out


def spherical_kn_complex_array(long n, double complex[::1] z not None):
    """
    Compute modified spherical Bessel function k_n(z) over a complex array.

    Parameters
    ----------
    n : int
        Order (non-negative).
    z : 1D array of complex128
        Arguments.

    Returns
    -------
    out : ndarray of complex128, shape (len(z),)
    """
    cdef Py_ssize_t N = z.shape[0]
    cdef cnp.ndarray[double complex, ndim=1] out = np.empty(N, dtype=np.complex128)
    cdef double complex[::1] out_view = out
    cdef Py_ssize_t i

    for i in range(N):
        out_view[i] = sph_bessel_k(n, z[i])

    return out
