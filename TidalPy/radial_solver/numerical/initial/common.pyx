# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from scipy.special.cython_special cimport spherical_jn

from TidalPy.utilities.math.complex cimport csqrt


cdef double complex z_calc(double complex x_squared, unsigned int degree_l) noexcept nogil:
    """ Calculates the z function using spherical Bessel function, see Eq. B14 of KMN15.

    References
    ----------
    TS72 Eqs. 96, 97
    KMN15 Eq. B14

    Parameters
    ----------
    x_squared : double complex
        Expression passed to the Bessel function.
    degree_l : unsigned int
        Tidal harmonic order.

    Returns
    -------
    z : double complex
        Result
    """

    cdef double complex x, numerator, denominator
    
    x = csqrt(x_squared)
    numerator   = x * spherical_jn(degree_l + 1, x)
    denominator = spherical_jn(degree_l, x)

    return numerator / denominator
