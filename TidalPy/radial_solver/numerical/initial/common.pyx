# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from scipy.special.cython_special cimport spherical_jn

from TidalPy.utilities.math.complex cimport csqrt, cipow
from TidalPy.utilities.math.special_x cimport double_factorial


cdef double complex z_calc(
        double complex x_squared,
        unsigned char degree_l) noexcept nogil:
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


cdef void takeuchi_phi_psi(
        double complex z,
        unsigned char degree_l,
        double complex* phi_ptr,
        double complex* phi_lplus1_ptr,
        double complex* psi_ptr,
        ) noexcept nogil:
    """ Calculate the two (plus one) functions used to find initial conditions for shooting method.

    References
    ----------
    TS72 Eq. 103

    """

    cdef double l_dbl_factorial   = double_factorial(2 * degree_l + 1)
    cdef double lp1_dbl_factorial = double_factorial(2 * degree_l + 2)
    
    cdef double complex zl   = cipow(z, degree_l)
    cdef double complex zlp1 = cipow(z, degree_l + 1)

    phi_ptr[0]        = l_dbl_factorial * (spherical_jn(degree_l, z) / zl)
    phi_lplus1_ptr[0] = lp1_dbl_factorial * (spherical_jn(degree_l + 1, z) / zlp1)
    psi_ptr[0]        = (2. * (2. * degree_l + 3.) / (z * z)) * (1. - phi_ptr[0])
