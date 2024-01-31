# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport abs

from scipy.special.cython_special cimport spherical_jn

from TidalPy.utilities.math.complex cimport cf_csqrt, cf_cipow
from TidalPy.utilities.math.special_x cimport cf_double_factorial

cdef double complex cf_z_calc(
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

    cdef double complex x, z, numerator, denominator
    
    x = cf_csqrt(x_squared)
    numerator   = x * spherical_jn(degree_l + 1, x)
    denominator = spherical_jn(degree_l, x)

    z = numerator / denominator

    # The convergence of this function depends on the absolute size of the real part of x^2.
    # The table provided outside this function (`Z_CALC_MAX_L`) was tested on 2021/12/02 to find the smallest
    #    max l that still allowed convergence. It is not comprehensive, thus the possibility of an error being
    #    thrown below.
    # cdef double[11][2] Z_CALC_MAX_L
    # Z_CALC_MAX_L[0][0] = 0.1
    # Z_CALC_MAX_L[0][1] = 7
    # Z_CALC_MAX_L[1][0] = 1.
    # Z_CALC_MAX_L[1][1] = 7
    # Z_CALC_MAX_L[2][0] = 10.
    # Z_CALC_MAX_L[2][1] = 8
    # Z_CALC_MAX_L[3][0] = 100.
    # Z_CALC_MAX_L[3][1] = 17.
    # Z_CALC_MAX_L[4][0] = 500.
    # Z_CALC_MAX_L[4][1] = 32.
    # Z_CALC_MAX_L[5][0] = 1000.
    # Z_CALC_MAX_L[5][1] = 44.
    # Z_CALC_MAX_L[6][0] = 5000.
    # Z_CALC_MAX_L[6][1] = 88.
    # Z_CALC_MAX_L[7][0] = 10000
    # Z_CALC_MAX_L[7][1] = 117
    # Z_CALC_MAX_L[8][0] = 50000
    # Z_CALC_MAX_L[8][1] = 245
    # Z_CALC_MAX_L[9][0] = 100000
    # Z_CALC_MAX_L[9][1] = 340
    # Z_CALC_MAX_L[10][0] = 500000
    # Z_CALC_MAX_L[10][1] = 737

    # cdef double x2_real = abs(x_squared.real)
    # cdef double min_val, max_l
    # cdef size_t max_l_to_use = 0
    # cdef size_t i, l_fake, l
    # for i in range(11):
    #     min_val = Z_CALC_MAX_L[i][0]
    #     max_l = Z_CALC_MAX_L[i][1]
    #     if x2_real <= min_val:
    #         max_l_to_use = <size_t> max_l + (degree_l - 2)
    #         break
    # if max_l_to_use == 0:
    #     max_l_to_use = 1000

    # cdef double complex z
    # z = x_squared / (2. * max_l_to_use + 3.)
    # for l_fake in range(degree_l, max_l_to_use + 1):
    #     l = max_l_to_use - l_fake + degree_l
    #     z = x_squared / ((2. * l + 1.) - z)

    return z


cdef void cf_takeuchi_phi_psi(
        double complex z2,
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

    # Floating point errors prevent us from using the exact definition of these functions (See Issue # 41)
    # Instead we use the limiting version of the functions.
    # However, we leave the full definition at the bottom of this function for reference.

    cdef double complex z4   = z2 * z2
    cdef double complex z6   = z4 * z2
    cdef double complex z8   = z6 * z2
    cdef double complex z10  = z8 * z2
    cdef double l    = <double> degree_l
    cdef double l_3  = 3. + 2. * l
    cdef double l_5  = 5. + 2. * l
    cdef double l_7  = 7. + 2. * l
    cdef double l_9  = 9. + 2. * l
    cdef double l_11 = 11. + 2. * l
    cdef double l_13 = 13. + 2. * l

    # Series expansion was done in `Takeuchi Starting Conditions.nb` on 2024-01-31 by JPR.
    phi_ptr[0] = (
         1. + 
        -z2  / (2.    * l_3) + 
         z4  / (8.    * l_3 * l_5) + 
        -z6  / (48.   * l_3 * l_5 * l_7) + 
         z8  / (384.  * l_3 * l_5 * l_7 * l_9) + 
        -z10 / (3840. * l_3 * l_5 * l_7 * l_9 * l_11)
    )
    
    # Note that the l+1 function has skips a denominator term than the function above.
    phi_lplus1_ptr[0] = (
         1. + 
        -z2  / (2.    * l_5) + 
         z4  / (8.    * l_5 * l_7) + 
        -z6  / (48.   * l_5 * l_7 * l_9) + 
         z8  / (384.  * l_5 * l_7 * l_9 * l_11) + 
        -z10 / (3840. * l_5 * l_7 * l_9 * l_11 * l_13)
    )
    
    psi_ptr[0] = (
         1. + 
        -z2  / (4.     * l_5) + 
        # RECORD: Error? TS72 quotes the next term as z4 / (12. * (5. + 2. * l) * (7 + 2. * l)); factor of 2 off in denom.
         z4  / (24.    * l_5 * l_7) + 
        -z6  / (192.   * l_5 * l_7 * l_9) + 
         z8  / (1920.  * l_5 * l_7 * l_9 * l_11) +
        -z10 / (23040. * l_5 * l_7 * l_9 * l_11 * l_13)
     )
    
    # Full version
    # cdef char lp1 = degree_l + 1
    # cdef double l_dbl_factorial   = cf_double_factorial(2 * degree_l + 1)
    # cdef double lp1_dbl_factorial = cf_double_factorial(2 * lp1 + 1)
    
    # cdef double complex z    = cf_csqrt(z2)
    # cdef double complex zl   = cf_cipow(z, degree_l)
    # cdef double complex zlp1 = cf_cipow(z, lp1)

    # phi_ptr[0]        = l_dbl_factorial * (spherical_jn(degree_l, z) / zl)
    # phi_lplus1_ptr[0] = lp1_dbl_factorial * (spherical_jn(lp1, z) / zlp1)
    # psi_ptr[0]        = (2. * (2. * degree_l + 3.) / (z * z)) * (1. - phi_ptr[0])

    # DEBUG: Uncomment the below and comment the above to use the method that TidalPy v0.4.0 used for this function.
    # Left here for comparison/debug purpouses.
    # cdef double degree_lp1 = degree_l + 1.
    # cdef double complex z4 = z2 * z2

    # phi_ptr[0] = 1. - \
    #       z2 / (2. * (2. * degree_l + 3.)) + \
    #       z4 / (8. * (2. * degree_l + 3.) * (2. * degree_l + 5.))
    # phi_lplus1_ptr[0] = 1. - \
    #              z2 / (2. * (2. * degree_lp1 + 3.)) + \
    #              z4 / (8. * (2. * degree_lp1 + 3.) * (2. * degree_lp1 + 5.))
    # psi_ptr[0] = 1. - \
    #       z2 / (4. * (2. * degree_l + 5.)) + \
    #       z4 / (12. * (2. * degree_l + 5.) * (2. * degree_l + 7.))
