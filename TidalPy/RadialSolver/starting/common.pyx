# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport abs

from scipy.special.cython_special cimport spherical_jn

from TidalPy.utilities.math.complex cimport cf_csqrt, cf_cipow, cf_cabs
from TidalPy.utilities.math.special cimport cf_double_factorial

cdef double complex cf_z_calc(
        double complex x_squared,
        int degree_l
        ) noexcept nogil:
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

    # QUESTION: (Issue #42) The recrusision formula shown in TS72 and KMN15 (See TS72 Eq. 97) does not reproduce the full version
    # of this function (TS72 Eq. 96). I am finding it to be off by about 30% at l = 2 and then the error improves to 
    # 10% at l=10. Perhaps thre is a problem with the implementation.
    # In any case, we will not use the recrusive formula here.

    # Have found that the Taylor expansion works well when x_square is small and is faster.
    cdef double l2            = <double>degree_l * 2.0
    
    # TODO: There is an issue with scipy 1.14.X+ on MacOS where the conversion between Py_ssize_t and long is not being accepted. See issue #65
    cdef long degree_l_ssizet = <long>degree_l

    cdef double complex x, z
    cdef double l2_3, l2_5, l2_7, l2_9, l2_11

    if cf_cabs(x_squared) > 0.1:
        # Use real function
        x = cf_csqrt(x_squared)
        z =  x * spherical_jn(degree_l_ssizet + 1, x) / spherical_jn(degree_l_ssizet, x)
    else:
        # Use Taylor series; JPR derived this on 2024-02-05
        l2_3  = l2 + 3.0
        l2_5  = l2 + 5.0
        l2_7  = l2 + 7.0
        l2_9  = l2 + 9.0
        l2_11 = l2 + 11.0
        #   Numerator                              | Denominator
        z = \
            (x_squared                             / l2_3) + \
            (x_squared**2                          / (l2_3**2 * l2_5)) + \
            (x_squared**4 * 2.                     / (l2_3**3 * l2_5 * l2_7)) + \
            (x_squared**6 * (27. + 10. * degree_l) / (l2_3**4 * l2_5**2 * l2_7 * l2_9)) + \
            (x_squared**8 * (90. + 28. * degree_l) / (l2_3**5 * l2_5**2 * l2_7 * l2_9 * l2_11))
    
    return z


cdef void cf_takeuchi_phi_psi(
        double complex z2,
        int degree_l,
        double complex* phi_ptr,
        double complex* phi_lplus1_ptr,
        double complex* psi_ptr,
        ) noexcept nogil:
    """ Calculate the two (plus one) functions used to find initial conditions for shooting method.

    References
    ----------
    TS72 Eq. 103

    """

    # Floating point errors prevent us from using the exact definition of these functions (See Issue #41)
    # Instead we use the limiting version of the functions.
    # However, we leave the full definition at the bottom of this function for reference.

    cdef double complex z4   = z2 * z2
    cdef double complex z6   = z4 * z2
    cdef double complex z8   = z6 * z2
    cdef double complex z10  = z8 * z2
    cdef double l    = <double>degree_l
    cdef double l2   = 2. * l
    cdef double l_3  = 3.  + l2
    cdef double l_5  = 5.  + l2
    cdef double l_7  = 7.  + l2
    cdef double l_9  = 9.  + l2
    cdef double l_11 = 11. + l2
    cdef double l_13 = 13. + l2

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
