""" Fundamental Matrix and its inverse as defined in Sabadini, Vermeerson, & Cambiotti (2016) (Hereafter cited as SVC16)

Two versions are defined, one for a generic order-l and one for order-l=2 (for performance).

Assumptions
-----------
These matrices assume an incompressible body.

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
"""

import numpy as np

from ....constants import pi, G
from ....utilities.performance import njit


@njit(cacheable=True)
def fundamental_matrix_generic(radius_array: np.ndarray, complex_shear_array: np.ndarray,
                               density_array: np.ndarray, gravity_array: np.ndarray, order_l: int = 2):
    """ Construct the fundamental matrix and its inverse for a generic order-l

    See Eq. 2.42 of SVC16

    Assumptions
    -----------
    - These matrices assume an incompressible body.

    Parameters
    ----------
    radius_array : np.ndarray
        Radius array of the planet [m]
    complex_shear_array : np.ndarray
        Complex Shear modulus at each radii [Pa]
    density_array : np.ndarray
        Density at each radii [kg m-3]
    gravity_array : np.ndarray
        Acceleration due to gravity at each radii [m s-2]
    order_l : int = 2
        Tidal Harmonic Degree Index

    Returns
    -------
    fundamental_matrix : np.ndarray
        Fundamental matrix used in the propagation technique
    fundamental_matrix_inverse : np.ndarray
        The inverse of the fundamental matrix used in the propagation technique
    derivative_mtx : np.ndarray
        The matrix, A, that satisfies the equation dy/dr = A * y

    """
    
    # First index: Rows of the propagation matrix
    # Second index: Columns of the propagation matrix
    # Third index: shell index (set outside of function)
    num_shells = radius_array.shape[0]
    fundamental_mtx = np.zeros((6, 6, num_shells), dtype=np.complex128)

    # Optimizations
    r_inv = 1. / radius_array
    rl = radius_array**order_l
    rlp1 = radius_array**(order_l + 1)
    rlp2 = radius_array**(order_l + 2)
    rlp3 = radius_array**(order_l + 3)
    rnl = radius_array**(-order_l)
    rnlm2 = radius_array**(-order_l - 2)
    rlm1 = radius_array**(order_l - 1)
    dlm1 = (2. * order_l - 1.)
    l2p3lm1 = (order_l**2 + 3. * order_l - 1.)
    l2mlm3 = (order_l**2 - order_l - 3.)
    lp1 = (order_l + 1.)
    lp2 = (order_l + 2.)
    lp3 = (order_l + 3.)
    l2m1 = (order_l**2 - 1.)
    dlp1 = (2. * order_l + 1.)
    dlp3 = (2. * order_l + 3.)
    
    rgp = radius_array * gravity_array * density_array
    rgp_s = rgp / complex_shear_array
    r_s = radius_array / complex_shear_array
    pr_s = density_array * r_s

    # Build Fundamental Matrix (zeros do not need to be specifically stated as they were put in at initialization)
    #     Eq. 2.42 in SVC
    ## Column 1
    fundamental_mtx[0, 0, :] = order_l * rlp1 / (2. * dlp3)
    fundamental_mtx[1, 0, :] = lp3 * rlp1 / (2. * dlp3 * lp1)
    fundamental_mtx[2, 0, :] = (order_l * rgp + 2. * l2mlm3 * complex_shear_array) * rl / (2. * dlp3)  # Believe there is a typo in HH14, they have the r^l only on one term instead of both.
    fundamental_mtx[3, 0, :] = order_l * lp2 * complex_shear_array * rl / (dlp3 * lp1)
    # fundamental_mtx[4, 0, :] = np.zeros(num_shells)
    fundamental_mtx[5, 0, :] = 2. * pi * G * density_array * order_l * rlp1 / dlp3

    ## Column 2
    fundamental_mtx[0, 1, :] = rlm1
    fundamental_mtx[1, 1, :] = rlm1 / order_l
    fundamental_mtx[2, 1, :] = (rgp + 2. * (order_l - 1.) * complex_shear_array) * radius_array**(order_l - 2.)
    fundamental_mtx[3, 1, :] = 2 * (order_l - 1.) * complex_shear_array * radius_array**(order_l - 2.) / order_l
    # fundamental_mtx[4, 1, :] = np.zeros(num_shells)
    fundamental_mtx[5, 1, :] = 4. * pi * G * density_array * rlm1

    ## Column 3
    # fundamental_mtx[0, 2, :] = np.zeros(num_shells)
    # fundamental_mtx[1, 2, :] = np.zeros(num_shells)
    fundamental_mtx[2, 2, :] = -density_array * rl
    # fundamental_mtx[3, 2, :] = np.zeros(num_shells)
    fundamental_mtx[4, 2, :] = -rl
    fundamental_mtx[5, 2, :] = -dlp1 * rlm1

    ## Column 4
    fundamental_mtx[0, 3, :] = lp1 * rnl / (2. * dlm1)
    fundamental_mtx[1, 3, :] = (2. - order_l) * rnl / (2. * order_l * dlm1)
    fundamental_mtx[2, 3, :] = (lp1 * rgp - 2. * l2p3lm1 * complex_shear_array) / (2. * dlm1 * rlp1)
    fundamental_mtx[3, 3, :] = l2m1 * complex_shear_array / (order_l * dlm1 * rlp1)
    # fundamental_mtx[4, 3, :] = np.zeros(num_shells)
    fundamental_mtx[5, 3, :] = 2 * pi * G * density_array * lp1 / (dlm1 * rl)

    ## Column 5
    fundamental_mtx[0, 4, :] = rnlm2
    fundamental_mtx[1, 4, :] = -rnlm2 / lp1
    fundamental_mtx[2, 4, :] = (rgp - 2. * lp2 * complex_shear_array) / rlp3
    fundamental_mtx[3, 4, :] = 2. * lp2 * complex_shear_array / (lp1 * rlp3)
    # fundamental_mtx[4, 4, :] = np.zeros(num_shells)
    fundamental_mtx[5, 4, :] = 4. * pi * G * density_array / rlp2

    ## Column 6
    # fundamental_mtx[0, 5, :] = np.zeros(num_shells)
    # fundamental_mtx[1, 5, :] = np.zeros(num_shells)
    fundamental_mtx[2, 5, :] = -density_array / rlp1
    # fundamental_mtx[3, 5, :] = np.zeros(num_shells)
    fundamental_mtx[4, 5, :] = -1. / rlp1
    # fundamental_mtx[5, 5, :] = np.zeros(num_shells)
    
    # Inverse of the Fundamental Matrix
    # This function manually defines the inverse matrix which is about 15--30% faster than a version that uses
    #    np.linalg.inv() to calculate the inverse matrix.
    #
    # From SVC16 Eq. 2.45: Fundamental Inverse = D_Mtx * Y^Bar_Mtx
    # D_Mtx is a diagonal matrix with
    # 1/(2l+1) * [ ... ]
    # We are going to multiple first and just write down the fundamental matrix inverse to avoid the additional D*Ybar
    #    calculation.
    inverse_fundamental_mtx = np.zeros((6, 6, num_shells), dtype=np.complex128)

    # D Coefficients
    coeff = (1. / dlp1)
    d_coeff_1 = coeff * lp1 / rlp1
    d_coeff_2 = coeff * order_l * lp1 / (2. * dlm1 * rlm1)
    d_coeff_3 = coeff * 1. / rlm1
    d_coeff_4 = coeff * order_l * rl
    d_coeff_5 = coeff * rlp2 * order_l * lp1 / (2. * dlp3)
    d_coeff_6 = coeff * -rlp1

    ## Column 1
    inverse_fundamental_mtx[0, 0, :] = d_coeff_1 * (rgp_s - 2. * lp2)
    inverse_fundamental_mtx[1, 0, :] = d_coeff_2 * (-rgp_s + 2. * l2p3lm1 / lp1)
    inverse_fundamental_mtx[2, 0, :] = d_coeff_3 * (4. * pi * G * density_array)
    inverse_fundamental_mtx[3, 0, :] = d_coeff_4 * (rgp_s + 2. * (order_l - 1.))
    inverse_fundamental_mtx[4, 0, :] = d_coeff_5 * (-rgp_s - 2. * l2mlm3 / order_l)
    inverse_fundamental_mtx[5, 0, :] = d_coeff_6 * (4. * pi * G * density_array * radius_array)

    ## Column 2
    inverse_fundamental_mtx[0, 1, :] = d_coeff_1 * (2. * order_l * lp2)
    inverse_fundamental_mtx[1, 1, :] = d_coeff_2 * (-2. * l2m1)
    # inverse_fundamental_mtx[2, 1, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[3, 1, :] = d_coeff_4 * (2. * l2m1)
    inverse_fundamental_mtx[4, 1, :] = d_coeff_5 * (-2. * order_l * lp2)
    # inverse_fundamental_mtx[5, 1, :] = np.zeros(num_shells)

    ## Column 3
    inverse_fundamental_mtx[0, 2, :] = d_coeff_1 * (-r_s)
    inverse_fundamental_mtx[1, 2, :] = d_coeff_2 * (r_s)
    # inverse_fundamental_mtx[2, 2, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[3, 2, :] = d_coeff_4 * (-r_s)
    inverse_fundamental_mtx[4, 2, :] = d_coeff_5 * (r_s)
    # inverse_fundamental_mtx[5, 2, :] = np.zeros(num_shells)

    ## Column 4
    inverse_fundamental_mtx[0, 3, :] = d_coeff_1 * (order_l * r_s)
    inverse_fundamental_mtx[1, 3, :] = d_coeff_2 * ((2. - order_l) * r_s)
    # inverse_fundamental_mtx[2, 3, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[3, 3, :] = d_coeff_4 * (-lp1 * r_s)
    inverse_fundamental_mtx[4, 3, :] = d_coeff_5 * (lp3 * r_s)
    # inverse_fundamental_mtx[5, 3, :] = np.zeros(num_shells)

    ## Column 5
    inverse_fundamental_mtx[0, 4, :] = d_coeff_1 * (pr_s)
    inverse_fundamental_mtx[1, 4, :] = d_coeff_2 * (-pr_s)
    # inverse_fundamental_mtx[2, 4, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[3, 4, :] = d_coeff_4 * (pr_s)
    inverse_fundamental_mtx[4, 4, :] = d_coeff_5 * (-pr_s)
    inverse_fundamental_mtx[5, 4, :] = d_coeff_6 * dlp1

    ## Column 6
    # inverse_fundamental_mtx[0, 5, :] = np.zeros(num_shells)
    # inverse_fundamental_mtx[1, 5, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[2, 5, :] =  d_coeff_3 * (-np.ones(num_shells, dtype=np.complex128))
    # inverse_fundamental_mtx[3, 5, :] = np.zeros(num_shells)
    # inverse_fundamental_mtx[4, 5, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[5, 5, :] = d_coeff_6 * (-radius_array)

    # Build derivative matrix
    # Defined in SV04 -- Only valid for the incompressible case.
    # See SVC16 Eq. 1.95
    #    Note: the lambda in SVC16 is defined as bulk_mod - (2. / 3.) * shear (Eq. 1.77; 2nd Lame parameter),
    #    for the incompressible assumption we will assume the ratio that SVC16 use (lambda / beta) -> 1 as K -> inf
    #    See SVC16 Eq. 1.95 for a compressible version. Take limit as K->inf to find below.
    derivative_mtx = np.zeros((6, 6, num_shells), dtype=np.complex128)
    ## Column 1
    derivative_mtx[0, 0, :] = -2. * r_inv
    derivative_mtx[1, 0, :] = -1. * r_inv
    derivative_mtx[2, 0, :] = (4. * r_inv) * (3. * complex_shear_array * r_inv - density_array * gravity_array)
    derivative_mtx[3, 0, :] = (-1. * r_inv) * (6. * complex_shear_array * r_inv - density_array * gravity_array)
    derivative_mtx[4, 0, :] = -4. * np.pi * G * density_array
    derivative_mtx[5, 0, :] = -4. * np.pi * G * density_array * (order_l + 1) * r_inv

    ## Column 2
    derivative_mtx[0, 1, :] = order_l * lp1 * r_inv
    derivative_mtx[1, 1, :] = r_inv
    derivative_mtx[2, 1, :] = -order_l * lp1 * r_inv * (6. * complex_shear_array * r_inv -
                                                                  density_array * gravity_array)
    derivative_mtx[3, 1, :] = 2. * (2. * order_l**2 + 2. * order_l - 1.) * complex_shear_array * (r_inv**2)
    # derivative_mtx[4, 1, :] = np.zeros(num_shells)
    derivative_mtx[5, 1, :] = 4. * np.pi * G * density_array * order_l * lp1 * r_inv

    ## Column 3
    # derivative_mtx[0, 2, :] = np.zeros(num_shells)
    # derivative_mtx[1, 2, :] = np.zeros(num_shells)
    # derivative_mtx[2, 2, :] = np.zeros(num_shells)
    derivative_mtx[3, 2, :] = -r_inv
    # derivative_mtx[4, 2, :] = np.zeros(num_shells)
    # derivative_mtx[5, 2, :] = np.zeros(num_shells)

    ## Column 4
    # derivative_mtx[0, 3, :] = np.zeros(num_shells)
    derivative_mtx[1, 3, :] = 1. / complex_shear_array
    derivative_mtx[2, 3, :] = order_l * lp1 * r_inv
    derivative_mtx[3, 3, :] = -3. * r_inv
    # derivative_mtx[4, 3, :]= np.zeros(num_shells)
    # derivative_mtx[5, 3, :] = np.zeros(num_shells)

    ## Column 5
    # derivative_mtx[0, 4, :] = np.zeros(num_shells)
    # derivative_mtx[1, 4, :] = np.zeros(num_shells)
    derivative_mtx[2, 4, :] = -density_array * lp1 * r_inv
    derivative_mtx[3, 4, :] = density_array * r_inv
    derivative_mtx[4, 4, :] = -lp1 * r_inv
    # derivative_mtx[5, 4, :] = np.zeros(num_shells)

    ## Column 6
    # derivative_mtx[0, 5, :] = np.zeros(num_shells)
    # derivative_mtx[1, 5, :] = np.zeros(num_shells)
    derivative_mtx[2, 5, :] = density_array
    # derivative_mtx[3, 5, :] = np.zeros(num_shells)
    derivative_mtx[4, 5, :] = np.ones(num_shells, dtype=np.complex128)
    derivative_mtx[5, 5, :] = (order_l - 1.) * r_inv

    return fundamental_mtx, inverse_fundamental_mtx, derivative_mtx

@njit(cacheable=True)
def fundamental_matrix_orderl2(radius_array: np.ndarray, complex_shear_array: np.ndarray,
                               density_array: np.ndarray, gravity_array: np.ndarray):
    """ Construct the fundamental matrix and its inverse for a order-l = 2

    This is a restricted version of the generic fundamental matrix that is only valid for tidal order l = 2. The
        purpose in using this over the fundamental is an increase in performance.

    See Eq. 2.42 of SVC16

    Compare fundamental matrix to Eq. A4 of HH14 and ID variable "Ypropmtx"

    Assumptions
    -----------
    - These matrices assume an incompressible body.
    - This function is restricted to order-l = 2

    Parameters
    ----------
    radius_array : np.ndarray
        Radius array of the planet [m]
    complex_shear_array : np.ndarray
        Complex Shear modulus at each radii [Pa]
    density_array : np.ndarray
        Density at each radii [kg m-3]
    gravity_array : np.ndarray
        Acceleration due to gravity at each radii [m s-2]

    Returns
    -------
    fundamental_matrix : np.ndarray
        Fundamental matrix used in the propagation technique
    fundamental_matrix_inverse : np.ndarray
        The inverse of the fundamental matrix used in the propagation technique
    derivative_mtx : np.ndarray
        The matrix, A, that satisfies the equation dy/dr = A * y

    See Also
    --------
    TidalPy.tides.multilayer.fundamental.fundamental_matrix_generic

    """

    # First index: Rows of the propagation matrix
    # Second index: Columns of the propagation matrix
    # Third index: shell index (set outside of function)
    num_shells = radius_array.shape[0]
    fundamental_mtx = np.zeros((6, 6, num_shells), dtype=np.complex128)

    # Optimizations
    r_inv = 1. / radius_array
    rgp = radius_array * gravity_array * density_array
    rgp_s = rgp / complex_shear_array
    r_s = radius_array / complex_shear_array
    pr_s = density_array * r_s
    r2 = radius_array * radius_array
    r3 = radius_array * r2
    r4 = radius_array * r3
    r5 = radius_array * r4

    # Build Fundamental Matrix (zeros do not need to be specifically stated as they were put in at initialization)
    #     Eq. 2.42 in SVC
    ## Column 1
    fundamental_mtx[0, 0, :] = r3 / 7.
    fundamental_mtx[1, 0, :] = 5. * r3 / 42.
    fundamental_mtx[2, 0, :] = (rgp - complex_shear_array) * r2 / 7.
    fundamental_mtx[3, 0, :] = 8. * complex_shear_array * r2 / 21.
    # fundamental_mtx[4, 0, :] = np.zeros(num_shells)
    fundamental_mtx[5, 0, :] = 4. * pi * G * density_array * r3 / 7.

    ## Column 2
    fundamental_mtx[0, 1, :] = radius_array
    fundamental_mtx[1, 1, :] = radius_array / 2.
    fundamental_mtx[2, 1, :] = rgp + 2. * complex_shear_array
    fundamental_mtx[3, 1, :] = complex_shear_array
    # fundamental_mtx[4, 1, :] = np.zeros(num_shells)
    fundamental_mtx[5, 1, :] = 4. * pi * G * density_array * radius_array

    ## Column 3
    # fundamental_mtx[0, 2, :] = np.zeros(num_shells)
    # fundamental_mtx[1, 2, :] = np.zeros(num_shells)
    fundamental_mtx[2, 2, :] = -density_array * r2
    # fundamental_mtx[3, 2, :] = np.zeros(num_shells)
    fundamental_mtx[4, 2, :] = -r2
    fundamental_mtx[5, 2, :] = -5. * radius_array

    ## Column 4
    fundamental_mtx[0, 3, :] = 1. / (2. * r2)
    # fundamental_mtx[1, 3, :] = np.zeros(num_shells)
    fundamental_mtx[2, 3, :] = (rgp - 6. * complex_shear_array) / (2. * r3)
    fundamental_mtx[3, 3, :] = complex_shear_array / (2. * r3)
    # fundamental_mtx[4, 3, :] = np.zeros(num_shells)
    fundamental_mtx[5, 3, :] = 2. * pi * G * density_array / r2

    ## Column 5
    fundamental_mtx[0, 4, :] = 1. / r4
    fundamental_mtx[1, 4, :] = -1. / (3. * r4)
    fundamental_mtx[2, 4, :] = (rgp - 8. * complex_shear_array) / r5
    fundamental_mtx[3, 4, :] = 8. * complex_shear_array / (3. * r5)
    # fundamental_mtx[4, 4, :] = np.zeros(num_shells)
    fundamental_mtx[5, 4, :] = 4. * pi * G * density_array / r4

    ## Column 6
    # fundamental_mtx[0, 5, :] = np.zeros(num_shells)
    # fundamental_mtx[1, 5, :] = np.zeros(num_shells)
    fundamental_mtx[2, 5, :] = -density_array / r3
    # fundamental_mtx[3, 5, :] = np.zeros(num_shells)
    fundamental_mtx[4, 5, :] = -1. / r3
    # fundamental_mtx[5, 5, :] = np.zeros(num_shells)

    # Inverse of the Fundamental Matrix
    # This function manually defines the inverse matrix which is about 15--30% faster than a version that uses
    #    np.linalg.inv() to calculate the inverse matrix.
    #
    # From SVC16 Eq. 2.45: Fundamental Inverse = D_Mtx * Y^Bar_Mtx
    # D_Mtx is a diagonal matrix with
    # 1/(2l+1) * [ ... ]
    # We are going to multiple first and just write down the fundamental matrix inverse to avoid the additional D*Ybar
    #    calculation.
    inverse_fundamental_mtx = np.zeros((6, 6, num_shells), dtype=np.complex128)

    # D Coefficients
    d_coeff_1 = 3. / (5. * r3)
    d_coeff_2 = 1. / (5. * radius_array)
    d_coeff_3 = 1. / (5. * radius_array)
    d_coeff_4 = (2. / 5.) * r2
    d_coeff_5 = (3. / 35.) * r4
    d_coeff_6 = (1. / 5.) * -r3

    ## Column 1
    inverse_fundamental_mtx[0, 0, :] = d_coeff_1 * (rgp_s - 8.)
    inverse_fundamental_mtx[1, 0, :] = d_coeff_2 * (-rgp_s + 6.)
    inverse_fundamental_mtx[2, 0, :] = d_coeff_3 * (4. * pi * G * density_array)
    inverse_fundamental_mtx[3, 0, :] = d_coeff_4 * (rgp_s + 2.)
    inverse_fundamental_mtx[4, 0, :] = d_coeff_5 * (-rgp_s + 1.)
    inverse_fundamental_mtx[5, 0, :] = d_coeff_6 * (4. * pi * G * density_array * radius_array)

    ## Column 2
    inverse_fundamental_mtx[0, 1, :] = d_coeff_1 * 16.
    inverse_fundamental_mtx[1, 1, :] = d_coeff_2 * (-6.)
    # inverse_fundamental_mtx[2, 1, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[3, 1, :] = d_coeff_4 * 6.
    inverse_fundamental_mtx[4, 1, :] = d_coeff_5 * (-16.)
    # inverse_fundamental_mtx[5, 1, :] = np.zeros(num_shells)

    ## Column 3
    inverse_fundamental_mtx[0, 2, :] = d_coeff_1 * (-r_s)
    inverse_fundamental_mtx[1, 2, :] = d_coeff_2 * (r_s)
    # inverse_fundamental_mtx[2, 2, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[3, 2, :] = d_coeff_4 * (-r_s)
    inverse_fundamental_mtx[4, 2, :] = d_coeff_5 * (r_s)
    # inverse_fundamental_mtx[5, 2, :] = np.zeros(num_shells)

    ## Column 4
    inverse_fundamental_mtx[0, 3, :] = d_coeff_1 * (2. * r_s)
    # inverse_fundamental_mtx[1, 3, :] = np.zeros(num_shells)
    # inverse_fundamental_mtx[2, 3, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[3, 3, :] = d_coeff_4 * (-3. * r_s)
    inverse_fundamental_mtx[4, 3, :] = d_coeff_5 * (5. * r_s)
    # inverse_fundamental_mtx[5, 3, :] = np.zeros(num_shells)

    ## Column 5
    inverse_fundamental_mtx[0, 4, :] = d_coeff_1 * (pr_s)
    inverse_fundamental_mtx[1, 4, :] = d_coeff_2 * (-pr_s)
    # inverse_fundamental_mtx[2, 4, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[3, 4, :] = d_coeff_4 * (pr_s)
    inverse_fundamental_mtx[4, 4, :] = d_coeff_5 * (-pr_s)
    inverse_fundamental_mtx[5, 4, :] = d_coeff_6 * 5.

    ## Column 6
    # inverse_fundamental_mtx[0, 5, :] = np.zeros(num_shells)
    # inverse_fundamental_mtx[1, 5, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[2, 5, :] = d_coeff_3 * (-np.ones(num_shells, dtype=np.complex128))
    # inverse_fundamental_mtx[3, 5, :] = np.zeros(num_shells)
    # inverse_fundamental_mtx[4, 5, :] = np.zeros(num_shells)
    inverse_fundamental_mtx[5, 5, :] = d_coeff_6 * (-radius_array)

    # Build derivative matrix
    # Defined in SV04 -- Only valid for the incompressible case.
    # See SVC16 Eq. 1.95
    #    Note: the lambda in SVC16 is defined as bulk_mod - (2. / 3.) * shear (Eq. 1.77; 2nd Lame parameter),
    #    for the incompressible assumption we will assume the ratio that SVC16 use (lambda / beta) -> 1 as K -> inf
    #    See SVC16 Eq. 1.95 for a compressible version. Take limit as K->inf to find below.
    derivative_mtx = np.zeros((6, 6, num_shells), dtype=np.complex128)
    ## Column 1
    derivative_mtx[0, 0, :] = -2. * r_inv
    derivative_mtx[1, 0, :] = -1. * r_inv
    derivative_mtx[2, 0, :] = (4. * r_inv) * (3. * complex_shear_array * r_inv - density_array * gravity_array)
    derivative_mtx[3, 0, :] = (-1. * r_inv) * (6. * complex_shear_array * r_inv - density_array * gravity_array)
    derivative_mtx[4, 0, :] = -4. * np.pi * G * density_array
    derivative_mtx[5, 0, :] = -12. * np.pi * G * density_array * r_inv

    ## Column 2
    derivative_mtx[0, 1, :] = 6. * r_inv
    derivative_mtx[1, 1, :] = r_inv
    derivative_mtx[2, 1, :] = -6. * r_inv * (6. * complex_shear_array * r_inv - density_array * gravity_array)
    derivative_mtx[3, 1, :] = 22. * complex_shear_array * (r_inv**2)
    # derivative_mtx[4, 1, :] = np.zeros(num_shells)
    derivative_mtx[5, 1, :] = 24. * np.pi * G * density_array * r_inv

    ## Column 3
    # derivative_mtx[0, 2, :] = np.zeros(num_shells)
    # derivative_mtx[1, 2, :] = np.zeros(num_shells)
    # derivative_mtx[2, 2, :] = np.zeros(num_shells)
    derivative_mtx[3, 2, :] = -r_inv
    # derivative_mtx[4, 2, :] = np.zeros(num_shells)
    # derivative_mtx[5, 2, :] = np.zeros(num_shells)

    ## Column 4
    # derivative_mtx[0, 3, :] = np.zeros(num_shells)
    derivative_mtx[1, 3, :] = 1. / complex_shear_array
    derivative_mtx[2, 3, :] = 6. * r_inv
    derivative_mtx[3, 3, :] = -3. * r_inv
    # derivative_mtx[4, 3, :]= np.zeros(num_shells)
    # derivative_mtx[5, 3, :] = np.zeros(num_shells)

    ## Column 5
    # derivative_mtx[0, 4, :] = np.zeros(num_shells)
    # derivative_mtx[1, 4, :] = np.zeros(num_shells)
    derivative_mtx[2, 4, :] = -3. * density_array * r_inv
    derivative_mtx[3, 4, :] = density_array * r_inv
    derivative_mtx[4, 4, :] = -3. * r_inv
    # derivative_mtx[5, 4, :] = np.zeros(num_shells)

    ## Column 6
    # derivative_mtx[0, 5, :] = np.zeros(num_shells)
    # derivative_mtx[1, 5, :] = np.zeros(num_shells)
    derivative_mtx[2, 5, :] = density_array
    # derivative_mtx[3, 5, :] = np.zeros(num_shells)
    derivative_mtx[4, 5, :] = np.ones(num_shells, dtype=np.complex128)
    derivative_mtx[5, 5, :] = r_inv

    return fundamental_mtx, inverse_fundamental_mtx, derivative_mtx
