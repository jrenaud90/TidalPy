""" Functions to help calculate initial guess for radial functions at the bottom of a solid or liquid layer

References
----------
KMN15 : Kamata+ (2015; JGR-P; DOI: 10.1002/2015JE004821)
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

from typing import Tuple

import numpy as np
from scipy.special import gamma, spherical_jn

from .....utilities.performance import njit
from .....utilities.types import NumArray

# Pre-calculate as much as we can
l2p1_double_factorials = list()
for l in range(25):
    factorial_one = gamma(2. * l + 1. + 1.)

    dbl_factorial = factorial_one / ((2.**l) * gamma(l + 1.))
    l2p1_double_factorials.append(dbl_factorial)

l2p1_double_factorials = tuple(l2p1_double_factorials)


# OPT: This can not be njited because it depends on the non-njited spherical bessel functions.
#    subsequent functions that depend on this must use the approximate `takeuchi_phi_psi`.
def takeuchi_phi_psi_general(z: NumArray, order_l: int = 2) -> Tuple[NumArray, NumArray, NumArray]:
    """ Calculate the two (plus one) functions used to find initial conditions for shooting method.

    References
    ----------
    TS72 Eq. 103

    Parameters
    ----------
    z : NumArray
        Input scalar or array.
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    phi : NumArray
        Phi function.
    phi_lplus1 : NumArray
        Phi function for l+1.
    psi : NumArray
        Psi function.

    """

    phi = l2p1_double_factorials[order_l] * (spherical_jn(order_l, z) / z**order_l)
    phi_lplus1 = l2p1_double_factorials[order_l + 1] * (spherical_jn(order_l + 1, z) / z**(order_l + 1.))
    psi = (2. * (2. * order_l + 3.) / (z * z)) * (1. - phi)

    return phi, phi_lplus1, psi


@njit(cacheable=True)
def takeuchi_phi_psi(z_squared: NumArray, order_l: int = 2) -> Tuple[NumArray, NumArray, NumArray]:
    """ Calculate the two (plus one) functions used to find initial conditions for shooting method.

    This version of the function uses a Taylor expansion on the bessel function and only requires the even powers of
        z. Thus only z^2 is provided.

    References
    ----------
    TS72 Eq. 103

    Parameters
    ----------
    z_squared : NumArray
        Input scalar or array. squared.
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    phi : NumArray
        Phi function.
    phi_lplus1 : NumArray
        Phi function for l+1.
    psi : NumArray
        Psi function.

    """

    # Optimizations
    order_lp1 = order_l + 1.
    z_fourth = z_squared * z_squared

    phi = 1. - \
          z_squared / (2. * (2. * order_l + 3.)) + \
          z_fourth / (8. * (2. * order_l + 3.) * (2. * order_l + 5.))
    phi_lplus1 = 1. - \
                 z_squared / (2. * (2. * order_lp1 + 3.)) + \
                 z_fourth / (8. * (2. * order_lp1 + 3.) * (2. * order_lp1 + 5.))
    psi = 1. - \
          z_squared / (4. * (2. * order_l + 5.)) + \
          z_fourth / (12. * (2. * order_l + 5.) * (2. * order_l + 7.))

    return phi, phi_lplus1, psi


# Table was updated March 2022 by testing positive inputs from 0.1 to 1e9 and ensuring convergence to within
#  0.001% (0.00001)
Z_CALC_MAX_L = (
    (0.1, 7),
    (1., 7),
    (10., 8),
    (100., 17),
    (500., 32),
    (1000., 44),  # 1.23us
    (5000., 88),
    (10000., 117),  # 1.75us
    (50000., 245),
    (100000., 340),  # 2.87us
    (500000., 737),
    (1000000., 1034),  # 6.09us
    (5000000., 2279),
    # After this point these are taking a several seconds to calculate.
    (10000000., 3212),  # 15.7us
    (50000000., 7132),
    (100000000., 10071),  # 46us
    (500000000., 22450),
    (1000000000., 31721)  # 140us
    )


@njit(cacheable=True)
def z_calc(x_squared: NumArray, order_l: int = 2, init_l: int = 0, raise_l_error: bool = True) -> NumArray:
    """ Calculates the z function used in the calculations of initial guesses for radial functions.
    Simplification (recursion calculation) of the spherical Bessel function, see Eq. B16 of KMN15.

    OPT: This function is actually faster if a (too) large init_l is provided vs. it trying to find it on its own
        but that comes with a big risk that you have the wrong answer.
        Look at the table above for approx run times when x_squared is a float.

    References
    ----------
    TS72 Eqs. 96, 97
    KMN15 Eq. B16

    Parameters
    ----------
    x_squared : NumArray
        Expression passed to the Bessel function.
    order_l : int = 2
        Tidal harmonic order.
    init_l : int = 0
        Max integer to start the calculation from.
        If set to 0 then the function will try to determine the smallest value that still provides convergence.
    raise_l_error : bool = True
        If `True` then an exception will be raised if max l appears too small for convergence.

    Returns
    -------
    z : NumArray
        Result of the recursive calculation

    """

    max_l_might_be_too_small = False
    if init_l == 0:
        # The convergence of this function depends on the absolute size of the real part of x^2.
        # The table provided outside this function (`Z_CALC_MAX_L`) was tested on 2021/12/02 to find the smallest
        #    max l that still allowed convergence. It is not comprehensive, thus the possibility of an error being
        #    thrown below.
        x2_real = np.abs(np.real(x_squared))
        for min_val, max_l in Z_CALC_MAX_L:
            if np.all(np.asarray(x2_real <= min_val)):
                max_l_to_use = max_l + (order_l - 2)
                break
        else:
            # The value is outside of the predefined ranges. Go with something large and cross your fingers.
            max_l_to_use = 1000
            max_l_might_be_too_small = True
    else:
        max_l_to_use = max(order_l + 3, init_l)

    if max_l_might_be_too_small and raise_l_error:
        # Max l may be too small for convergence on large values of x_squared.
        #   Perform tests and then you can set raise_l_error=False to avoid this error.
        raise Exception

    z = x_squared / (2. * max_l_to_use + 3.)
    for l_fake in range(order_l, max_l_to_use + 1):
        l = max_l_to_use - l_fake + order_l
        # OPT: The above range is nicely written as range(order_l, init_l)[::-1]; but njit does not support this atm.
        z = x_squared / ((2. * l + 1.) - z)

    # print(max_l_to_use)

    return z
