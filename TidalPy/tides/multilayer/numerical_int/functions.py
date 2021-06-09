from typing import Union, Tuple

from scipy.special import gamma, spherical_jn

from ....utilities.types import FloatArray, ComplexArray
from ....utilities.performance import njit

FltCmplxArray = Union[FloatArray, ComplexArray]

# Pre-calculate as much as we can
l2p1_double_factorials = list()
for l in range(25):
    factorial_one = gamma(2. * l + 1. + 1.)

    dbl_factorial = factorial_one / ((2.**l) * gamma(l + 1.))
    l2p1_double_factorials.append(dbl_factorial)

l2p1_double_factorials = tuple(l2p1_double_factorials)

# OPT: This can not be njited because it depends on the non-njited spherical bessel functions.
#    subsequent functions that depend on this must use the approximate `takeuchi_phi_psi`.
def takeuchi_phi_psi_general(z : FltCmplxArray, order_l: int = 2) -> Tuple[FltCmplxArray, FltCmplxArray, FltCmplxArray]:
    """ Calculate the two functions used to find initial conditions for shooting method.

    See Eq. 103 in Takeuchi & Saito 1972

    Parameters
    ----------
    z : FltCmplxArray
        Input float or array.
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    phi : FltCmplxArray
    phi_lplus1 : FltCmplxArray
    psi : FltCmplxArray

    """

    phi = l2p1_double_factorials[order_l] * (spherical_jn(order_l, z) / z**order_l)
    phi_lplus1 = l2p1_double_factorials[order_l + 1] * (spherical_jn(order_l + 1, z) / z**(order_l + 1.))
    psi = (2. * (2. * order_l + 3.) / (z * z)) * (1. - phi)

    return phi, phi_lplus1, psi

@njit(cacheable=True)
def takeuchi_phi_psi(z2 : FltCmplxArray, order_l: int = 2) -> Tuple[FltCmplxArray, FltCmplxArray, FltCmplxArray]:
    """ Calculate the two functions used to find initial conditions for shooting method.

    See Eq. 103 in Takeuchi & Saito 1972

    Parameters
    ----------
    z2 : FltCmplxArray
        Input float or array. squared.
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    phi : FltCmplxArray
    phi_lplus1 : FltCmplxArray
    psi : FltCmplxArray

    """
    order_lp1 = order_l + 1.

    phi = 1. - \
          z2 / (2. * (2. * order_l + 3.)) + \
          (z2 * z2) / (8. * (2. * order_l + 3.) * (2. * order_l + 5.))
    phi_lplus1 = 1. - \
                 z2 / (2. * (2. * order_lp1 + 3.)) + \
                 (z2 * z2) / (8. * (2. * order_lp1 + 3.) * (2. * order_lp1 + 5.))
    psi = 1. - \
          z2 / (4. * (2. * order_lp1 + 5.)) + \
          (z2 * z2) / (12. * (2. * order_lp1 + 5.) * (2. * order_lp1 + 7.))

    return phi, phi_lplus1, psi