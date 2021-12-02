from typing import Tuple, Union

from scipy.special import gamma, spherical_jn

from ....utilities.performance import njit
from ....utilities.types import ComplexArray, FloatArray

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
def takeuchi_phi_psi_general(z: FltCmplxArray, order_l: int = 2) -> Tuple[FltCmplxArray, FltCmplxArray, FltCmplxArray]:
    """ Calculate the two (plus one) functions used to find initial conditions for shooting method.

    See Eq. 103 in Takeuchi & Saito 1972

    Parameters
    ----------
    z : FltCmplxArray
        Input scalar or array.
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    phi : FltCmplxArray
        Phi function.
    phi_lplus1 : FltCmplxArray
        Phi function for l+1.
    psi : FltCmplxArray
        Psi function.

    """

    phi = l2p1_double_factorials[order_l] * (spherical_jn(order_l, z) / z**order_l)
    phi_lplus1 = l2p1_double_factorials[order_l + 1] * (spherical_jn(order_l + 1, z) / z**(order_l + 1.))
    psi = (2. * (2. * order_l + 3.) / (z * z)) * (1. - phi)

    return phi, phi_lplus1, psi


@njit(cacheable=True)
def takeuchi_phi_psi(z_squared: FltCmplxArray, order_l: int = 2) -> Tuple[FltCmplxArray, FltCmplxArray, FltCmplxArray]:
    """ Calculate the two (plus one) functions used to find initial conditions for shooting method.

    This version of the function uses a Taylor expansion on the bessel function and only requires the even powers of
        z. Thus only z^2 is provided.

    See Eq. 103 in Takeuchi & Saito 1972

    Parameters
    ----------
    z_squared : FltCmplxArray
        Input scalar or array. squared.
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    phi : FltCmplxArray
        Phi function.
    phi_lplus1 : FltCmplxArray
        Phi function for l+1.
    psi : FltCmplxArray
        Psi function.

    """
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
