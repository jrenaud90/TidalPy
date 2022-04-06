""" Functions to calculate the initial guess for radial functions at the bottom of a solid or liquid layer

These functions are the most general and allow for:
    - Compressibility
    - Dynamic (non-static) Tides
    - Liquid layer propagation

References
----------
KMN15 : Kamata+ (2015; JGR-P; DOI: 10.1002/2015JE004821)
B15   : Beuthe+ (2015; Icarus; DOI: 10.1016/j.icarus.2015.06.008)
TMS05 : Tobie+ (2005; Icarus; DOI: 10.1016/j.icarus.2005.04.006)
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

from typing import List

import numpy as np

from .functions import takeuchi_phi_psi, z_calc
from .....constants import G, pi
from .....utilities.math.special import sqrt_neg
from .....utilities.performance import njit, nbList
from .....utilities.types import ComplexArray, FloatArray, NumArray

SolidDynamicGuess = List[ComplexArray]
LiquidDynamicGuess = List[ComplexArray]


@njit(cacheable=True)
def solid_guess_kamata(
    radius: FloatArray, shear_modulus: NumArray, bulk_modulus: NumArray,
    density: FloatArray, frequency: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> SolidDynamicGuess:
    """ Calculate the initial guess at the bottom of a solid layer using the dynamic assumption.

    This function uses the Kamata et al (2015; JGR:P) equations (Eq. B1-B16).

    Using the dynamic assumption in a solid layer results in three independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
       bulk and shear dissipation.

    References
    ----------
    KMN15 Eqs. B1-B16

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    shear_modulus : NumArray
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : NumArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at  at `radius` [kg m-3]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : SolidDynamicGuess
        The three independent solid guesses (sn1, sn2, sn3)

    """

    # Convert compressibility parameters
    lame = bulk_modulus - (2. / 3.) * shear_modulus

    # Constants (See Eqs. B13-B16 of KMN15)
    dynamic_term = frequency * frequency
    alpha2 = (lame + 2. * shear_modulus) / density
    beta2 = shear_modulus / density
    gamma = 4. * pi * G_to_use * density / 3.

    # Optimizations
    r_inverse = 1. / radius
    r2_inverse = r_inverse * r_inverse
    r2 = radius * radius

    # Helper functions
    k2_quad_pos = (dynamic_term / beta2) + ((dynamic_term + 4. * gamma) / alpha2)
    k2_quad_neg = (dynamic_term / beta2) - ((dynamic_term + 4. * gamma) / alpha2)
    k2_quad = k2_quad_neg**2 + ((4. * order_l * (order_l + 1) * gamma**2) / (alpha2 * beta2))

    # TODO: TS74 (eq. 99) has these flipped compared to KMN15. Going with KMN for this func.
    #    [GitHub Issue](https://github.com/jrenaud90/TidalPy/issues/31)
    k2_pos = (1. / 2.) * (k2_quad_pos + sqrt_neg(k2_quad, is_real=True))
    k2_neg = (1. / 2.) * (k2_quad_pos - sqrt_neg(k2_quad, is_real=True))

    f_k2_pos = (beta2 * k2_pos - dynamic_term) / gamma
    f_k2_neg = (beta2 * k2_neg - dynamic_term) / gamma

    h_k2_pos = f_k2_pos - (order_l + 1.)
    h_k2_neg = f_k2_neg - (order_l + 1.)

    z_k2_pos = z_calc(k2_pos * r2, order_l=order_l, init_l=0, raise_l_error=True)
    z_k2_neg = z_calc(k2_neg * r2, order_l=order_l, init_l=0, raise_l_error=True)

    # See Eqs. B1-B12 of KMN15
    y1_s1 = -f_k2_pos * z_k2_pos * r_inverse
    y1_s2 = -f_k2_neg * z_k2_neg * r_inverse
    y1_s3 = order_l * r_inverse

    y2_s1 = -density * f_k2_pos * alpha2 * k2_pos + \
            (2. * shear_modulus * r2_inverse) * (2. * f_k2_pos + order_l * (order_l + 1.)) * z_k2_pos
    y2_s2 = -density * f_k2_neg * alpha2 * k2_neg + \
            (2. * shear_modulus * r2_inverse) * (2. * f_k2_neg + order_l * (order_l + 1.)) * z_k2_neg
    y2_s3 = 2. * shear_modulus * order_l * (order_l - 1) * r2_inverse

    y3_s1 = z_k2_pos * r_inverse
    y3_s2 = z_k2_neg * r_inverse
    y3_s3 = r_inverse

    y4_s1 = shear_modulus * k2_pos - (2. * shear_modulus * r2_inverse) * (f_k2_pos + 1.) * z_k2_pos
    y4_s2 = shear_modulus * k2_neg - (2. * shear_modulus * r2_inverse) * (f_k2_neg + 1.) * z_k2_neg
    y4_s3 = 2. * shear_modulus * (order_l - 1.) * r2_inverse

    y5_s1 = 3. * gamma * f_k2_pos - h_k2_pos * (order_l * gamma - dynamic_term)
    y5_s2 = 3. * gamma * f_k2_neg - h_k2_neg * (order_l * gamma - dynamic_term)
    y5_s3 = order_l * gamma - dynamic_term

    y6_s1 = (2. * order_l + 1.) * y5_s1 * r_inverse
    y6_s2 = (2. * order_l + 1.) * y5_s2 * r_inverse
    y6_s3 = (2. * order_l + 1.) * y5_s3 * r_inverse - (3. * order_l * gamma * r_inverse)

    # TODO: Right now numba does not support np.stack for purely scalar inputs. A temp fix is to make sure all the
    #    inputs are cast into arrays. See the github issue here: https://github.com/numba/numba/issues/7002
    y1_s1 = np.asarray(y1_s1, dtype=np.complex128)
    y2_s1 = np.asarray(y2_s1, dtype=np.complex128)
    y3_s1 = np.asarray(y3_s1, dtype=np.complex128)
    y4_s1 = np.asarray(y4_s1, dtype=np.complex128)
    y5_s1 = np.asarray(y5_s1, dtype=np.complex128)
    y6_s1 = np.asarray(y6_s1, dtype=np.complex128)

    y1_s2 = np.asarray(y1_s2, dtype=np.complex128)
    y2_s2 = np.asarray(y2_s2, dtype=np.complex128)
    y3_s2 = np.asarray(y3_s2, dtype=np.complex128)
    y4_s2 = np.asarray(y4_s2, dtype=np.complex128)
    y5_s2 = np.asarray(y5_s2, dtype=np.complex128)
    y6_s2 = np.asarray(y6_s2, dtype=np.complex128)

    y1_s3 = np.asarray(y1_s3, dtype=np.complex128)
    y2_s3 = np.asarray(y2_s3, dtype=np.complex128)
    y3_s3 = np.asarray(y3_s3, dtype=np.complex128)
    y4_s3 = np.asarray(y4_s3, dtype=np.complex128)
    y5_s3 = np.asarray(y5_s3, dtype=np.complex128)
    y6_s3 = np.asarray(y6_s3, dtype=np.complex128)

    # Combine the three solutions
    tidaly_s1 = np.stack((y1_s1, y2_s1, y3_s1, y4_s1, y5_s1, y6_s1))
    tidaly_s2 = np.stack((y1_s2, y2_s2, y3_s2, y4_s2, y5_s2, y6_s2))
    tidaly_s3 = np.stack((y1_s3, y2_s3, y3_s3, y4_s3, y5_s3, y6_s3))

    return nbList([tidaly_s1, tidaly_s2, tidaly_s3])


@njit(cacheable=True)
def liquid_guess_kamata(
    radius: FloatArray, bulk_modulus: NumArray,
    density: FloatArray, frequency: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> LiquidDynamicGuess:
    """  Calculate the initial guess at the bottom of a liquid layer using the dynamic assumption.

    This function uses the Kamata et al (2015; JGR:P) equations (Eq. B29-B37).

    Using the dynamic assumption in a liquid layer results in two independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
       bulk dissipation (no shear dissipation within liquid layers).

    References
    ----------
    KMN15 Eq. B29-B37

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    bulk_modulus : NumArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at  at `radius` [kg m-3]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : LiquidDynamicGuess
        The two independent liquid guesses (sn1, sn2)

    """

    # Convert compressibility parameters
    # For the liquid layer the shear modulus is zero so the 1st Lame parameter = bulk modulus
    lame = bulk_modulus

    # Optimizations
    dynamic_term = frequency * frequency
    r_inverse = 1. / radius
    r2 = radius * radius

    # Helper functions
    gamma = (4. * pi * G_to_use * density / 3.)
    alpha2 = lame / density
    k2 = (1. / alpha2) * (dynamic_term + 4. * gamma - (order_l * (order_l + 1) * gamma**2 / dynamic_term))
    f = -dynamic_term / gamma
    h = f - (order_l + 1.)

    # See Eqs. B33--B36 in KMN15
    y1_s1 = -f * r_inverse * z_calc(k2 * r2, order_l=order_l)
    y1_s2 = order_l * r_inverse

    y2_s1 = -density * (f * (dynamic_term + 4 * gamma) + order_l * (order_l + 1) * gamma)
    y2_s2 = 0. * radius

    y5_s1 = 3. * gamma * f - h * (order_l * gamma - dynamic_term)
    y5_s2 = order_l * gamma - dynamic_term

    y6_s1 = (2. * order_l + 1.) * y5_s1 * r_inverse
    y6_s2 = ((2. * order_l + 1.) * y5_s2 * r_inverse) - ((3. * order_l * gamma) * r_inverse)

    # TODO: Right now numba does not support np.stack for purely scalar inputs. A temp fix is to make sure all the
    #    inputs are cast into arrays. See the github issue here: https://github.com/numba/numba/issues/7002
    y1_s1 = np.asarray(y1_s1, dtype=np.complex128)
    y2_s1 = np.asarray(y2_s1, dtype=np.complex128)
    y5_s1 = np.asarray(y5_s1, dtype=np.complex128)
    y6_s1 = np.asarray(y6_s1, dtype=np.complex128)

    y1_s2 = np.asarray(y1_s2, dtype=np.complex128)
    y2_s2 = np.asarray(y2_s2, dtype=np.complex128)
    y5_s2 = np.asarray(y5_s2, dtype=np.complex128)
    y6_s2 = np.asarray(y6_s2, dtype=np.complex128)

    # Combine the two solutions
    tidaly_s1 = np.stack((y1_s1, y2_s1, y5_s1, y6_s1))
    tidaly_s2 = np.stack((y1_s2, y2_s2, y5_s2, y6_s2))

    return nbList([tidaly_s1, tidaly_s2])


@njit(cacheable=True)
def solid_guess_takeuchi(
    radius: FloatArray, shear_modulus: NumArray, bulk_modulus: NumArray,
    density: FloatArray, frequency: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> SolidDynamicGuess:
    """ Calculate the initial guess at the bottom of a solid layer using the dynamic assumption.

    This function uses the Takeuchi and Saito 1972 equations (Eq. 95-101).

    Using the dynamic assumption in a solid layer results in three independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
       bulk and shear dissipation.

    References
    ----------
    TS72

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    shear_modulus : NumArray
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : NumArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at  at `radius` [kg m-3]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : SolidDynamicGuess
        The three independent solid guesses (sn1, sn2, sn3)

    """

    # Convert compressibility parameters
    lame = bulk_modulus - (2. / 3.) * shear_modulus

    # Constants (See Eqs. B13-B16 of KMN15)
    dynamic_term = frequency * frequency
    alpha2 = (lame + 2. * shear_modulus) / density
    beta2 = shear_modulus / density
    gamma = 4. * pi * G_to_use * density / 3.

    # Optimizations
    r_inverse = 1. / radius
    r2 = radius * radius

    # Helper functions
    k2_quad_pos = (dynamic_term / beta2) + ((dynamic_term + 4. * gamma) / alpha2)
    k2_quad_neg = (dynamic_term / beta2) - ((dynamic_term + 4. * gamma) / alpha2)
    k2_quad = k2_quad_neg**2 + ((4. * order_l * (order_l + 1.) * gamma**2) / (alpha2 * beta2))

    # TODO: TS74 has these flipped compared to KMN15. Going with TS74 for this func.
    k2_pos = (1. / 2.) * (k2_quad_pos - sqrt_neg(k2_quad, is_real=True))
    k2_neg = (1. / 2.) * (k2_quad_pos + sqrt_neg(k2_quad, is_real=True))

    f_k2_pos = (beta2 * k2_pos - dynamic_term) / gamma
    f_k2_neg = (beta2 * k2_neg - dynamic_term) / gamma

    h_k2_pos = f_k2_pos - (order_l + 1.)
    h_k2_neg = f_k2_neg - (order_l + 1.)

    # Calculate Takeuchi and Saito functions
    # TODO: do we need to worry about the plus/minus on this sqrt?
    #    [GitHub Issue](https://github.com/jrenaud90/TidalPy/issues/31)
    # z_k2_pos = sqrt_neg(k2_pos, is_real=False) * radius
    # z_k2_neg = sqrt_neg(k2_neg, is_real=False) * radius
    z_k2_pos = k2_pos * r2
    z_k2_neg = k2_neg * r2

    phi_k2_pos, phi_lp1_k2_pos, psi_k2_pos = takeuchi_phi_psi(z_k2_pos, order_l)
    phi_k2_neg, phi_lp1_k2_neg, psi_k2_neg = takeuchi_phi_psi(z_k2_neg, order_l)

    # See Eq. 102 in TS72
    # # y1 solutions
    y1_s1 = ((-radius**(order_l + 1.)) / (2. * order_l + 3.)) * (0.5 * order_l * h_k2_pos * psi_k2_pos +
                                                                 f_k2_pos * phi_lp1_k2_pos)
    y1_s2 = ((-radius**(order_l + 1.)) / (2. * order_l + 3.)) * (0.5 * order_l * h_k2_neg * psi_k2_neg +
                                                                 f_k2_neg * phi_lp1_k2_neg)
    y1_s3 = order_l * radius**(order_l - 1.)

    # # y2 solutions
    y2_s1 = -(lame + 2. * shear_modulus) * radius**order_l * f_k2_pos * phi_k2_pos + \
            (shear_modulus * radius**order_l / (2. * order_l + 3.)) * \
            (-order_l * (order_l - 1.) * h_k2_pos * psi_k2_pos +
             2. * (2. * f_k2_pos + order_l * (order_l + 1.)) * phi_lp1_k2_pos)
    y2_s2 = -(lame + 2. * shear_modulus) * radius**order_l * f_k2_neg * phi_k2_neg + \
            (shear_modulus * radius**order_l / (2. * order_l + 3.)) * \
            (-order_l * (order_l - 1.) * h_k2_neg * psi_k2_neg +
             2. * (2. * f_k2_neg + order_l * (order_l + 1.)) * phi_lp1_k2_neg)
    y2_s3 = 2. * shear_modulus * order_l * (order_l - 1.) * radius**(order_l - 2.)

    # # y3 solutions
    y3_s1 = (-radius**(order_l + 1.) / (2. * order_l + 3.)) * (0.5 * h_k2_pos * psi_k2_pos - phi_lp1_k2_pos)
    y3_s2 = (-radius**(order_l + 1.) / (2. * order_l + 3.)) * (0.5 * h_k2_neg * psi_k2_neg - phi_lp1_k2_neg)
    y3_s3 = radius**(order_l - 1.)

    # # y4 solutions
    y4_s1 = shear_modulus * radius**order_l * \
            (phi_k2_pos - (1. / (2. * order_l + 3.)) * ((order_l - 1.) * h_k2_pos * psi_k2_pos +
                                                        2. * (f_k2_pos + 1.) * phi_lp1_k2_pos))
    y4_s2 = shear_modulus * radius**order_l * \
            (phi_k2_neg - (1. / (2. * order_l + 3.)) * ((order_l - 1.) * h_k2_neg * psi_k2_neg +
                                                        2. * (f_k2_neg + 1.) * phi_lp1_k2_neg))
    y4_s3 = 2. * shear_modulus * (order_l - 1.) * radius**(order_l - 2.)

    # # y5 solutions
    y5_s1 = radius**(order_l + 2.) * ((alpha2 * f_k2_pos - (order_l + 1.) * beta2) / r2 -
                                      (3. * gamma * f_k2_pos / (2. * (2. * order_l + 3.))) * psi_k2_pos)
    y5_s2 = radius**(order_l + 2.) * ((alpha2 * f_k2_neg - (order_l + 1.) * beta2) / r2 -
                                      (3. * gamma * f_k2_neg / (2. * (2. * order_l + 3.))) * psi_k2_neg)
    y5_s3 = (order_l * gamma - dynamic_term) * radius**order_l

    # # y6 solutions
    y6_s1 = (2. * order_l + 1.) * r_inverse * y5_s1 + \
            (3. * order_l * gamma * h_k2_pos * radius**(order_l + 1.) / (2. * (2. * order_l + 3.))) * psi_k2_pos
    y6_s2 = (2. * order_l + 1.) * r_inverse * y5_s2 + \
            (3. * order_l * gamma * h_k2_neg * radius**(order_l + 1.) / (2. * (2. * order_l + 3.))) * psi_k2_neg
    y6_s3 = (2. * order_l + 1.) * r_inverse * y5_s3 - \
            3. * order_l * gamma * radius**(order_l - 1.)

    # TODO: Right now numba does not support np.stack for purely scalar inputs. A temp fix is to make sure all the
    #    inputs are cast into arrays. See the github issue here: https://github.com/numba/numba/issues/7002
    y1_s1 = np.asarray(y1_s1, dtype=np.complex128)
    y2_s1 = np.asarray(y2_s1, dtype=np.complex128)
    y3_s1 = np.asarray(y3_s1, dtype=np.complex128)
    y4_s1 = np.asarray(y4_s1, dtype=np.complex128)
    y5_s1 = np.asarray(y5_s1, dtype=np.complex128)
    y6_s1 = np.asarray(y6_s1, dtype=np.complex128)

    y1_s2 = np.asarray(y1_s2, dtype=np.complex128)
    y2_s2 = np.asarray(y2_s2, dtype=np.complex128)
    y3_s2 = np.asarray(y3_s2, dtype=np.complex128)
    y4_s2 = np.asarray(y4_s2, dtype=np.complex128)
    y5_s2 = np.asarray(y5_s2, dtype=np.complex128)
    y6_s2 = np.asarray(y6_s2, dtype=np.complex128)

    y1_s3 = np.asarray(y1_s3, dtype=np.complex128)
    y2_s3 = np.asarray(y2_s3, dtype=np.complex128)
    y3_s3 = np.asarray(y3_s3, dtype=np.complex128)
    y4_s3 = np.asarray(y4_s3, dtype=np.complex128)
    y5_s3 = np.asarray(y5_s3, dtype=np.complex128)
    y6_s3 = np.asarray(y6_s3, dtype=np.complex128)

    # Combine the three solutions
    tidaly_s1 = np.stack((y1_s1, y2_s1, y3_s1, y4_s1, y5_s1, y6_s1))
    tidaly_s2 = np.stack((y1_s2, y2_s2, y3_s2, y4_s2, y5_s2, y6_s2))
    tidaly_s3 = np.stack((y1_s3, y2_s3, y3_s3, y4_s3, y5_s3, y6_s3))

    return nbList([tidaly_s1, tidaly_s2, tidaly_s3])


@njit(cacheable=True)
def liquid_guess_takeuchi(
    radius: FloatArray, bulk_modulus: NumArray,
    density: FloatArray, frequency: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> LiquidDynamicGuess:
    """ Calculate the initial guess at the bottom of a liquid layer using the dynamic assumption.

    This function uses the Takeuchi and Saito 1972 equations (Eq. 95-101).

    Using the dynamic assumption in a liquid layer results in two independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
       bulk dissipation (no shear dissipation within liquid layers).

    References
    ----------
    TS72

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    bulk_modulus : NumArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at  at `radius` [kg m-3]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : LiquidDynamicGuess
        The two independent liquid guesses (sn1, sn2)

    """

    # Convert compressibility parameters
    lame = bulk_modulus

    # Constants (See Eqs. B13-B16 of KMN15)
    dynamic_term = frequency * frequency
    alpha2 = lame / density
    gamma = 4. * pi * G_to_use * density / 3.

    # Optimizations
    r_inverse = 1. / radius
    r2 = radius * radius

    # Helper functions
    k2 = (1. / alpha2) * (dynamic_term + 4. * gamma - order_l * (order_l + 1.) * gamma**2 / dynamic_term)
    f = -dynamic_term / gamma
    h = f - (order_l + 1.)

    # Calculate Takeuchi and Saito functions
    z = k2 * r2
    phi, phi_lp1, psi = takeuchi_phi_psi(z, order_l)

    # See Eq. 102 in TS72
    # # y1 solutions
    y1_s1 = ((-radius**(order_l + 1.)) / (2. * order_l + 3.)) * (0.5 * order_l * h * psi + f * phi_lp1)
    y1_s2 = order_l * radius**(order_l - 1.)

    # # y2 solutions
    y2_s1 = -lame * radius**order_l * f * phi
    y2_s2 = np.zeros_like(radius, dtype=np.complex128)

    # # y5 solutions
    y5_s1 = radius**(order_l + 2.) * ((alpha2 * f / r2) - (3. * gamma * f / (2. * (2. * order_l + 3.))) * psi)
    y5_s2 = (order_l * gamma - dynamic_term) * radius**order_l

    # # y6 solutions
    y6_s1 = (2. * order_l + 1.) * r_inverse * y5_s1 + \
            (3. * order_l * gamma * h * radius**(order_l + 1.) / (2. * (2. * order_l + 3.))) * psi
    y6_s2 = (2. * order_l + 1.) * r_inverse * y5_s2 - \
            3. * order_l * gamma * radius**(order_l - 1.)

    # TODO: Right now numba does not support np.stack for purely scalar inputs. A temp fix is to make sure all the
    #    inputs are cast into arrays. See the github issue here: https://github.com/numba/numba/issues/7002
    y1_s1 = np.asarray(y1_s1, dtype=np.complex128)
    y2_s1 = np.asarray(y2_s1, dtype=np.complex128)
    y5_s1 = np.asarray(y5_s1, dtype=np.complex128)
    y6_s1 = np.asarray(y6_s1, dtype=np.complex128)

    y1_s2 = np.asarray(y1_s2, dtype=np.complex128)
    y2_s2 = np.asarray(y2_s2, dtype=np.complex128)
    y5_s2 = np.asarray(y5_s2, dtype=np.complex128)
    y6_s2 = np.asarray(y6_s2, dtype=np.complex128)

    # Combine the three solutions
    tidaly_s1 = np.stack((y1_s1, y2_s1, y5_s1, y6_s1))
    tidaly_s2 = np.stack((y1_s2, y2_s2, y5_s2, y6_s2))

    return nbList([tidaly_s1, tidaly_s2])


# TODO: Development paused on the Martens power series method
"""
# @njit(cacheable=True)
def solid_guess_martens(
    radius: FloatArray, shear_modulus: NumArray, bulk_modulus: NumArray,
    density: FloatArray, frequency: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> SolidDynamicGuess:

    # TODO: This looks like a more stable solution for the inner core, but I can't figure out how to implement it right now.

    # Optimizations and shortcuts
    w2 = frequency * frequency
    grav_term = 4. * pi * G_to_use * density

    # Other functions
    k2 = order_l * (order_l + 1.)
    lame = (bulk_modulus - (2. / 3.) * shear_modulus)
    gamma = grav_term / 3.
    zeta = 1. / (lame + 2. * shear_modulus)
    delta = 2. * shear_modulus * (3. * lame + 2. * shear_modulus) * zeta
    epsilon = 4. * k2 * shear_modulus * (lame + shear_modulus) * zeta - 2. * shear_modulus

    # # See Martens Eqs. 4.98 and 4.99
    p1 = 2. * order_l * (order_l * (order_l + 2.) * lame + (order_l * (order_l + 2) - 1.) * shear_modulus)
    q1 = (order_l * (order_l + 1.) + order_l * (order_l + 3.)) * lame + 2. * order_l * (order_l + 1.) * shear_modulus
    p2 = order_l * (order_l + 5) + (lame * order_l * (order_l + 3.)) / shear_modulus
    q2 = 2. * (order_l + 1.) + (lame * (order_l + 3.) / shear_modulus)

    # # First solution
    # Set the values of the arbitrary coefficients
    A_11 = 1.
    A_61 = 1.
    A_40 = 1.

    # Solve for the dependent coefficients (Martens Eq. 4.97)
    A_33 = A_11 * density * ((3. - order_l) * gamma + w2) / p1
    A_13 = A_33 * -order_l
    A_22 = A_33 * -q1
    A_42 = 0.
    A_54 = grav_term * ((order_l + 3.) * A_13 - k2 * A_33) / (2. * (2. * order_l + 3.))
    A_63 = (order_l + 2.) * A_54 - grav_term * A_13

    # Eq. 4.102
    B_33 = A_61 * density / p1
    B_13 = B_33 * -order_l
    B_22 = B_33 * -q1
    B_42 = 0.
    B_54 = (grav_term / (2. * (2. * order_l + 3.))) * ((order_l + 3.) * B_13 - k2 * B_33)
    B_63 = (order_l + 2.) * B_54 - grav_term * B_13

    # Eq. 4.104
    p_ratio = p2 / p1
    C_11 = A_40 * (1. / shear_modulus - order_l * p_ratio)
    C_20 = A_40 * (q2 - q1 * p_ratio)
    C_31 = A_40 * p_ratio
    C_52 = (grav_term / (2. * (2. * order_l + 3.))) * ((order_l + 3.) * C_11 - k2 * C_31)
    C_61 = (order_l + 2.) * C_52 - grav_term * C_11

    # Solve for other coefficients which requires matrix operations
    A_matrix_1 = np.asarray(
        [
            [2. * lame * zeta + order_l + 3., -zeta, -k2 * lame * zeta, 0., 0., 0.],
            [-2. * delta, 4 * shear_modulus * zeta + order_l + 2., k2 * delta, -k2, 0., 0.],
            [1., 0., order_l + 2., -(1. / shear_modulus), 0., 0.],
            [delta, lame * zeta, -epsilon, order_l + 5., 0., 0.],
            [-3. * gamma, 0., 0., 0., order_l + 4., -1.],
            [0., 0., 3. * gamma * k2, 0., -k2, order_l + 5.],
            ]
        )
    A_matrix_2 = A_matrix_1 + 2. * np.eye(6, dtype=A_matrix_1.dtype)
    A_matrix_1_inv = np.linalg.inv(A_matrix_1)
    A_matrix_2_inv = np.linalg.inv(A_matrix_2)

    # Solve for first set of coeffs
    A_vector = np.asarray(
        [
            [0.],
            [density * (-(4. * gamma + w2) * A_13 + (k2 * gamma) * A_33 - A_63)],
            [0.],
            [density * (gamma * A_13 - w2 * A_33 - A_54)],
            [0.],
            [0.]
            ]
        )
    B_vector = np.asarray(
        [
            [0.],
            [density * (-(4. * gamma + w2) * A_13 + (k2 * gamma) * A_33 - A_63)],
            [0.],
            [density * (gamma * A_13 - w2 * A_33 - A_54)],
            [0.],
            [0.]
            ]
        )
    C_vector_1 = np.asarray(
        [
            [0.],
            [density * (-(4. * gamma + w2) * C_11 + (k2 * gamma) * C_31 - C_61)],
            [0.],
            [density * (gamma * C_11 - w2 * C_31 - C_52)],
            [0.],
            [0.]
            ]
        )
    C_vector_2 = np.asarray(
        [
            [0.],
            [density * (-(4. * gamma + w2) * C_13 + (k2 * gamma) * C_33 - C_63)],
            [0.],
            [density * (gamma * C_13 - w2 * C_33 - C_54)],
            [0.],
            [0.]
            ]
        )
    A_15, A_24, A_35, A_44, A_56, A_65 = (A_matrix_1_inv @ LHS_vector_1)[:, 0]


    A_Matrix_inv = np.linalg.inv(A_Matrix)

    A_1_coeff = A_Matrix_inv @ LHS_vector

    return A_1_coeff

"""