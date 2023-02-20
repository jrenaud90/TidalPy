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

from typing import Union

import numpy as np

from .functions import takeuchi_phi_psi, z_calc
from TidalPy.constants import G, pi
from TidalPy.utilities.math.special import sqrt_neg
from TidalPy.utilities.performance import njit, nbList
from TidalPy.utilities.types import ComplexArray

SolidDynamicGuess = nbList[ComplexArray]
LiquidDynamicGuess = nbList[ComplexArray]


@njit(cacheable=True)
def solid_guess_kamata(
    radius: float, shear_modulus: Union[float, complex], bulk_modulus: Union[float, complex],
    density: float, frequency: float,
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
    radius : float
        Radius where the radial functions are calculated. [m]
    shear_modulus : Union[float, complex]
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : Union[float, complex]
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : float
        Density at `radius` [kg m-3]
    frequency : float
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

    tidaly_s1 = np.empty(6, dtype=np.complex128)
    tidaly_s1[0] = y1_s1
    tidaly_s1[1] = y2_s1
    tidaly_s1[2] = y3_s1
    tidaly_s1[3] = y4_s1
    tidaly_s1[4] = y5_s1
    tidaly_s1[5] = y6_s1

    tidaly_s2 = np.empty(6, dtype=np.complex128)
    tidaly_s2[0] = y1_s2
    tidaly_s2[1] = y2_s2
    tidaly_s2[2] = y3_s2
    tidaly_s2[3] = y4_s2
    tidaly_s2[4] = y5_s2
    tidaly_s2[5] = y6_s2

    tidaly_s3 = np.empty(6, dtype=np.complex128)
    tidaly_s3[0] = y1_s3
    tidaly_s3[1] = y2_s3
    tidaly_s3[2] = y3_s3
    tidaly_s3[3] = y4_s3
    tidaly_s3[4] = y5_s3
    tidaly_s3[5] = y6_s3

    return nbList([tidaly_s1, tidaly_s2, tidaly_s3])


@njit(cacheable=True)
def liquid_guess_kamata(
    radius: float, bulk_modulus: Union[float, complex],
    density: float, frequency: float,
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
    radius : float
        Radius where the radial functions are calculated. [m]
    bulk_modulus : Union[float, complex]
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : float
        Density at `radius` [kg m-3]
    frequency : float
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

    tidaly_s1 = np.empty(4, dtype=np.complex128)
    tidaly_s1[0] = y1_s1
    tidaly_s1[1] = y2_s1
    tidaly_s1[2] = y5_s1
    tidaly_s1[3] = y6_s1

    tidaly_s2 = np.empty(4, dtype=np.complex128)
    tidaly_s2[0] = y1_s2
    tidaly_s2[1] = y2_s2
    tidaly_s2[2] = y5_s2
    tidaly_s2[3] = y6_s2

    return nbList([tidaly_s1, tidaly_s2])


@njit(cacheable=True)
def solid_guess_takeuchi(
    radius: float, shear_modulus: Union[float, complex], bulk_modulus: Union[float, complex],
    density: float, frequency: float,
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
    radius : float
        Radius where the radial functions are calculated. [m]
    shear_modulus : Union[float, complex]
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : Union[float, complex]
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : float
        Density at `radius` [kg m-3]
    frequency : float
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

    tidaly_s1 = np.empty(6, dtype=np.complex128)
    tidaly_s1[0] = y1_s1
    tidaly_s1[1] = y2_s1
    tidaly_s1[2] = y3_s1
    tidaly_s1[3] = y4_s1
    tidaly_s1[4] = y5_s1
    tidaly_s1[5] = y6_s1

    tidaly_s2 = np.empty(6, dtype=np.complex128)
    tidaly_s2[0] = y1_s2
    tidaly_s2[1] = y2_s2
    tidaly_s2[2] = y3_s2
    tidaly_s2[3] = y4_s2
    tidaly_s2[4] = y5_s2
    tidaly_s2[5] = y6_s2

    tidaly_s3 = np.empty(6, dtype=np.complex128)
    tidaly_s3[0] = y1_s3
    tidaly_s3[1] = y2_s3
    tidaly_s3[2] = y3_s3
    tidaly_s3[3] = y4_s3
    tidaly_s3[4] = y5_s3
    tidaly_s3[5] = y6_s3

    return nbList([tidaly_s1, tidaly_s2, tidaly_s3])


@njit(cacheable=True)
def liquid_guess_takeuchi(
    radius: float, bulk_modulus: Union[float, complex],
    density: float, frequency: float,
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
    radius : float
        Radius where the radial functions are calculated. [m]
    bulk_modulus : Union[float, complex]
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : float
        Density at `radius` [kg m-3]
    frequency : float
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

    tidaly_s1 = np.empty(4, dtype=np.complex128)
    tidaly_s1[0] = y1_s1
    tidaly_s1[1] = y2_s1
    tidaly_s1[2] = y5_s1
    tidaly_s1[3] = y6_s1

    tidaly_s2 = np.empty(4, dtype=np.complex128)
    tidaly_s2[0] = y1_s2
    tidaly_s2[1] = y2_s2
    tidaly_s2[2] = y5_s2
    tidaly_s2[3] = y6_s2

    return nbList([tidaly_s1, tidaly_s2])
#
# @njit(cacheable=True)
# def solid_guess_power(
#     radius: FloatArray, shear_modulus: NumArray, bulk_modulus: NumArray,
#     density: FloatArray, frequency: FloatArray,
#     order_l: int = 2, G_to_use: float = G
#     ) -> SolidDynamicGuess:
#     """ Calculate the initial guess at the bottom of a solid layer using the dynamic assumption.
#
#     This function uses the power series described in MartensThesis
#
#     Using the dynamic assumption in a solid layer results in three independent solutions for the radial derivatives.
#
#     These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
#        bulk and shear dissipation.
#
#     References
#     ----------
#     MartensThesis
#
#     Parameters
#     ----------
#     radius : FloatArray
#         Radius where the radial functions are calculated. [m]
#     shear_modulus : NumArray
#         Shear modulus (can be complex for dissipation) at `radius` [Pa]
#     bulk_modulus : NumArray
#         Bulk modulus (can be complex for dissipation) at `radius` [Pa]
#     density : FloatArray
#         Density at `radius` [kg m-3]
#     frequency : FloatArray
#         Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
#     order_l : int = 2
#         Tidal harmonic order.
#     G_to_use : float = G
#         Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
#
#     Returns
#     -------
#     solid_guesses : SolidDynamicGuess
#         The three independent solid guesses (sn1, sn2, sn3)
#
#     """
#
#     # Check for model problems
#     if order_l < 2:
#         raise Exception('Model not designed for l=0,1. See Smylie 2013.')
#
#     # Convert compressibility parameters
#     lame = bulk_modulus - (2. / 3.) * shear_modulus
#
#     # Constants (See Eqs. B13-B16 of KMN15)
#     dynamic_term = frequency * frequency
#     alpha2 = (lame + 2. * shear_modulus) / density
#     beta2 = shear_modulus / density
#     gamma = 4. * pi * G_to_use * density / 3.
#
#     # Optimizations
#     r_inverse = 1. / radius
#     r2 = radius * radius
#     r4 = r2 * r2
#
#     # Other functions
#     k2 = order_l * (order_l + 1.)
#     lame = (bulk_modulus - (2. / 3.) * shear_modulus)
#     grav_term = G_to_use * 4. * pi * density
#     gamma = grav_term / 3.
#     zeta = 1. / (lame + 2. * shear_modulus)
#     delta = 2. * shear_modulus * (3. * lame + 2. * shear_modulus) * zeta
#     epsilon = 4. * k2 * shear_modulus * (lame + shear_modulus) * zeta - 2. * shear_modulus
#
#     # # See Martens Eqs. 4.98 and 4.99
#     p1 = 2. * order_l * (order_l * (order_l + 2.) * lame + (order_l * (order_l + 2) - 1.) * shear_modulus)
#     q1 = (order_l * (order_l + 1.) + order_l * (order_l + 3.)) * lame + 2. * order_l * (order_l + 1.) * shear_modulus
#     p2 = order_l * (order_l + 5) + (lame * order_l * (order_l + 3.)) / shear_modulus
#     q2 = 2. * (order_l + 1.) + (lame * (order_l + 3.) / shear_modulus)
#
#     # # First solution
#     # Set the values of the arbitrary coefficients
#     A_11 = 1.
#     B_61 = 1.
#     C_40 = 1.
#
#     # Solve for the dependent coefficients (Martens Eq. 4.97)
#     A_33 = A_11 * density * ((3. - order_l) * gamma + dynamic_term) / p1
#     A_13 = A_33 * -order_l
#     A_22 = A_33 * -q1
#     A_42 = 0.
#     A_54 = grav_term * ((order_l + 3.) * A_13 - k2 * A_33) / (2. * (2. * order_l + 3.))
#     A_63 = (order_l + 2.) * A_54 - grav_term * A_13
#
#     # Eq. 4.102
#     B_33 = B_61 * density / p1
#     B_13 = B_33 * -order_l
#     B_22 = B_33 * -q1
#     B_42 = 0.
#     B_54 = (grav_term / (2. * (2. * order_l + 3.))) * ((order_l + 3.) * B_13 - k2 * B_33)
#     B_63 = (order_l + 2.) * B_54 - grav_term * B_13
#
#     # Eq. 4.104
#     p_ratio = p2 / p1
#     C_11 = C_40 * (1. / shear_modulus - order_l * p_ratio)
#     C_20 = C_40 * (q2 - q1 * p_ratio)
#     C_31 = C_40 * p_ratio
#     C_52 = (grav_term / (2. * (2. * order_l + 3.))) * ((order_l + 3.) * C_11 - k2 * C_31)
#     C_61 = (order_l + 2.) * C_52 - grav_term * C_11
#
#     # Solve for other coefficients which requires matrix operations
#     # Matrix 1 defined in Eq. 4.100
#     matrix_1 = np.asarray(
#         [
#             [2. * lame * zeta + order_l + 3., -zeta, -k2 * lame * zeta, 0., 0., 0.        ],
#             [-2. * delta, 4 * shear_modulus * zeta + order_l + 2., k2 * delta, -k2, 0., 0.],
#             [1., 0., order_l + 2., -(1. / shear_modulus), 0., 0.                          ],
#             [delta, lame * zeta, -epsilon, order_l + 5., 0., 0.],
#             [-3. * gamma, 0., 0., 0., order_l + 4., -1.],
#             [0., 0., 3. * gamma * k2, 0., -k2, order_l + 5.],
#             ]
#         )
#     matrix_2 = matrix_1 + 2. * np.eye(6, dtype=matrix_1.dtype)
#
#     # Perform matrix inversions
#     matrix_1_inv = np.linalg.inv(matrix_1)
#     matrix_2_inv = np.linalg.inv(matrix_2)
#
#     # Solve for Solution 1, r^4 coefficients. Eq. 4.100
#     A4_vector = np.asarray(
#         [
#             [0.],
#             [density * (-(4. * gamma + dynamic_term) * A_13 + (k2 * gamma) * A_33 - A_63)],
#             [0.],
#             [density * (gamma * A_13 - dynamic_term * A_33 - A_54)],
#             [0.],
#             [0.]
#             ]
#         )
#     A_15, A_24, A_35, A_44, A_56, A_65 = (matrix_1_inv @ A4_vector)[:, 0]
#
#     # Solve for Solution 2, r^4 coefficients - The LHS has the same form as solution 1s. Eq. 4.100
#     B4_vector = np.asarray(
#         [
#             [0.],
#             [density * (-(4. * gamma + dynamic_term) * B_13 + (k2 * gamma) * B_33 - B_63)],
#             [0.],
#             [density * (gamma * B_13 - dynamic_term * B_33 - B_54)],
#             [0.],
#             [0.]
#         ]
#     )
#     B_15, B_24, B_35, B_44, B_56, B_65 = (matrix_1_inv @ B4_vector)[:, 0]
#
#     # Solve for Solution 3, r^2 coefficients. Eq. 4.107
#     C2_vector = np.asarray(
#         [
#             [0.],
#             [density * (-(4. * gamma + dynamic_term) * C_11 + (k2 * gamma) * C_31 - C_61)],
#             [0.],
#             [density * (gamma * C_11 - dynamic_term * C_31 - C_52)],
#             [0.],
#             [0.]
#         ]
#     )
#     C_13, C_22, C_33, C_42, C_54, C_63 = (matrix_1_inv @ C2_vector)[:, 0]
#
#     # Solve for Solution 3, r^4 coefficients. Eq. 4.108
#     C4_vector = np.asarray(
#         [
#             [0.],
#             [density * (-(4. * gamma + dynamic_term) * C_13 + (k2 * gamma) * C_33 - C_63)],
#             [0.],
#             [density * (gamma * C_13 - dynamic_term * C_33 - C_54)],
#             [0.],
#             [0.]
#         ]
#     )
#     C_15, C_24, C_35, C_44, C_56, C_65 = (matrix_2_inv @ C4_vector)[:, 0]
#
#     # Find z functions based on the coefficients
#
#     # Solution 1
#     z1_sol1 = A_11 + A_13 * r2 + A_15 * r4
#     z2_sol1 = 2. * (order_l - 1.) * shear_modulus * A_11 + A_22 * r2 + A_24 * r4
#     z3_sol1 = (1. / order_l) * A_11 + A_33 * r2 + A_35 * r4
#     # LEFTOFF
#
#     # Solution 2
#     z1_sol2 = -order_l * (density / p1) * B_61 * r2 + B_15 * r4
#     z2_sol2 = -(q1 / p1) * density * B_61 * r2 + B_24 * r4
#     z3_sol2 = (density / p1) * B_61 * r2 + B_35 * r4
#     z4_sol2 = B_44 * r4
#     z5_sol2 = (1. / order_l) * B_61 + B_54 * r2 + B_56 * r4
#     z6_sol2 = B_61 + B_63 * r2 + B_65 * r4
#
#     # Solution 3
#     z1_sol3 = ((1. / shear_modulus) - order_l * p_ratio) * C_40 + C_13 * r2 + C_15 * r4
#     z2_sol3 = (q2 - q1 * p_ratio) * C_40 + C_22 * r2 + C_24 * r4
#     z3_sol3 = p_ratio * C_40 + C_33 * r2 + C_35 * r4
#     z4_sol3 = C_40 + C_42 * r2 + C_44 * r4
#     z5_sol3 = C_52 + C_54 * r2 + C_56 * r4
#     z6_sol3 = C_61 + C_63 * r2 + C_65 * r4
#
#     return A_1_coeff
#
# """