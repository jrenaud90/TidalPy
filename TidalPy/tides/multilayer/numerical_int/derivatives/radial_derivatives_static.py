""" Functions to calculate the Rayleigh wave propagation differential equations

These functions do not allow for dynamic tides (w=0):
    - Compressibility
    - Liquid layer propagation

References
----------
KMN15 : Kamata+ (2015; JGR-P; DOI: 10.1002/2015JE004821)
B15   : Beuthe+ (2015; Icarus; DOI: 10.1016/j.icarus.2015.06.008)
TMS05 : Tobie+ (2005; Icarus; DOI: 10.1016/j.icarus.2005.04.006)
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

import numpy as np

from .....constants import G, pi
from .....utilities.performance import njit
from .....utilities.types import ComplexArray, FloatArray, NumArray


@njit(cacheable=False)
def radial_derivatives_solid_general(
    radius: FloatArray,
    radial_functions: np.ndarray,
    shear_modulus: NumArray, bulk_modulus: NumArray, density: FloatArray,
    gravity: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> np.ndarray:
    """ Calculates the derivatives of the radial functions using the static assumption - for solid layers.

    Allows for compressibility.
    Dynamic tides are not permitted (w=0)
    Tidal harmonic l is allowed to be >= 2.

    References
    ----------
    KMN15; B15; TS72

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m; or dimensionless]
    radial_functions : np.ndarray
        Tuple of radial functions for a solid layer (y1, y2, y3, y4, y5, y6)
    shear_modulus : NumArray
        Shear modulus (can be complex for dissipation) at `radius` [Pa; or dimensionless]
    bulk_modulus : NumArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa; or dimensionless]
    density : FloatArray
        Density at `radius` [kg m-3; or dimensionless]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2; or dimensionless]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Newton's gravitational constant. This can be provided as in its MKS units or dimensionless to match the other
        inputs.

    Returns
    -------
    radial_derivatives : np.ndarray
        The radial derivatives of the radial functions

    """

    y1_real = radial_functions[0]
    y1_imag = radial_functions[1]
    y2_real = radial_functions[2]
    y2_imag = radial_functions[3]
    y3_real = radial_functions[4]
    y3_imag = radial_functions[5]
    y4_real = radial_functions[6]
    y4_imag = radial_functions[7]
    y5_real = radial_functions[8]
    y5_imag = radial_functions[9]
    y6_real = radial_functions[10]
    y6_imag = radial_functions[11]

    # Convert floats to complex
    y1 = y1_real + 1.0j * y1_imag
    y2 = y2_real + 1.0j * y2_imag
    y3 = y3_real + 1.0j * y3_imag
    y4 = y4_real + 1.0j * y4_imag
    y5 = y5_real + 1.0j * y5_imag
    y6 = y6_real + 1.0j * y6_imag

    # Convert compressibility parameters (the first lame parameter can be complex)
    lame = (bulk_modulus - (2. / 3.) * shear_modulus)

    # Optimizations
    lp1              = order_l + 1.
    lm1              = order_l - 1.
    llp1             = order_l * lp1
    lame_2mu         = lame + 2. * shear_modulus
    lame_2mu_inverse = 1. / lame_2mu
    r_inverse        = 1. / radius
    two_shear_r_inv  = 2. * shear_modulus * r_inverse
    density_gravity  = density * gravity
    grav_term        = 4. * pi * G_to_use * density
    y1_y3_term       = 2. * y1 - llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # The static case just sets all frequency dependence in these equations to zero.
    # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    dy1 = lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
    )

    dy2 = r_inverse * (
            y1 * -2. * density_gravity +
            y2 * -2. +
            y4 * llp1 +
            y5 * density * lp1 +
            y6 * -density * radius +
            dy1 * 2. * lame +
            y1_y3_term * (2. * (lame + shear_modulus) * r_inverse - density_gravity)
    )

    dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / shear_modulus)

    dy4 = r_inverse * (
            y1 * (density_gravity + two_shear_r_inv) +
            y3 * -two_shear_r_inv +
            y4 * -3. +
            y5 * -density +
            dy1 * -lame +
            y1_y3_term * -lame_2mu * r_inverse
    )

    dy5 = \
        y1 * grav_term + \
        y5 * -lp1 * r_inverse + \
        y6

    dy6 = r_inverse * (
            y1 * grav_term * lm1 +
            y6 * lm1 +
            y1_y3_term * grav_term
    )

    # Convert back to floats
    dy1_real = np.real(dy1)
    dy1_imag = np.imag(dy1)
    dy2_real = np.real(dy2)
    dy2_imag = np.imag(dy2)
    dy3_real = np.real(dy3)
    dy3_imag = np.imag(dy3)
    dy4_real = np.real(dy4)
    dy4_imag = np.imag(dy4)
    dy5_real = np.real(dy5)
    dy5_imag = np.imag(dy5)
    dy6_real = np.real(dy6)
    dy6_imag = np.imag(dy6)

    dy = np.asarray(
            (dy1_real, dy1_imag, dy2_real, dy2_imag, dy3_real, dy3_imag,
             dy4_real, dy4_imag, dy5_real, dy5_imag, dy6_real, dy6_imag),
            dtype=np.float64)

    return dy


@njit(cacheable=False)
def radial_derivatives_liquid_general(
    radius: FloatArray, radial_functions: np.ndarray,
    density: FloatArray, gravity: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> np.ndarray:
    """ Calculates the derivatives of the radial functions using the static assumption - for liquid layers (mu = 0).

    Allows for compressibility (technically, but bulk mod is not used in this equations).
    Dynamic tides are not permitted (w=0)
    Tidal harmonic l is allowed to be >= 2.

    References
    ----------
    KMN15; B15; TS72; S74

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    radial_functions : np.ndarray
        Tuple of radial functions for a solid layer (y5, y7)
    density : FloatArray
        Density at `radius` [kg m-3]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Newton's gravitational constant. This can be provided as in its MKS units or dimensionless to match the other
        inputs.

    Returns
    -------
    radial_derivatives : np.ndarray
        The radial derivatives of the radial functions

    """

    # For the static version, only y5 and y7 are defined.
    y5_real = radial_functions[0]
    y5_imag = radial_functions[1]
    y7_real = radial_functions[2]
    y7_imag = radial_functions[3]

    # Convert floats to complex
    y5 = y5_real + 1.0j * y5_imag
    y7 = y7_real + 1.0j * y7_imag

    # Optimizations
    r_inverse = 1. / radius
    grav_term = 4. * pi * G_to_use * density / gravity

    # See Eq. 18 in S75
    dy5 = \
        y5 * (grav_term - (order_l + 1.) * r_inverse) + \
        y7

    dy7 = \
        y5 * 2. * (order_l - 1.) * r_inverse * grav_term + \
        y7 * ((order_l - 1.) * r_inverse - grav_term)

    # Convert back to floats
    dy5_real = np.real(dy5)
    dy5_imag = np.imag(dy5)
    dy7_real = np.real(dy7)
    dy7_imag = np.imag(dy7)

    dy = np.asarray(
            (dy5_real, dy5_imag, dy7_real, dy7_imag),
            dtype=np.float64)

    return dy
