""" Functions to calculate the Rayleigh wave propagation differential equations

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

from typing import Tuple

from .....constants import G, pi
from .....utilities.performance import njit
from .....utilities.types import FloatArray, NumArray


@njit(cacheable=True)
def radial_derivatives_solid_general(
    radius: FloatArray,
    radial_functions: Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray],
    shear_modulus: NumArray, bulk_modulus: NumArray, density: FloatArray,
    gravity: FloatArray, frequency: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]:
    """ Calculates the radial derivative of the radial functions in the most general form - for solid layers.

    Allows for compressibility and dynamic tides.
    Tidal harmonic l is allowed to be an integer >= 2.

    OPT: This and the other radial derivative functions are called many times during integration.
        As of TidalPy V0.3.4 this function call for floats is around 2.3 micro seconds.

    References
    ----------
    KMN15; B15; TS72

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    radial_functions : Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]
        Tuple of radial functions for a solid layer (y1, y2, y3, y4, y5, y6)
    shear_modulus : NumArray
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : NumArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at `radius` [kg m-3]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Newton's gravitational constant. This can be provided as in its MKS units or dimensionless to match the other
        inputs.

    Returns
    -------
    radial_derivatives : Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]
        The derivatives of the radial functions for a solid layer (dynamic assumption)

    """

    y1, y2, y3, y4, y5, y6 = radial_functions

    # Convert compressibility parameters (the first lame parameter can be complex)
    lame = (bulk_modulus - (2. / 3.) * shear_modulus)

    # Optimizations
    lp1 = order_l + 1.
    lm1 = order_l - 1.
    llp1 = order_l * lp1
    lame_2mu = lame + 2. * shear_modulus
    lame_2mu_inverse = 1. / lame_2mu
    r_inverse = 1. / radius
    two_shear_r_inv = 2. * shear_modulus * r_inverse
    density_gravity = density * gravity
    dynamic_term = -frequency * frequency * density * radius
    grav_term = 4. * pi * G_to_use * density
    y1_y3_term = 2. * y1 - llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    dy1 = lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
    )

    dy2 = r_inverse * (
            y1 * (dynamic_term - 2. * density_gravity) +
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
            y3 * (dynamic_term - two_shear_r_inv) +
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

    return dy1, dy2, dy3, dy4, dy5, dy6


@njit(cacheable=True)
def radial_derivatives_liquid_general(
    radius: FloatArray,
    radial_functions: Tuple[NumArray, NumArray, NumArray, NumArray],
    bulk_modulus: NumArray, density: FloatArray,
    gravity: FloatArray, frequency: FloatArray,
    order_l: int = 2, G_to_use: float = G
    ) -> Tuple[NumArray, NumArray, NumArray, NumArray]:
    """ Calculates the radial derivative of the radial functions in the most general form - for liquid layers.

    Allows for compressibility and dynamic tides.
    Tidal harmonic l is allowed to be an integer >= 2.

    References
    ----------
    KMN15; B15; TS72

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    radial_functions : Tuple[ComplexArray, ComplexArray, ComplexArray, ComplexArray]
        Tuple of radial functions for a solid layer (y1, y2, y5, y6)
    bulk_modulus : NumArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at `radius` [kg m-3]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Newton's gravitational constant. This can be provided as in its MKS units or dimensionless to match the other
        inputs.

    Returns
    -------
    radial_derivatives : Tuple[NumArray, NumArray, NumArray, NumArray]
        The derivatives of the radial functions for a liquid layer (dynamic assumption)

    """

    # For the dynamic version, y4 = 0 always in a liquid layer and y3 is defined by y1, y2, and y5 analytically
    y1, y2, y5, y6 = radial_functions

    # Convert compressibility parameters (the first lame parameter can be complex)
    # For the liquid layer it is assumed that the shear modulus is zero so the lame parameter simply
    #    equals the bulk modulus
    lame = bulk_modulus

    # Optimizations
    lp1 = order_l + 1.
    lm1 = order_l - 1.
    llp1 = order_l * lp1
    lame_inverse = 1. / lame
    r_inverse = 1. / radius
    density_gravity = density * gravity
    dynamic_term = -frequency * frequency * density * radius
    grav_term = 4. * pi * G_to_use * density

    # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
    #   other derivatives.
    y3 = (1. / dynamic_term) * (y2 - density_gravity * y1 + density * y5)
    y1_y3_term = 2. * y1 - llp1 * y3

    # Eqs. 11--14 in KMN15 equations look like they don't match TS72 because they applied the rheology already.
    #    and substituted y3.
    # We will use TS72 eq. 87 to allow for a generic rheology and bulk dissipation.
    # # dy2 contain all three of: dynamic, viscoelastic, and gravitational terms.
    dy1 = \
        y2 * lame_inverse + \
        y1_y3_term * -r_inverse

    dy2 = r_inverse * (
            y1 * (dynamic_term - 2. * density_gravity) +
            y5 * density * lp1 +
            y6 * -density * radius +
            # TODO: In the solid version there is a [2. * (lame + shear_modulus) * r_inverse] coefficient for y1_y3_term
            #   In TS72 the first term is gone. Shouldn't Lame + mu = Lame = Bulk for liquid layer?
            y1_y3_term * -density_gravity
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

    return dy1, dy2, dy5, dy6
