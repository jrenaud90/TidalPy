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

from typing import Tuple, Union

from ....constants import pi, G
from ....utilities.performance import njit
from ....utilities.types import FloatArray, ComplexArray

CmplxFltArray = Union[FloatArray, ComplexArray]
RadialFuncSolidStaticType = Tuple[CmplxFltArray, CmplxFltArray, CmplxFltArray, CmplxFltArray, CmplxFltArray, ComplexArray]
RadialFuncLiquidStaticType = Tuple[CmplxFltArray, CmplxFltArray]


@njit(cacheable=True)
def radial_derivatives_solid_general(radius : FloatArray, radial_functions : RadialFuncSolidStaticType,
                                     shear_modulus : CmplxFltArray, bulk_modulus : CmplxFltArray, density : FloatArray,
                                     gravity : FloatArray, order_l: int = 2) -> RadialFuncSolidStaticType:
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
        Radius where the radial functions are calculated. [m]
    radial_functions : RadialFuncSolidStaticType
        Tuple of radial functions for a solid layer (y1, y2, y3, y4, y5, y6)
    shear_modulus : CmplxFltArray
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : CmplxFltArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at  at `radius` [kg m-3]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2]
    order_l : int = 2
        Tidal harmonic

    Returns
    -------
    radial_derivatives : RadialFuncSolidStaticType
        The radial derivatives of the radial functions

    """

    y1, y2, y3, y4, y5, y6 = radial_functions

    # Convert compressibility parameters (the first lame parameter can be complex)
    lame = (bulk_modulus - (2. / 3.) * shear_modulus)

    # Optimizations
    lame_2mu = 1. / (lame + 2. * shear_modulus)
    r_inverse = 1. / radius
    r2_inverse = r_inverse * r_inverse
    llp1 = order_l * (order_l + 1.)
    grav_term = 4. * pi * G * density

    # See Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    # # dy2 and dy4 contain: viscoelastic, and gravitational terms.
    dy1 = \
        y1 * -2. * lame * lame_2mu * r_inverse + \
        y2 * lame_2mu + \
        y3 * llp1 * lame * lame_2mu * r_inverse

    dy2 = \
        y1 * (12. * bulk_modulus * shear_modulus * lame_2mu * r2_inverse - 4. * density * gravity * r_inverse) - \
        y2 * (4. * shear_modulus * lame_2mu * r_inverse) + \
        y3 * r_inverse * llp1 * (density * gravity - 6. * bulk_modulus * shear_modulus * lame_2mu * r_inverse) + \
        y4 * llp1 * r_inverse + \
        y5 * (order_l + 1.) * r_inverse * density - \
        y6 * density

    dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / shear_modulus)

    dy4 = \
        y1 * (density * gravity * r_inverse - 6. * bulk_modulus * shear_modulus * lame_2mu * r2_inverse) - \
        y2 * lame * lame_2mu * r_inverse + \
        y3 * (4. * llp1 * (lame + shear_modulus) * shear_modulus * lame_2mu * r2_inverse -
              2. * shear_modulus * r2_inverse) - \
        y4 * 3. * r_inverse - \
        y5 * density * r_inverse

    dy5 = \
        y1 * grav_term - \
        y5 * (order_l + 1.) * r_inverse + \
        y6

    dy6 = \
        y1 * grav_term * (order_l + 1.) * r_inverse - \
        y3 * grav_term * llp1 * r_inverse + \
        y6 * (order_l - 1.) * r_inverse

    return dy1, dy2, dy3, dy4, dy5, dy6

@njit(cacheable=True)
def radial_derivatives_liquid_general(radius: FloatArray, radial_functions: RadialFuncLiquidStaticType,
                                      density: FloatArray, gravity: FloatArray,
                                      order_l: int = 2) -> RadialFuncLiquidStaticType:
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
    radial_functions : RadialFuncLiquidStaticType
        Tuple of radial functions for a solid layer (y5, y7)
    density : FloatArray
        Density at  at `radius` [kg m-3]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2]
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    radial_derivatives : RadialFuncLiquidStaticType
        The radial derivatives of the radial functions

    """

    # For the dynamic version, y4 = 0 always in a liquid layer and y1 and y2 are not defined uniquely
    y5, y7 = radial_functions

    # Optimizations
    r_inverse = 1. / radius
    grav_term = 4. * pi * G * density

    # See Eq. 18 in S75
    dy5 = \
        y5 * (grav_term / gravity - (order_l + 1.) * r_inverse) + \
        y7

    dy7 = \
        y5 * 2. * (order_l - 1.) * r_inverse * grav_term / gravity + \
        y7 * ((order_l - 1.) * r_inverse - grav_term / gravity)

    return dy5, dy7