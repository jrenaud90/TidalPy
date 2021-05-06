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

from typing import Tuple, Union

from ....constants import pi, G
from ....utilities.performance import njit
from ....utilities.types import FloatArray, ComplexArray

CmplxFltArray = Union[FloatArray, ComplexArray]
RadialFuncSolidDynamicType = Tuple[CmplxFltArray, CmplxFltArray, CmplxFltArray, CmplxFltArray, CmplxFltArray, CmplxFltArray]
RadialFuncLiquidDynamicType = Tuple[CmplxFltArray, CmplxFltArray, CmplxFltArray, CmplxFltArray]


@njit(cacheable=True)
def radial_derivatives_solid_general(radius : FloatArray, radial_functions : RadialFuncSolidDynamicType,
                                     shear_modulus : CmplxFltArray, bulk_modulus : CmplxFltArray, density : FloatArray,
                                     gravity : FloatArray, frequency: FloatArray,
                                     order_l: int = 2) -> RadialFuncSolidDynamicType:
    """ Calculates the radial derivative of the radial functions in the most general form - for solid layers.

    Allows for compressibility and dynamic tides.
    Tidal harmonic l is allowed to be an integer >= 2.

    References
    ----------
    KMN15; B15; TS72

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    radial_functions : RadialFuncSolidType
        Tuple of radial functions for a solid layer (y1, y2, y3, y4, y5, y6)
    shear_modulus : CmplxFltArray
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : CmplxFltArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at  at `radius` [kg m-3]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    radial_derivatives : RadialFuncSolidType
        The derivatives of the radial functions for a solid layer (dynamic assumption)

    """

    y1, y2, y3, y4, y5, y6 = radial_functions

    # Convert compressibility parameters (the first lame parameter can be complex)
    lame = (bulk_modulus - (2. / 3.) * shear_modulus)

    # Optimizations
    lame_2mu = 1. / (lame + 2. * shear_modulus)
    dynamic_term = frequency * frequency
    r_inverse = 1. / radius
    r2_inverse = r_inverse * r_inverse
    llp1 = order_l * (order_l + 1.)
    grav_term = 4. * pi * G * density

    # See Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    # # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    dy1 = \
        y1 * (-2. * lame * lame_2mu) * r_inverse + \
        y2 * lame_2mu + \
        y3 * llp1 * lame * lame_2mu * r_inverse

    dy2 = \
        y1 * (-dynamic_term * density + 12. * bulk_modulus * shear_modulus * lame_2mu * r2_inverse -
              4. * density * gravity * r_inverse) - \
        y2 * (4. * shear_modulus * lame_2mu * r_inverse) + \
        y3 * r_inverse * llp1 * (density * gravity - 6. * bulk_modulus * shear_modulus * lame_2mu * r_inverse) + \
        y4 * llp1 * r_inverse + \
        y5 * (order_l + 1.) * r_inverse * density - \
        y6 * density

    dy3 = \
        y1 * (-r_inverse) + \
        y3 * r_inverse + \
        y4 * (1. / shear_modulus)

    dy4 = \
        y1 * (density * gravity * r_inverse - 6. * bulk_modulus * shear_modulus * lame_2mu * r2_inverse) - \
        y2 * lame * lame_2mu * r_inverse + \
        y3 * (-dynamic_term * density + 4. * llp1 * (lame + shear_modulus) * shear_modulus * lame_2mu * r2_inverse -
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
def radial_derivatives_liquid_general(radius : FloatArray, radial_functions : RadialFuncLiquidDynamicType,
                                      bulk_modulus : CmplxFltArray, density : FloatArray,
                                      gravity : FloatArray, frequency: FloatArray,
                                      order_l: int = 2) -> RadialFuncLiquidDynamicType:
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
    radial_functions : RadialFuncLiquidDynamicType
        Tuple of radial functions for a solid layer (y1, y2, y5, y6)
    bulk_modulus : CmplxFltArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at  at `radius` [kg m-3]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.

    Returns
    -------
    radial_derivatives : RadialFuncLiquidDynamicType
        The derivatives of the radial functions for a liquid layer (dynamic assumption)

    """

    # For the dynamic version, y4 = 0 always in a liquid layer and y3 is defined by y1, y2, and y5 analytically
    y1, y2, y5, y6 = radial_functions


    # Convert compressibility parameters (the first lame parameter can be complex)
    # For the liquid layer it is assumed that the shear modulus is zero so the lame parameter simply
    #    equals the bulk modulus
    lame = bulk_modulus

    # Optimizations
    dynamic_term = frequency * frequency
    dynamic_term_inverse = 1. / dynamic_term
    r_inverse = 1. / radius
    r2_inverse = r_inverse * r_inverse
    llp1 = order_l * (order_l + 1.)
    dynamic_grav = (1. - order_l * gravity * dynamic_term_inverse * r_inverse)
    grav_term = 4. * pi * G * density


    # See Eqs. 11--14 in KMN15
    # # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    dy1 = \
        y1 * (-2. * r_inverse + llp1 * gravity * r2_inverse * dynamic_term_inverse) + \
        y2 * ((1. / lame) - llp1 * dynamic_term_inverse * r2_inverse / density) - \
        y5 * llp1 * r2_inverse * dynamic_term_inverse

    dy2 = \
        y1 * (-dynamic_term * density - 4. * density * gravity * r_inverse +
              llp1 * density * gravity**2 * dynamic_term_inverse * r2_inverse) - \
        y2 * llp1 * gravity * dynamic_term_inverse * r2_inverse + \
        y5 * (order_l + 1.) * density * r_inverse * dynamic_grav - \
        y6 * density

    dy5 = \
        y1 * grav_term - \
        y5 * (order_l + 1.) * r_inverse + \
        y6

    dy6 = \
        y1 * grav_term * (order_l + 1.) * dynamic_grav * r_inverse + \
        y2 * grav_term * llp1 * r2_inverse * dynamic_term_inverse / density + \
        y5 * grav_term * llp1 * r2_inverse * dynamic_term_inverse + \
        y6 * (order_l - 1.) * r_inverse

    return dy1, dy2, dy5, dy6