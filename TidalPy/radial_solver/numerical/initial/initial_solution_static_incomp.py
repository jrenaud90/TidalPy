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

from TidalPy.constants import G
from TidalPy.utilities.performance import njit, nbList
from TidalPy.utilities.types import ComplexArray

SolidStaticGuess = nbList[ComplexArray]
LiquidStaticGuess = nbList[ComplexArray]


@njit(cacheable=True)
def solid_guess_kamata(
    radius: float, shear_modulus: Union[float, complex],
    density: float,
    order_l: int = 2, G_to_use: float = G
    ) -> SolidStaticGuess:
    """ Calculate the initial guess at the bottom of a solid layer using the static and incompressible assumption.

    This function uses the Kamata et al (2015; JGR:P) equations (Eq. B17-B28).

    Using the dynamic assumption in a solid layer results in three independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), incompressible, and
       bulk and shear dissipation.

    References
    ----------
    KMN15 Eqs. B17-28

    Parameters
    ----------
    radius : float
        Radius where the radial functions are calculated. [m]
    shear_modulus : Union[float, complex]
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    density : float
        Density at `radius` [kg m-3]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : SolidStaticGuess
        The three independent solid guesses (sn1, sn2, sn3)

    """

    raise Exception('Not Implemented for the Incompressible Assumption')


@njit(cacheable=True)
def solid_guess_takeuchi(
    radius: float, shear_modulus: Union[float, complex],
    density: float,
    order_l: int = 2, G_to_use: float = G
    ) -> SolidStaticGuess:
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
    density : float
        Density at `radius` [kg m-3]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : SolidStaticGuess
        The three independent solid guesses (sn1, sn2, sn3)

    """

    # TODO
    raise Exception('Not Implemented for the Incompressible Assumption')


@njit(cacheable=True)
def liquid_guess_saito(radius: float, order_l: int = 2, G_to_use: float = G) -> LiquidStaticGuess:
    """ Calculate the initial guess at the bottom of a liquid layer using the static assumption.

    This function uses the Saito 1974 equations (Eq. 19).

    Using the static assumption in a liquid layer results in one independent solutions for the radial derivative.

    These independent solution allow for a general tidal harmonic l, for static tides (w = 0).
    However, compressibility and all dissipation dependence is lost due to no dependence on bulk or shear moduli.


    References
    ----------
    S74

    Parameters
    ----------
    radius : float
        Radius where the radial functions are calculated. [m]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : LiquidStaticGuess
        The one independent liquid guess (sn1)

    """

    # See Eq. 19 in Saito 1974
    # # y5 solutions
    y5_s1 = radius**order_l

    # # y7 solutions
    y7_s1 = 2. * (order_l - 1.) * radius**(order_l - 1.)

    # Since there is no bulk or shear dependence then the y's in this function will be strictly real
    tidaly_s1 = np.empty(2, dtype=np.complex128)
    tidaly_s1[0] = y5_s1
    tidaly_s1[1] = y7_s1

    # There is only one solution for the static liquid layer initial condition. However, that outputting a single value
    # does not match the form of the other initial condition functions. So we will wrap the value in a tuple.

    return nbList([tidaly_s1])