""" Functions to calculate the initial conditions for an overlying liquid layer above another liquid layer.

For liquid-liquid layer interfaces: all the radial functions are continuous expect for if you are moving
from a dynamic layer to static layer or vice versa.

References
----------
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import njit, nbList

from ..initial.initial_solution_dynamic import LiquidDynamicGuess
from ..initial.initial_solution_static import LiquidStaticGuess


@njit(cacheable=True)
def both_dynamic(liquid_layer_ys: LiquidDynamicGuess) -> LiquidDynamicGuess:
    """ Find the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid layer. Assumes dynamic tides in both layers.

    References
    ----------
    TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is liquid and dynamic.

    Returns
    -------
    initial_solutions_liquid : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be two independent solutions for the upper layer.
    """

    initial_solutions_liquid = nbList([
        np.ascontiguousarray(liquid_layer_ys[0]),
        np.ascontiguousarray(liquid_layer_ys[1])
        ])
    return initial_solutions_liquid


@njit(cacheable=True)
def static_dynamic(
        liquid_layer_ys: LiquidStaticGuess,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> LiquidDynamicGuess:
    """ Find the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid layer. Assumes static tides in the lower layer and dynamic tides the upper.

    References
    ----------
    TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidStaticGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is liquid and static.
    interface_gravity : float
        Acceleration due to gravity at the interface [m s-2].
    liquid_density : float
        The density at the top of the liquid layer (liquid's density) [kg m-3].
        For this method we assume that the density provided is for the _static_ layer.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    initial_solutions_liquid : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be two independent solutions for the upper layer.
    """

    # JPR decided to follow a similar approach as Eq. 20 in S74:
    #   Treat the lower static liquid as normal.
    #   The upper dynamic liquid layer is treated like the solid layer in Eq. 20 except
    #    that y_3 is undefined as is "set 3" solution mentioned in that text.

    # For a static lower liquid layer there will be one independent solutions at the top of the layer
    #   however for consistency this single solution is still stored in a tuple (of size 1) so we need to pull it out.
    liquid_static_sol = liquid_layer_ys[0]

    # For an upper layer that is liquid and dynamic there will be two independent solutions that need an initial guess.
    # For a dynamic liquid layer there will be four y values used.
    initial_solutions_liquid = nbList([
        np.empty(4, dtype=np.complex128),
        np.empty(4, dtype=np.complex128)
        ])

    for solution in (0, 1):
        solution_ys = initial_solutions_liquid[solution]

        if solution == 0:
            # y_1_dynamic = 0
            solution_ys[0] = 0.
            # y_2_dynamic = -rho * y_5_static
            solution_ys[1] = -liquid_density * liquid_static_sol[0]
            # y_5_dynamic = y_5_static
            solution_ys[2] = liquid_static_sol[0]
            # y_6_dynamic = y_7_static + (4 pi G rho / g) y_5_static
            solution_ys[3] = liquid_static_sol[1] + \
                             (4. * np.pi * G_to_use * liquid_density / interface_gravity) * liquid_static_sol[0]
        else:
            # y_1_dynamic = 1.
            solution_ys[0] = 1.
            # y_2_dynamic = rho * g * y_1_dynamic
            solution_ys[1] = liquid_density * interface_gravity * solution_ys[0]
            # y_5_dynamic = 0.
            solution_ys[2] = 0.
            # y_6_dynamic = -4 pi G rho y_1_dynamic
            solution_ys[3] = -4. * np.pi * G_to_use * liquid_density * solution_ys[0]

    return initial_solutions_liquid


@njit(cacheable=True)
def dynamic_static(
        liquid_layer_ys: LiquidDynamicGuess,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> LiquidStaticGuess:
    """ Find the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid layer. Assumes static tides in the upper layer and dynamic tides the lower.

    References
    ----------
    TS72
    S74

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is liquid and dynamic.
    interface_gravity : float
        Acceleration due to gravity at the interface [m s-2].
    liquid_density : float
        The density at the top of the liquid layer (liquid's density) [kg m-3].
        For this method we assume that the density provided is for the _static_ layer.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    initial_solutions_liquid : LiquidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be one independent solutions for the upper layer.
    """

    # JPR decided to follow a similar approach as Eq. 20 in S74:
    #   Treat the lower static liquid as normal.
    #   The upper dynamic liquid layer is treated like the solid layer in Eq. 20 except
    #    that y_3 is undefined as is "set 3" solution mentioned in that text.

    # Pull out ys
    # # Solution 1
    lower_s1y1 = liquid_layer_ys[0][0]
    lower_s1y2 = liquid_layer_ys[0][1]
    lower_s1y5 = liquid_layer_ys[0][2]
    lower_s1y6 = liquid_layer_ys[0][3]
    # # Solution 2
    lower_s2y1 = liquid_layer_ys[1][0]
    lower_s2y2 = liquid_layer_ys[1][1]
    lower_s2y5 = liquid_layer_ys[1][2]
    lower_s2y6 = liquid_layer_ys[1][3]

    # For a static liquid layer there will be one independent solution with 2 y's
    initial_solutions_liquid = np.empty(2, dtype=np.complex128)

    # lambda_j = (y_2j - rho * ( g * y_1j - y_5j))
    lambda_1 = lower_s1y2 - liquid_density * (interface_gravity * lower_s1y1 - lower_s1y5)
    lambda_2 = lower_s2y2 - liquid_density * (interface_gravity * lower_s2y1 - lower_s2y5)

    # Set the first coefficient to 1. It will be solved for later on during the collapse phase.
    coeff_1 = 1.
    # The other coefficient which is related to 1 via...
    coeff_2 = -(lambda_1 / lambda_2) * coeff_1

    y_7_const = (4. * np.pi * G_to_use / interface_gravity)
    y_7_IC_0 = lower_s1y6 + y_7_const * lower_s1y2
    y_7_IC_1 = lower_s2y6 + y_7_const * lower_s2y2

    # y^liq(st)_5 = C^liq(dy)_1 * y^liq(dy)_5,1 + C^liq(dy)_2 * y^liq(dy)_5,2
    initial_solutions_liquid[0] = coeff_1 * lower_s1y5 + coeff_2 * lower_s2y5
    # y^liq(st)_7 = C^liq(dy)_1 * y^liq(dy)_7,1 + C^liq(dy)_2 * y^liq(dy)_7,2
    initial_solutions_liquid[1] = coeff_1 * y_7_IC_0 + coeff_2 * y_7_IC_1

    return nbList([initial_solutions_liquid])


@njit(cacheable=True)
def both_static(liquid_layer_ys: LiquidStaticGuess) -> LiquidStaticGuess:
    """ Find the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid layer. Assumes static tides in both the upper and lower layers.

    References
    ----------
    S74

    Parameters
    ----------
    liquid_layer_ys : LiquidStaticGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is liquid and static.

    Returns
    -------
    initial_solutions_liquid : LiquidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be one independent solutions for the upper layer.
    """

    initial_solutions_liquid = nbList([np.ascontiguousarray(liquid_layer_ys[0])])

    return initial_solutions_liquid
