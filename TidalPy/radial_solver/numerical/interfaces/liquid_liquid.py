""" Functions to calculate the initial conditions for an overlying liquid layer above another liquid layer.

References
----------
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import njit, nbList

from ..initial_conditions.initial_solution_dynamic import LiquidDynamicGuess, SolidDynamicGuess
from ..initial_conditions.initial_solution_static import LiquidStaticGuess, SolidStaticGuess


# For liquid-liquid layer interfaces: all the radial functions are continuous expect for if you are moving
#    from a dynamic to static or vice versa.

@njit(cacheable=True)
def both_dynamic(liquid_layer_ys: LiquidDynamicGuess) -> LiquidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid surface. Assumes dynamic tides in both solid layers.

    References
    ----------
    TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic liquid lower layer

    Returns
    -------
    base_liquid_ys : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be two independent solutions.
    """

    base_liquid_ys = [
        np.ascontiguousarray(liquid_layer_ys[0][:, -1]),
        np.ascontiguousarray(liquid_layer_ys[1][:, -1])
        ]
    return nbList(base_liquid_ys)


@njit(cacheable=True)
def static_dynamic(
        liquid_layer_ys: LiquidStaticGuess,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> LiquidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid surface. Assumes dynamic tides in the upper layer and static tides in the lower.

    References
    ----------
    TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static liquid lower layer
    interface_gravity : float
        Acceleration due to gravity at the interface [m s-2]
    liquid_density : float
        The density at the top of the liquid layer (liquid's density) [kg m-3]
        For this method we assume that the density provided is for the _static_ layer.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    base_liquid_ys : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be two independent solutions.
    """

    # JPR decided to follow a similar approach as Eq. 20 in S74:
    #   Treat the lower static liquid as normal.
    #   The upper dynamic liquid layer is treated like the solid layer in Eq. 20 except
    #    that y_3 is undefined as is "set 3" solution mentioned in that text.

    # For a dynamic upper liquid layer there will be three independent solutions that we need an initial guess for.
    base_liquid_list = nbList()

    # For a static lower liquid layer there will be one independent solutions at the top of the layer
    #   however for consistency this single solution is still stored in a tuple (of size 1) so we need to pull it out.
    liquid_static_sol = liquid_layer_ys[0][:, -1]

    for solution in (0, 1):
        # For a dynamic liquid layer there will be four y values used.
        solution_values = np.empty(4, dtype=liquid_static_sol.dtype)

        if solution == 0:
            # y_1_dynamic = 0
            solution_values[0] = 0.
            # y_2_dynamic = -rho * y_5_static
            solution_values[1] = -liquid_density * liquid_static_sol[0]
            # y_5_dynamic = y_5_static
            solution_values[2] = liquid_static_sol[0]
            # y_6_dynamic = y_7_static + (4 pi G rho / g) y_5_static
            solution_values[3] = liquid_static_sol[1] + \
                                 (4. * np.pi * G_to_use * liquid_density / interface_gravity) * liquid_static_sol[0]
        else:
            # y_1_dynamic = 1.
            solution_values[0] = 1.
            # y_2_dynamic = rho * g * y_1_dynamic
            solution_values[1] = liquid_density * interface_gravity * solution_values[0]
            # y_5_dynamic = 0.
            solution_values[2] = 0.
            # y_6_dynamic = -4 pi G rho y_1_dynamic
            solution_values[3] = -4. * np.pi * G_to_use * liquid_density * solution_values[0]

        base_liquid_list.append(solution_values)

    return base_liquid_list


@njit(cacheable=True)
def dynamic_static(
        liquid_layer_ys: LiquidDynamicGuess,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> LiquidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid surface. Assumes dynamic tides in the upper layer and static tides in the lower.

    References
    ----------
    TS72
    S74

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic liquid lower layer
    interface_gravity : float
        Acceleration due to gravity at the interface [m s-2]
    liquid_density : float
        The density at the top of the liquid layer (liquid's density) [kg m-3]
        For this method we assume that the density provided is for the _static_ layer.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    base_liquid_ys : LiquidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be one independent solutions.
    """

    # JPR decided to follow a similar approach as Eq. 20 in S74:
    #   Treat the lower static liquid as normal.
    #   The upper dynamic liquid layer is treated like the solid layer in Eq. 20 except
    #    that y_3 is undefined as is "set 3" solution mentioned in that text.

    # For a static liquid layer there will be one independent solution with 2 y's
    base_liquid_ys = np.empty(2, dtype=liquid_layer_ys[0].dtype)

    # lambda_j = (y_2j - rho * ( g * y_1j - y_5j))
    lambda_1 = liquid_layer_ys[0][1, -1] - \
               liquid_density * (interface_gravity * liquid_layer_ys[0][0, -1] - liquid_layer_ys[0][4, -1])
    lambda_2 = liquid_layer_ys[1][1, -1] - \
               liquid_density * (interface_gravity * liquid_layer_ys[1][0, -1] - liquid_layer_ys[1][4, -1])

    # Set the first coefficient to 1. It will be solved for later on during the collapse phase.
    coeff_1 = 1.
    # The other coefficient which is related to 1 via...
    coeff_2 = -(lambda_1 / lambda_2) * coeff_1

    y_7_const = (4. * np.pi * G_to_use / interface_gravity)
    y_7_IC_0 = liquid_layer_ys[0][5, -1] + y_7_const * liquid_layer_ys[0][1, -1]
    y_7_IC_1 = liquid_layer_ys[1][5, -1] + y_7_const * liquid_layer_ys[1][1, -1]

    # y^liq_5 = C^sol_1 * y^sol_5,1 + C^sol_2 * y^sol_5,2
    base_liquid_ys[0] = coeff_1 * liquid_layer_ys[0][4, -1] + coeff_2 * liquid_layer_ys[1][4, -1]
    # y^liq_7 = C^sol_1 * y^sol_7,1 + C^sol_2 * y^sol_7,2
    base_liquid_ys[1] = coeff_1 * y_7_IC_0 + coeff_2 * y_7_IC_1

    return nbList([base_liquid_ys])


@njit(cacheable=True)
def both_static(liquid_layer_ys: LiquidStaticGuess) -> LiquidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid surface. Assumes static tides in both solid layers.

    References
    ----------
    S74

    Parameters
    ----------
    liquid_layer_ys : LiquidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static liquid lower layer

    Returns
    -------
    base_liquid_ys : LiquidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be one independent solution.
    """

    base_liquid_ys = liquid_layer_ys[0][:, -1]
    base_liquid_ys = nbList([np.ascontiguousarray(base_liquid_ys)])

    return base_liquid_ys
