""" Functions to calculate the initial conditions for an overlying solid layer above a liquid layer.

References
----------
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

import numpy as np

from ..initial_solution_dynamic import SolidDynamicGuess, LiquidDynamicGuess
from ..initial_solution_static import SolidStaticGuess, LiquidStaticGuess
from .....constants import G
from .....utilities.performance import njit


# Note that there are 2 or 4 y's for the liquid layers and 6 for solid layers.
#    y_sol_1 = sol index 0 <--> y_liq_1 = liq index 0 (for dynamic tides)
#    y_sol_2 = sol index 1 <--> y_liq_2 = liq index 1 (for dynamic tides)
#    y_sol_5 = sol index 4 <--> y_liq_5 = liq index 2 (for dynamic tides)
#    y_sol_6 = sol index 5 <--> y_liq_6 = liq index 3 (for dynamic tides)
#    y_sol_5 = sol index 4 <--> y_liq_5 = liq index 0 (for static tides)
#    y_sol_6 = sol index 5 <--> y_liq_7 = liq index 1 (for static tides)

@njit(cacheable=True)
def both_dynamic(liquid_layer_ys: LiquidDynamicGuess) -> SolidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a solid layer that is above a liquid
    surface. Assumes dynamic tides in both the liquid and solid layers.

    References
    ----------
    Eqs. 148-149 in TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic liquid lower layer

    Returns
    -------
    base_solid_ys : SolidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the liquid layer.
        For the assumptions used in this model there will be three independent solutions.
    """

    # For a dynamic solid layer there will be three independent solutions that we need an initial guess for.
    base_solid_list = list()

    for solution in range(3):
        # For a dynamic solid layer there will be six y values used.
        solution_values = np.zeros(6, dtype=liquid_layer_ys[0][:, -1].dtype)

        if solution in (0, 1):
            # For a dynamic liquid layer there will be two independent solutions at the top of the layer
            liquid_sol = liquid_layer_ys[solution][:, -1]
            solution_values[0] = liquid_sol[0]
            solution_values[1] = liquid_sol[1]
            solution_values[4] = liquid_sol[2]
            solution_values[5] = liquid_sol[3]

            # For solutions 1 and 2 y_3 and y_4 for the solid layer are zero (already set by np.zeros)
        else:
            # For the third solid solution all the y's are set to zero (already set by np.zeros) except y_3
            solution_values[2] = 1.

        base_solid_list.append(solution_values)

    base_solid_ys = (base_solid_list[0], base_solid_list[1], base_solid_list[2])
    return base_solid_ys


@njit(cacheable=True)
def static_dynamic(liquid_layer_ys: LiquidStaticGuess,
                   interface_gravity: float, liquid_density: float, G_to_use: float = G) -> SolidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a solid layer that is above a liquid
    surface. Assumes dynamic tides in the solid layer and static tides in the lower liquid layer.

    References
    ----------
    Eqs. 20 in S74

    Parameters
    ----------
    liquid_layer_ys : LiquidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static liquid lower layer
    interface_gravity : float
        Acceleration due to gravity at the interface [m s-2]
    liquid_density : float
        The density at the top of the liquid layer (liquid's density) [kg m-3]
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    base_solid_ys : SolidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the liquid layer.
        For the assumptions used in this model there will be three independent solutions.
    """

    # For a dynamic solid layer there will be three independent solutions that we need an initial guess for.
    base_solid_list = list()

    # For a dynamic liquid layer there will be one independent solutions at the top of the layer
    liquid_sol = liquid_layer_ys[:, -1]

    for solution in range(3):
        # For a dynamic solid layer there will be six y values used.
        solution_values = np.zeros(6, dtype=liquid_sol.dtype)

        if solution == 0:
            # y_2_sol = -rho * y_5_liq
            solution_values[1] = -liquid_density * liquid_sol[0]
            # y_5_sol = y_5_liq
            solution_values[4] = liquid_sol[0]
            # y_6_sol = y_7_liq + (4 pi G rho / g) y_5_liq
            solution_values[5] = liquid_sol[1] + \
                                 (4. * np.pi * G_to_use * liquid_density / interface_gravity) * liquid_sol[0]
        elif solution == 1:
            # y_1_sol = 1.
            solution_values[0] = 1.
            # y_2_sol = rho * g * y_1_sol
            solution_values[1] = liquid_density * interface_gravity * solution_values[0]
            # y_6_sol = -4 pi G rho y_1_sol
            solution_values[5] = -4. * np.pi * G_to_use * liquid_density * solution_values[0]

        else:
            # y_3_sol = 1.
            solution_values[2] = 1.

        base_solid_list.append(solution_values)

    base_solid_ys = (base_solid_list[0], base_solid_list[1], base_solid_list[2])
    return base_solid_ys


@njit(cacheable=True)
def dynamic_static(liquid_layer_ys: LiquidDynamicGuess) -> SolidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a solid layer that is above a liquid
    surface. Assumes static tides in the solid layer and dynamic tides in the lower liquid layer.

    References
    ----------
    Eqs. 148-149 in TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic liquid lower layer

    Returns
    -------
    base_solid_ys : SolidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the liquid layer.
        For the assumptions used in this model there will be three independent solutions.
    """

    # As far as I am aware, this should be the same as the both dynamic solutions.
    return both_dynamic(liquid_layer_ys)


@njit(cacheable=True)
def both_static(liquid_layer_ys: LiquidStaticGuess,
                interface_gravity: float, liquid_density: float, G_to_use: float = G) -> SolidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a solid layer that is above a liquid
    surface. Assumes static tides in both layers.

    References
    ----------
    Eqs. 20 in S74

    Parameters
    ----------
    liquid_layer_ys : LiquidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static liquid lower layer
    interface_gravity : float
        Acceleration due to gravity at the interface [m s-2]
    liquid_density : float
        The density at the top of the liquid layer (liquid's density) [kg m-3]
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    base_solid_ys : SolidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the liquid layer.
        For the assumptions used in this model there will be three independent solutions.
    """

    # As far as I am aware, this should be the same as the static-dynamic solutions.
    return static_dynamic(liquid_layer_ys, interface_gravity, liquid_density, G_to_use=G_to_use)
