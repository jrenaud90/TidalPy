""" Functions to calculate the initial conditions for an overlying liquid layer above a solid layer.

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
def both_dynamic(solid_layer_ys: SolidDynamicGuess) -> LiquidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above a solid
    surface. Assumes dynamic tides in both the liquid and solid layers.

    References
    ----------
    Eqs. 140-144 in TS72

    Parameters
    ----------
    solid_layer_ys : SolidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic solid lower layer

    Returns
    -------
    base_liquid_ys : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the liquid layer.
        For the assumptions used in this model there will be two independent solutions.
    """

    # For a dynamic liquid layer there will be two independent solutions
    base_liquid_y_list = list()

    for solution in range(2):
        lower_layer_top_solution = solid_layer_ys[solution][:, -1]
        # For a dynamic liquid layer there will be four y values used.
        solution_values = np.zeros(4, dtype=lower_layer_top_solution.dtype)

        # Solve for y_1, y_2, y_5, y_6
        #    Note that the liquid solution does not have y_3, y_4 which are index 2, 3 for solid solution.
        solution_frac = lower_layer_top_solution[3] / solid_layer_ys[2][3, -1]
        solution_values[0] = lower_layer_top_solution[0] - solution_frac * solid_layer_ys[2][0, -1]
        solution_values[1] = lower_layer_top_solution[1] - solution_frac * solid_layer_ys[2][1, -1]
        solution_values[2] = lower_layer_top_solution[4] - solution_frac * solid_layer_ys[2][4, -1]
        solution_values[3] = lower_layer_top_solution[5] - solution_frac * solid_layer_ys[2][5, -1]

        base_liquid_y_list.append(solution_values)

    base_liquid_ys = (base_liquid_y_list[0], base_liquid_y_list[1])
    return base_liquid_ys


@njit(cacheable=True)
def static_dynamic(solid_layer_ys: SolidStaticGuess) -> LiquidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above a solid
    surface. Assumes dynamic tides in the liquid layer and static tides in solid layers.

    References
    ----------
    Eqs. 140-144 in TS72

    Parameters
    ----------
    solid_layer_ys : SolidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static solid lower layer

    Returns
    -------
    base_liquid_ys : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the liquid layer.
        For the assumptions used in this model there will be two independent solutions.
    """

    # As far as I am aware, this should work the same as the dynamic in both layers.
    return both_dynamic(solid_layer_ys)


@njit(cacheable=True)
def dynamic_static(solid_layer_ys: SolidDynamicGuess,
                   interface_gravity: float, liquid_density: float, G_to_use: float = G) -> LiquidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above a solid
    surface. Assumes dynamic tides in the lower solid layer and static tides in the liquid layer.

    References
    ----------
    Eq. 20 in S74

    Parameters
    ----------
    solid_layer_ys : SolidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic solid lower layer
    interface_gravity : float
        Acceleration due to gravity at the interface [m s-2]
    liquid_density : float
        The density at the base of the liquid layer [kg m-3]
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    base_liquid_ys : LiquidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the liquid layer.
        For the assumptions used in this model there will be one independent solution.
    """

    # For a static liquid layer there will be one independent solution with 2 y's
    base_liquid_ys = np.zeros(2, dtype=solid_layer_ys[0].dtype)

    y4_frac_1 = solid_layer_ys[0][3, -1] / solid_layer_ys[2][3, -1]
    y4_frac_2 = solid_layer_ys[1][3, -1] / solid_layer_ys[2][3, -1]

    # gamma_j = (y_2j - f_j y_23) - rho( g(y_1j - f_j y_13) - (y_5j - f_j y_53))
    gamma_1 = (solid_layer_ys[0][1, -1] - y4_frac_1 * solid_layer_ys[2][1, -1]) - \
              liquid_density * (interface_gravity * (solid_layer_ys[0][0, -1] - y4_frac_1 * solid_layer_ys[2][0, -1]) -
                                (solid_layer_ys[0][4, -1] - y4_frac_1 * solid_layer_ys[2][4, -1]))
    gamma_2 = (solid_layer_ys[1][1, -1] - y4_frac_2 * solid_layer_ys[2][1, -1]) - \
              liquid_density * (interface_gravity * (solid_layer_ys[1][0, -1] - y4_frac_2 * solid_layer_ys[2][0, -1]) -
                                (solid_layer_ys[1][4, -1] - y4_frac_2 * solid_layer_ys[2][4, -1]))

    base_liquid_ys[0] = solid_layer_ys[0][4, -1] - (gamma_1 / gamma_2) * solid_layer_ys[1][4, -1] - \
                        (y4_frac_1 - (gamma_1 / gamma_2) * y4_frac_2) * solid_layer_ys[2][4, -1]

    y_7_IC_0 = solid_layer_ys[0][5, -1] + (4. * np.pi * G_to_use / interface_gravity) * solid_layer_ys[0][1, -1]
    y_7_IC_1 = solid_layer_ys[1][5, -1] + (4. * np.pi * G_to_use / interface_gravity) * solid_layer_ys[1][1, -1]
    y_7_IC_2 = solid_layer_ys[2][5, -1] + (4. * np.pi * G_to_use / interface_gravity) * solid_layer_ys[2][1, -1]

    base_liquid_ys[1] = y_7_IC_0 - (gamma_1 / gamma_2) * y_7_IC_1 - \
                        (y4_frac_1 - (gamma_1 / gamma_2) * y4_frac_2) * y_7_IC_2

    return base_liquid_ys


@njit(cacheable=True)
def both_static(solid_layer_ys: SolidStaticGuess,
                interface_gravity: float, liquid_density: float, G_to_use: float = G) -> LiquidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above a solid
    surface. Assumes static tides in both the liquid and solid layers.

    References
    ----------
    Eq. 20 in S74

    Parameters
    ----------
    solid_layer_ys : SolidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static solid lower layer
    interface_gravity : float
        Acceleration due to gravity at the interface [m s-2]
    liquid_density : float
        The density at the base of the liquid layer [kg m-3]
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    base_liquid_ys : LiquidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the liquid layer.
        For the assumptions used in this model there will be one independent solution.
    """

    # As far as I am aware, this should work the same as the dynamic-static function.
    return dynamic_static(solid_layer_ys, interface_gravity, liquid_density, G_to_use=G_to_use)
