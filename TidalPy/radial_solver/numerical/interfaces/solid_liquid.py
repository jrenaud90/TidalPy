""" Functions to calculate the initial conditions for an overlying liquid layer above a solid layer.

Note that there are 2 or 4 y's for the liquid layers and 6 for solid layers.
* y_sol_1 = sol index 0 <--> y_liq_1 = liq index 0 (for dynamic tides)
* y_sol_2 = sol index 1 <--> y_liq_2 = liq index 1 (for dynamic tides)
* y_sol_5 = sol index 4 <--> y_liq_5 = liq index 2 (for dynamic tides)
* y_sol_6 = sol index 5 <--> y_liq_6 = liq index 3 (for dynamic tides)
* y_sol_5 = sol index 4 <--> y_liq_5 = liq index 0 (for static tides)
* y_sol_6 = sol index 5 <--> y_liq_7 = liq index 1 (for static tides)

References
----------
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import njit, nbList

from ..initial.initial_solution_dynamic import LiquidDynamicGuess, SolidDynamicGuess
from ..initial.initial_solution_static import LiquidStaticGuess, SolidStaticGuess


@njit(cacheable=True)
def both_dynamic(solid_layer_ys: SolidDynamicGuess) -> LiquidDynamicGuess:
    """ Find the starting values for the radial functions at the bottom of a liquid layer that is above a
    solid layer. Assumes dynamic tides in both layers.

    References
    ----------
    Eqs. 140-144 in TS72

    Parameters
    ----------
    solid_layer_ys : SolidDynamicGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and dynamic.

    Returns
    -------
    initial_solutions_liquid : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be two independent solutions for the upper layer.
    """

    # For a dynamic liquid layer there will be two independent solutions
    # For a dynamic liquid layer there will be 4 y values used.
    initial_solutions_liquid = nbList([
        np.empty(4, dtype=np.complex128),
        np.empty(4, dtype=np.complex128),
        ])

    solid_layer_solution2 = solid_layer_ys[2]

    for solution in (0, 1):
        solution_ys = initial_solutions_liquid[solution]
        lower_layer_top_solution = solid_layer_ys[solution]

        # Solve for y^liq_1, y^liq_2, y^liq_5, y^liq_6 (TS72 Eq. 143)
        #    Note that the liquid solution does not have y_3, y_4 which are index 2, 3 for solid solution.
        solution_frac  = lower_layer_top_solution[3] / solid_layer_solution2[3]
        solution_ys[0] = lower_layer_top_solution[0] - solution_frac * solid_layer_solution2[0]
        solution_ys[1] = lower_layer_top_solution[1] - solution_frac * solid_layer_solution2[1]
        solution_ys[2] = lower_layer_top_solution[4] - solution_frac * solid_layer_solution2[4]
        solution_ys[3] = lower_layer_top_solution[5] - solution_frac * solid_layer_solution2[5]

    return initial_solutions_liquid


@njit(cacheable=True)
def static_dynamic(solid_layer_ys: SolidStaticGuess) -> LiquidDynamicGuess:
    """ Find the starting values for the radial functions at the bottom of a liquid layer that is above a
    solid layer. Assumes static tides in the lower layer and dynamic in the upper.

    References
    ----------
    Eqs. 140-144 in TS72

    Parameters
    ----------
    solid_layer_ys : SolidStaticGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and static.

    Returns
    -------
    initial_solutions_liquid : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be two independent solutions for the upper layer.
    """

    # As far as I am aware, this should work the same as the dynamic in both layers.
    return both_dynamic(solid_layer_ys)


@njit(cacheable=True)
def dynamic_static(
        solid_layer_ys: SolidDynamicGuess,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> LiquidStaticGuess:
    """ Find the starting values for the radial functions at the bottom of a liquid layer that is above a
    solid layer. Assumes static tides in the upper layer and dynamic in the lower.

    References
    ----------
    Eq. 21 in S74

    Parameters
    ----------
    solid_layer_ys : SolidDynamicGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and dynamic.
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

    # For a static liquid layer there will be one independent solution with 2 y's
    initial_solutions_liquid = np.empty(2, dtype=np.complex128)

    # Pull out ys
    # # Solution 1
    lower_s1y1 = solid_layer_ys[0][0]
    lower_s1y2 = solid_layer_ys[0][1]
    lower_s1y3 = solid_layer_ys[0][2]
    lower_s1y4 = solid_layer_ys[0][3]
    lower_s1y5 = solid_layer_ys[0][4]
    lower_s1y6 = solid_layer_ys[0][5]
    # # Solution 2
    lower_s2y1 = solid_layer_ys[1][0]
    lower_s2y2 = solid_layer_ys[1][1]
    lower_s2y3 = solid_layer_ys[1][2]
    lower_s2y4 = solid_layer_ys[1][3]
    lower_s2y5 = solid_layer_ys[1][4]
    lower_s2y6 = solid_layer_ys[1][5]
    # # Solution 3
    lower_s3y1 = solid_layer_ys[2][0]
    lower_s3y2 = solid_layer_ys[2][1]
    lower_s3y3 = solid_layer_ys[2][2]
    lower_s3y4 = solid_layer_ys[2][3]
    lower_s3y5 = solid_layer_ys[2][4]
    lower_s3y6 = solid_layer_ys[2][5]

    y4_frac_1 = -lower_s1y4 / lower_s3y4
    y4_frac_2 = -lower_s2y4 / lower_s3y4

    # gamma_j = (y_2j - f_j y_23) - rho( g(y_1j - f_j y_13) - (y_5j - f_j y_53))
    gamma_1 = (lower_s1y2 + y4_frac_1 * lower_s3y2) - \
              liquid_density * (interface_gravity * (lower_s1y1 + y4_frac_1 * lower_s3y1) -
                                (lower_s1y5 + y4_frac_1 * lower_s3y5))
    gamma_2 = (lower_s2y2 + y4_frac_2 * lower_s3y2) - \
              liquid_density * (interface_gravity * (lower_s2y1 + y4_frac_2 * lower_s3y1) -
                                (lower_s2y5 + y4_frac_2 * lower_s3y5))

    # Set the first coefficient to 1. It will be solved for later on during the collapse phase.
    coeff_1 = 1.
    # The other two coefficients are related to 1 via...
    coeff_2 = -(gamma_1 / gamma_2) * coeff_1
    coeff_3 = y4_frac_1 * coeff_1 + y4_frac_2 * coeff_2

    y_7_const = (4. * np.pi * G_to_use / interface_gravity)
    y_7_IC_0 = lower_s1y6 + y_7_const * lower_s1y2
    y_7_IC_1 = lower_s2y6 + y_7_const * lower_s2y2
    y_7_IC_2 = lower_s3y6 + y_7_const * lower_s3y2

    # y^liq_5 = C^sol_1 * y^sol_5,1 + C^sol_2 * y^sol_5,2 + C^sol_3 * y^sol_5,3
    initial_solutions_liquid[0] = coeff_1 * lower_s1y5 + coeff_2 * lower_s2y5 + coeff_3 * lower_s3y5
    # y^liq_7 = C^sol_1 * y^sol_7,1 + C^sol_2 * y^sol_7,2 + C^sol_3 * y^sol_7,3
    initial_solutions_liquid[1] = coeff_1 * y_7_IC_0 + coeff_2 * y_7_IC_1 + coeff_3 * y_7_IC_2

    return nbList([initial_solutions_liquid])


@njit(cacheable=True)
def both_static(
        solid_layer_ys: SolidStaticGuess,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> LiquidStaticGuess:
    """ Find the starting values for the radial functions at the bottom of a liquid layer that is above a
    solid layer. Assumes static tides in both the upper and lower layers.

    References
    ----------
    Eq. 20 in S74

    Parameters
    ----------
    solid_layer_ys : SolidStaticGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and static.
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

    # As far as I am aware, this should work the same as the dynamic-static function.
    return dynamic_static(solid_layer_ys, interface_gravity, liquid_density, G_to_use=G_to_use)
