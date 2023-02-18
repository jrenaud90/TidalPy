""" Functions to calculate the initial conditions for an overlying solid layer above a liquid layer.

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
def both_dynamic(liquid_layer_ys: LiquidDynamicGuess) -> SolidDynamicGuess:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    liquid layer. Assumes dynamic tides in both layers.

    References
    ----------
    Eqs. 148-149 in TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is liquid and dynamic.

    Returns
    -------
    initial_solutions_solid : SolidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    # For a dynamic solid layer there will be three independent solutions that we need an initial guess for.
    # For a dynamic solid layer there will be six y values used.
    initial_solutions_solid = nbList([
        np.empty(6, dtype=np.complex128),
        np.empty(6, dtype=np.complex128),
        np.empty(6, dtype=np.complex128)
        ])

    for solution in (0, 1, 2):
        solution_ys = initial_solutions_solid[solution]

        if solution in (0, 1):
            # For a dynamic liquid layer there will be two independent solutions at the top of the layer
            liquid_sol = liquid_layer_ys[solution]
            solution_ys[0] = liquid_sol[0]
            solution_ys[1] = liquid_sol[1]
            solution_ys[4] = liquid_sol[2]
            solution_ys[5] = liquid_sol[3]

            # For solutions 1 and 2 y_3 and y_4 for the solid layer are zero
            solution_ys[2] = 0.
            solution_ys[3] = 0.
        else:
            # For the third solid solution all the y's are set to zero (already set by np.zeros) except y_3
            solution_ys[0] = 0.
            solution_ys[1] = 0.
            solution_ys[2] = 1.
            solution_ys[3] = 0.
            solution_ys[4] = 0.
            solution_ys[5] = 0.

    return initial_solutions_solid


@njit(cacheable=True)
def static_dynamic(
        liquid_layer_ys: LiquidStaticGuess,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> SolidDynamicGuess:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    liquid layer. Assumes static tides in the lower layer and dynamic in the upper.

    References
    ----------
    Eqs. 20 in S74

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
    initial_solutions_solid : SolidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    # For a dynamic solid layer there will be three independent solutions that we need an initial guess for.
    # For a dynamic solid layer there will be six y values used.
    initial_solutions_solid = nbList([
        np.empty(6, dtype=np.complex128),
        np.empty(6, dtype=np.complex128),
        np.empty(6, dtype=np.complex128)
        ])

    # For a static liquid layer there will be one independent solutions at the top of the layer
    #   however for consistency this single solution is still stored in a tuple (of size 1) so we need to pull it out.
    liquid_sol = liquid_layer_ys[0]

    for solution in range(3):
        solution_ys = initial_solutions_solid[solution]

        if solution == 0:
            solution_ys[0] = 0.
            # y_2_sol = -rho * y_5_liq
            solution_ys[1] = -liquid_density * liquid_sol[0]
            solution_ys[2] = 0.
            solution_ys[3] = 0.
            # y_5_sol = y_5_liq
            solution_ys[4] = liquid_sol[0]
            # y_6_sol = y_7_liq + (4 pi G rho / g) y_5_liq
            solution_ys[5] = liquid_sol[1] + \
                                 (4. * np.pi * G_to_use * liquid_density / interface_gravity) * liquid_sol[0]
        elif solution == 1:
            # y_1_sol = 1.
            solution_ys[0] = 1.
            # y_2_sol = rho * g * y_1_sol
            solution_ys[1] = liquid_density * interface_gravity * solution_ys[0]
            solution_ys[2] = 0.
            solution_ys[3] = 0.
            solution_ys[4] = 0.
            # y_6_sol = -4 pi G rho y_1_sol
            solution_ys[5] = -4. * np.pi * G_to_use * liquid_density * solution_ys[0]

        else:
            solution_ys[0] = 0.
            solution_ys[1] = 0.
            # y_3_sol = 1.
            solution_ys[2] = 1.
            solution_ys[3] = 0.
            solution_ys[4] = 0.
            solution_ys[5] = 0.

    return initial_solutions_solid


@njit(cacheable=True)
def dynamic_static(liquid_layer_ys: LiquidDynamicGuess) -> SolidStaticGuess:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    liquid layer. Assumes static tides in the upper layer and dynamic in the lower.

    References
    ----------
    Eqs. 148-149 in TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is liquid and dynamic.

    Returns
    -------
    initial_solutions_solid : SolidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    # As far as I am aware, this should be the same as the both dynamic solutions.
    return both_dynamic(liquid_layer_ys)


@njit(cacheable=True)
def both_static(
        liquid_layer_ys: LiquidStaticGuess,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> SolidStaticGuess:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    liquid layer. Assumes static tides in both the lower and upper layers.

    References
    ----------
    Eqs. 20 in S74

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
    initial_solutions_solid : SolidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    # As far as I am aware, this should be the same as the static-dynamic solutions.
    return static_dynamic(liquid_layer_ys, interface_gravity, liquid_density, G_to_use=G_to_use)
