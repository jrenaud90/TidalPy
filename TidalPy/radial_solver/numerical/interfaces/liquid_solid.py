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
from TidalPy.utilities.performance import njit


@njit(cacheable=False)
def both_dynamic(liquid_layer_ys: np.ndarray) -> np.ndarray:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    liquid layer. Assumes dynamic tides in both layers.

    References
    ----------
    Eqs. 148-149 in TS72

    Parameters
    ----------
    liquid_layer_ys : np.ndarray
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is liquid and dynamic.

    Returns
    -------
    initial_solutions_solid : np.ndarray
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    # For a dynamic solid layer there will be three independent solutions that we need an initial guess for.
    # For a dynamic solid layer there will be six y values used.
    initial_solutions = np.empty((3, 6), dtype=np.complex128)

    for solution in (0, 1, 2):
        if solution in (0, 1):
            # For a dynamic liquid layer there will be two independent solutions at the top of the layer
            initial_solutions[solution, 0] = liquid_layer_ys[solution, 0]
            initial_solutions[solution, 1] = liquid_layer_ys[solution, 1]
            initial_solutions[solution, 4] = liquid_layer_ys[solution, 2]
            initial_solutions[solution, 5] = liquid_layer_ys[solution, 3]

            # For solutions 1 and 2 y_3 and y_4 for the solid layer are zero
            initial_solutions[solution, 2] = 0.
            initial_solutions[solution, 3] = 0.
        else:
            # For the third solid solution all the y's are set to zero except y_3
            initial_solutions[solution, 0] = 0.
            initial_solutions[solution, 1] = 0.
            initial_solutions[solution, 2] = 1.
            initial_solutions[solution, 3] = 0.
            initial_solutions[solution, 4] = 0.
            initial_solutions[solution, 5] = 0.

    return initial_solutions


@njit(cacheable=False)
def static_dynamic(
        liquid_layer_ys: np.ndarray,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> np.ndarray:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    liquid layer. Assumes static tides in the lower layer and dynamic in the upper.

    References
    ----------
    Eqs. 20 in S74

    Parameters
    ----------
    liquid_layer_ys : np.ndarray
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
    initial_solutions_solid : np.ndarray
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    # For a dynamic solid layer there will be three independent solutions that we need an initial guess for.
    # For a dynamic solid layer there will be six y values used.
    initial_solutions = np.empty((3, 6), dtype=np.complex128)

    # Solution 1
    initial_solutions[0, 0] = 0.
    # y_2_sol = -rho * y_5_liq
    initial_solutions[0, 1] = -liquid_density * liquid_layer_ys[0, 0]
    initial_solutions[0, 2] = 0.
    initial_solutions[0, 3] = 0.
    # y_5_sol = y_5_liq
    initial_solutions[0, 4] = liquid_layer_ys[0, 0]
    # y_6_sol = y_7_liq + (4 pi G rho / g) y_5_liq
    initial_solutions[0, 5] = liquid_layer_ys[0, 1] + \
                                (4. * np.pi * G_to_use * liquid_density / interface_gravity) * liquid_layer_ys[0, 0]
    # Solution 2
    # y_1_sol = 1.
    initial_solutions[1, 0] = 1.
    # y_2_sol = rho * g * y_1_sol
    initial_solutions[1, 1] = liquid_density * interface_gravity * initial_solutions[1, 0]
    initial_solutions[1, 2] = 0.
    initial_solutions[1, 3] = 0.
    initial_solutions[1, 4] = 0.
    # y_6_sol = -4 pi G rho y_1_sol
    initial_solutions[1, 5] = -4. * np.pi * G_to_use * liquid_density * initial_solutions[1, 0]

    # Solution 3
    initial_solutions[2, 0] = 0.
    initial_solutions[2, 1] = 0.
    # y_3_sol = 1.
    initial_solutions[2, 2] = 1.
    initial_solutions[2, 3] = 0.
    initial_solutions[2, 4] = 0.
    initial_solutions[2, 5] = 0.

    return initial_solutions


@njit(cacheable=False)
def dynamic_static(liquid_layer_ys: np.ndarray) -> np.ndarray:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    liquid layer. Assumes static tides in the upper layer and dynamic in the lower.

    References
    ----------
    Eqs. 148-149 in TS72

    Parameters
    ----------
    liquid_layer_ys : np.ndarray
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is liquid and dynamic.

    Returns
    -------
    initial_solutions_solid : np.ndarray
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    # As far as I am aware, this should be the same as the both dynamic solutions.
    return both_dynamic(liquid_layer_ys)


@njit(cacheable=False)
def both_static(
        liquid_layer_ys: np.ndarray,
        interface_gravity: float, liquid_density: float, G_to_use: float = G
        ) -> np.ndarray:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    liquid layer. Assumes static tides in both the lower and upper layers.

    References
    ----------
    Eqs. 20 in S74

    Parameters
    ----------
    liquid_layer_ys : np.ndarray
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
    initial_solutions_solid : np.ndarray
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    # As far as I am aware, this should be the same as the static-dynamic solutions.
    return static_dynamic(liquid_layer_ys, interface_gravity, liquid_density, G_to_use=G_to_use)
