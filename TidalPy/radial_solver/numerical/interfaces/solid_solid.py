""" Functions to calculate the initial conditions for an overlying solid layer above another solid layer.

For solid-solid layer interfaces, all radial functions are continuous.

Since the solid solutions do not lose a y or an independent solution when moving from dynamic to static: Then the
interfaces between static-dynamic or dynamic-static will also be fully continuous.

In all likely-hood you probably will not want to use these as true interfaces. Instead, just combining all adjacent
"solid" layers into one super solid layer.

References
----------
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""
import numpy as np

from TidalPy.utilities.performance import njit


@njit(cacheable=False)
def both_dynamic(solid_layer_ys: np.ndarray) -> np.ndarray:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    solid layer. Assumes dynamic tides in both of the layers.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : np.ndarray
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and dynamic.

    Returns
    -------
    initial_solutions_solid : np.ndarray
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    initial_solutions = np.empty((3, 6), dtype=np.complex128)

    for solution in range(3):
        for yi in range(6):
            # Upper layer equals lower layer for all solutions and ys.
            initial_solutions[solution, yi] = solid_layer_ys[solution, yi]

    return initial_solutions


@njit(cacheable=False)
def static_dynamic(solid_layer_ys: np.ndarray) -> np.ndarray:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    solid layer. Assumes static tides in the lower layer and dynamic in the upper.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : np.ndarray
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and static.

    Returns
    -------
    initial_solutions_solid : np.ndarray
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    initial_solutions = np.empty((3, 6), dtype=np.complex128)

    for solution in range(3):
        for yi in range(6):
            # Upper layer equals lower layer for all solutions and ys.
            initial_solutions[solution, yi] = solid_layer_ys[solution, yi]

    return initial_solutions


@njit(cacheable=False)
def dynamic_static(solid_layer_ys: np.ndarray) -> np.ndarray:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    solid layer. Assumes static tides in the upper layer and dynamic in the lower.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : np.ndarray
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and dynamic.

    Returns
    -------
    initial_solutions_solid : np.ndarray
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    initial_solutions = np.empty((3, 6), dtype=np.complex128)

    for solution in range(3):
        for yi in range(6):
            # Upper layer equals lower layer for all solutions and ys.
            initial_solutions[solution, yi] = solid_layer_ys[solution, yi]

    return initial_solutions


@njit(cacheable=False)
def both_static(solid_layer_ys: np.ndarray) -> np.ndarray:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    solid layer. Assumes static tides in both of the layers.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : np.ndarray
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and static.

    Returns
    -------
    initial_solutions_solid : np.ndarray
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    initial_solutions = np.empty((3, 6), dtype=np.complex128)

    for solution in range(3):
        for yi in range(6):
            # Upper layer equals lower layer for all solutions and ys.
            initial_solutions[solution, yi] = solid_layer_ys[solution, yi]

    return initial_solutions
