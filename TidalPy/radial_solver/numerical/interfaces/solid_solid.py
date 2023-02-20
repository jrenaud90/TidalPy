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

from TidalPy.utilities.performance import njit, nbList

from ..initial.initial_solution_dynamic import SolidDynamicGuess
from ..initial.initial_solution_static import SolidStaticGuess

@njit(cacheable=True)
def both_dynamic(solid_layer_ys: SolidDynamicGuess) -> SolidDynamicGuess:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    solid layer. Assumes dynamic tides in both of the layers.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : SolidDynamicGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and dynamic.

    Returns
    -------
    initial_solutions_solid : SolidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    initial_solutions_solid = nbList([
        np.ascontiguousarray(solid_layer_ys[0]),
        np.ascontiguousarray(solid_layer_ys[1]),
        np.ascontiguousarray(solid_layer_ys[2])
        ])

    return initial_solutions_solid


@njit(cacheable=True)
def static_dynamic(solid_layer_ys: SolidStaticGuess) -> SolidDynamicGuess:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    solid layer. Assumes static tides in the lower layer and dynamic in the upper.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : SolidStaticGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and static.

    Returns
    -------
    initial_solutions_solid : SolidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    initial_solutions_solid = nbList([
        np.ascontiguousarray(solid_layer_ys[0]),
        np.ascontiguousarray(solid_layer_ys[1]),
        np.ascontiguousarray(solid_layer_ys[2])
        ])

    return initial_solutions_solid


@njit(cacheable=True)
def dynamic_static(solid_layer_ys: SolidDynamicGuess) -> SolidStaticGuess:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    solid layer. Assumes static tides in the upper layer and dynamic in the lower.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : SolidDynamicGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and dynamic.

    Returns
    -------
    initial_solutions_solid : SolidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    initial_solutions_solid = nbList([
        np.ascontiguousarray(solid_layer_ys[0]),
        np.ascontiguousarray(solid_layer_ys[1]),
        np.ascontiguousarray(solid_layer_ys[2])
        ])

    return initial_solutions_solid


@njit(cacheable=True)
def both_static(solid_layer_ys: SolidStaticGuess) -> SolidStaticGuess:
    """ Find the starting values for the radial functions at the bottom of a solid layer that is above a
    solid layer. Assumes static tides in both of the layers.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : SolidStaticGuess
        The solution for the radial functions in the layer below.
        This function assumes a lower layer that is solid and static.

    Returns
    -------
    initial_solutions_solid : SolidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the upper layer.
        For this function's assumptions, there will be three independent solutions for the upper layer.
    """

    initial_solutions_solid = nbList([
        np.ascontiguousarray(solid_layer_ys[0]),
        np.ascontiguousarray(solid_layer_ys[1]),
        np.ascontiguousarray(solid_layer_ys[2])
        ])

    return nbList(initial_solutions_solid)
