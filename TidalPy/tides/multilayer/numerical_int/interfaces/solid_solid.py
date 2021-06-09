""" Functions to calculate the initial conditions for an overlying solid layer above another solid layer.

References
----------
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

from ..initial_solution_dynamic import SolidDynamicGuess
from ..initial_solution_static import SolidStaticGuess
from .....utilities.performance import njit


# For solid-solid layer interfaces all of the radial functions are continuous.
# Since the solid solutions do not lose a y or an independent solution when moving from dynamic to static: Then the
#    interfaces between static-dynamic or dynamic-static will also fully continuous.
# TODO: Is this true for layers that experience a density jump/drop?
# In all likely hood you probably will not want to use these as true interfaces. Instead just combining all adjacent
#    "solid" layers into one super layer.

@njit(cacheable=True)
def both_dynamic(solid_layer_ys: SolidDynamicGuess) -> SolidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a solid layer that is above another
    solid surface. Assumes dynamic tides in both solid layers.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : SolidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic solid lower layer

    Returns
    -------
    base_solid_ys : SolidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be three independent solutions.
    """

    base_solid_ys = (
        solid_layer_ys[0][:, -1],
        solid_layer_ys[1][:, -1],
        solid_layer_ys[2][:, -1]
    )
    return base_solid_ys


@njit(cacheable=True)
def static_dynamic(solid_layer_ys: SolidStaticGuess) -> SolidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a solid layer that is above another
    solid surface. Assumes the lower layer is static and the upper is dynamic.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : SolidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static solid lower layer

    Returns
    -------
    base_solid_ys : SolidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be three independent solutions.
    """

    base_solid_ys = (
        solid_layer_ys[0][:, -1],
        solid_layer_ys[1][:, -1],
        solid_layer_ys[2][:, -1]
    )
    return base_solid_ys


@njit(cacheable=True)
def dynamic_static(solid_layer_ys: SolidDynamicGuess) -> SolidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a solid layer that is above another
    solid surface. Assumes the upper layer is static and the lower is dynamic.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : SolidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic solid lower layer

    Returns
    -------
    base_solid_ys : SolidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be three independent solutions.
    """

    base_solid_ys = (
        solid_layer_ys[0][:, -1],
        solid_layer_ys[1][:, -1],
        solid_layer_ys[2][:, -1]
    )
    return base_solid_ys


@njit(cacheable=True)
def both_static(solid_layer_ys: SolidStaticGuess) -> SolidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a solid layer that is above another
    solid surface. Assumes the upper and lower solid layers are both static.

    References
    ----------
    TS72

    Parameters
    ----------
    solid_layer_ys : SolidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static solid lower layer

    Returns
    -------
    base_solid_ys : SolidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be three independent solutions.
    """

    base_solid_ys = (
        solid_layer_ys[0][:, -1],
        solid_layer_ys[1][:, -1],
        solid_layer_ys[2][:, -1]
    )
    return base_solid_ys
