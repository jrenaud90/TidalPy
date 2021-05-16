""" Functions to calculate the initial conditions for an overlying liquid layer above another liquid layer.

References
----------
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

from ..initial_solution_dynamic import LiquidDynamicGuess
from ..initial_solution_static import LiquidStaticGuess
from .....utilities.performance import njit


# For liquid-liquid layer interfaces all of the radial functions are continuous expect for if you are moving
#    from a dynamic to static or vice versa.
# TODO: Is this true for layers that experience a density jump/drop?

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

    base_liquid_ys = (
        liquid_layer_ys[0][:, -1],
        liquid_layer_ys[1][:, -1]
    )
    return base_liquid_ys


# @njit(cacheable=True)
def static_dynamic(liquid_layer_ys: LiquidStaticGuess) -> LiquidDynamicGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid surface. Assumes dynamic tides in the upper layer and static tides in the lower.

    References
    ----------
    TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidStaticGuess
        The solution for the radial functions in the layer below, this function assumes a static liquid lower layer

    Returns
    -------
    base_liquid_ys : LiquidDynamicGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be two independent solutions.
    """

    # TODO: Not Implemented
    raise NotImplementedError


# @njit(cacheable=True)
def dynamic_static(liquid_layer_ys: LiquidDynamicGuess) -> LiquidStaticGuess:
    """ Calculated the starting values for the radial functions at the bottom of a liquid layer that is above another
    liquid surface. Assumes dynamic tides in the upper layer and static tides in the lower.

    References
    ----------
    TS72

    Parameters
    ----------
    liquid_layer_ys : LiquidDynamicGuess
        The solution for the radial functions in the layer below, this function assumes a dynamic liquid lower layer

    Returns
    -------
    base_liquid_ys : LiquidStaticGuess
        The base (initial) solutions used to calculate the radial functions in the solid layer.
        For the assumptions used in this model there will be one independent solutions.
    """

    # TODO: Not Implemented
    raise NotImplementedError


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

    base_liquid_ys = liquid_layer_ys[:, -1]

    return base_liquid_ys
