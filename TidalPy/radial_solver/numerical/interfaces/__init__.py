from typing import Tuple, Union

from .liquid_liquid import (both_dynamic as interface_LDy_LDy, both_static as interface_LSt_LSt,
                            dynamic_static as interface_LDy_LSt, static_dynamic as interface_LSt_LDy)
from .liquid_solid import (both_dynamic as interface_LDy_SDy, both_static as interface_LSt_SSt,
                           dynamic_static as interface_LDy_SSt, static_dynamic as interface_LSt_SDy)
from .solid_liquid import (both_dynamic as interface_SDy_LDy, both_static as interface_SSt_LSt,
                           dynamic_static as interface_SDy_LSt, static_dynamic as interface_SSt_LDy)
from .solid_solid import (both_dynamic as interface_SDy_SDy, both_static as interface_SSt_SSt,
                          dynamic_static as interface_SDy_SSt, static_dynamic as interface_SSt_SDy)
from TidalPy.constants import G

# Stored by is_lower_solid, is_upper_solid, is_lower_static, is_upper_static.
# Result is the correct function and a boolean if the function requires additional inputs.
known_interfaces = {
    # Solid-Solid
    (True, True, True, True)    : (interface_SSt_SSt, False),
    (True, True, True, False)   : (interface_SSt_SDy, False),
    (True, True, False, True)   : (interface_SDy_SSt, False),
    (True, True, False, False)  : (interface_SDy_SDy, False),

    # Liquid-Liquid
    (False, False, True, True)    : (interface_LSt_LSt, False),
    (False, False, True, False)   : (interface_LSt_LDy, True),
    (False, False, False, True)   : (interface_LDy_LSt, True),
    (False, False, False, False)  : (interface_LDy_LDy, False),

    # Solid-Liquid
    (True, False, True, True)   : (interface_SSt_LSt, True),
    (True, False, True, False)  : (interface_SSt_LDy, False),
    (True, False, False, True)  : (interface_SDy_LSt, True),
    (True, False, False, False) : (interface_SDy_LDy, False),

    # Liquid-Solid
    (False, True, True, True)   : (interface_LSt_SSt, True),
    (False, True, True, False)  : (interface_LSt_SDy, True),
    (False, True, False, True)  : (interface_LDy_SSt, False),
    (False, True, False, False) : (interface_LDy_SDy, False),
    }


def find_interface_func(
    lower_layer_is_solid: bool, lower_layer_is_static: bool,
    upper_layer_is_solid: bool, upper_layer_is_static: bool,
    static_liquid_density: float = None, interface_gravity: float = None, G_to_use: float = G
    ) -> Tuple[callable, Union[Tuple[()], Tuple[float, float, float]]]:
    """ Helper to find the correct interface function based on user provided assumptions.

    Used in `TidalPy.radial_solver.numerical`.

    Parameters
    ----------
    lower_layer_is_solid : bool
        If yes, the layer below the interface is assumed to be solid (rather than liquid).
    lower_layer_is_static : bool
        If yes, the layer below the interface is assumed to be static (rather than dynamic).
    upper_layer_is_solid : bool
        If yes, the layer above the interface is assumed to be solid (rather than liquid).
    upper_layer_is_static : bool
        If yes, the layer above the interface is assumed to be static (rather than dynamic).
    static_liquid_density : float = None
        Depending on the layer types and assumptions, the density of the static liquid layer may be required [kg m-3].
    interface_gravity : float = None
        Depending on the layer types and assumptions, the acceleration of gravity
        at the interface may be required [m s-2].
    G_to_use : float = G
        Gravitational constant (allowed to be different from standard G if non-dimensionalized)

    Returns
    -------
    interface_func : callable
        Interface function.
    inputs : Union[Tuple[()], Tuple[float, float, float]]
        If required, this will hold a tuple with additional required inputs for the interface function.

    """

    interface_func, needs_extra_input = known_interfaces[(lower_layer_is_solid,
                                                          upper_layer_is_solid,
                                                          lower_layer_is_static,
                                                          upper_layer_is_static)]

    inputs = tuple()
    if needs_extra_input:
        if static_liquid_density is None or interface_gravity is None:
            raise ValueError(f'{interface_func} requires additional interface parameters, none provided.')

        inputs = (interface_gravity, static_liquid_density, G_to_use)

    return interface_func, inputs
