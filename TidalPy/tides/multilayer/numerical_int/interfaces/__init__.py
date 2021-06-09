from .....utilities.performance import njit

from .liquid_liquid import both_static as interface_LSt_LSt, both_dynamic as interface_LDy_LDy

from .liquid_solid import both_static as interface_LSt_SSt, both_dynamic as interface_LDy_SDy, \
    static_dynamic as interface_LSt_SDy, dynamic_static as interface_LDy_SSt

from .solid_liquid import both_static as interface_SSt_LSt, both_dynamic as interface_SDy_LDy, \
    static_dynamic as interface_SSt_LDy, dynamic_static as interface_SDy_LSt

from .solid_solid import both_static as interface_SSt_SSt, both_dynamic as interface_SDy_SDy, \
    static_dynamic as interface_SSt_SDy, dynamic_static as interface_SDy_SSt

# Stored by is_lower_solid, is_upper_solid, is_lower_static, is_upper_static
known_interfaces = {
    # Solid-Solid
    (True, True, True, True): (interface_SSt_SSt, False),
    (True, True, True, False): (interface_SSt_SDy, False),
    (True, True, False, True): (interface_SDy_SSt, False),
    (True, True, False, False): (interface_SDy_SDy, False),

    # Liquid-Liquid
    (False, False, True, True): (interface_LSt_LSt, False),
    (False, False, False, False): (interface_LDy_LDy, False),
    # Currently the mixed dynamic and static for L-L is not implemented.

    # Solid-Liquid
    (True, False, True, True): (interface_SSt_LSt, True),
    (True, False, True, False): (interface_SSt_LDy, False),
    (True, False, False, True): (interface_SDy_LSt, True),
    (True, False, False, False): (interface_SDy_LDy, False),

    # Liquid-Solid
    (False, True, True, True): (interface_LSt_SSt, True),
    (False, True, True, False): (interface_LSt_SDy, True),
    (False, True, False, True): (interface_LDy_SSt, False),
    (False, True, False, False): (interface_LDy_SDy, False),
}


def find_interface_func(lower_layer_is_solid: bool, lower_layer_is_static: bool,
                        upper_layer_is_solid: bool, upper_layer_is_static: bool,
                        liquid_density: float = None, interface_gravity: float = None):

    if not lower_layer_is_solid and not upper_layer_is_solid:
        # Both layers are liquid
        if (lower_layer_is_static and not upper_layer_is_static) or \
            (not lower_layer_is_static and upper_layer_is_static):
            raise NotImplementedError('Mixes of static and dynamic are not implemented for liquid-liquid interfaces.')

    interface_func, needs_extra_input = known_interfaces[(lower_layer_is_solid,
                                                          upper_layer_is_solid,
                                                          lower_layer_is_static,
                                                          upper_layer_is_static)]

    if needs_extra_input:
        if liquid_density is None or interface_gravity is None:
            raise ValueError(f'{interface_func} requires additional interface parameters, none provided.')

        @njit(cacheable=False)
        def interface_func_partial(lower_layer_ys):
            return interface_func(lower_layer_ys, interface_gravity, liquid_density)

        return interface_func_partial
    return interface_func