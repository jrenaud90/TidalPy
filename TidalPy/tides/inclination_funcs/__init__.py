from typing import Dict, Tuple

from ...utilities.types import FloatArray

InclinOutput = Dict[Tuple[int, int], FloatArray]

from .orderl2 import calc_inclination as calc_inclin_l2
from .orderl2 import calc_inclination_off as calc_inclin_l2_off
from .orderl3 import calc_inclination as calc_inclin_l3
from .orderl3 import calc_inclination_off as calc_inclin_l3_off
from .orderl4 import calc_inclination as calc_inclin_l4
from .orderl4 import calc_inclination_off as calc_inclin_l4_off
from .orderl5 import calc_inclination as calc_inclin_l5
from .orderl5 import calc_inclination_off as calc_inclin_l5_off
from .orderl6 import calc_inclination as calc_inclin_l6
from .orderl6 import calc_inclination_off as calc_inclin_l6_off
from .orderl7 import calc_inclination as calc_inclin_l7
from .orderl7 import calc_inclination_off as calc_inclin_l7_off

# Build Truncation Tables
#    Two different tables for the case that inclination (obliquity) is used or not. In the case that obliquity is not
#    used (or is always = 0) then it is much more efficient to use the "_off" version.
inclination_functions_on = {
    2: calc_inclin_l2,
    3: calc_inclin_l3,
    4: calc_inclin_l4,
    5: calc_inclin_l5,
    6: calc_inclin_l6,
    7: calc_inclin_l7
}

inclination_functions_off = {
    2: calc_inclin_l2_off,
    3: calc_inclin_l3_off,
    4: calc_inclin_l4_off,
    5: calc_inclin_l5_off,
    6: calc_inclin_l6_off,
    7: calc_inclin_l7_off
}

inclination_functions = {
    True: inclination_functions_on,
    False: inclination_functions_off
}


# @njit
def get_inclination_func(tidal_order_lvl: int = 2, inclination_nonzero: bool = True):
    inclination_functions_on_infunc = {
        2: calc_inclin_l2,
        3: calc_inclin_l3,
        4: calc_inclin_l4,
        5: calc_inclin_l5,
        6: calc_inclin_l6,
        7: calc_inclin_l7
    }

    inclination_functions_off_infunc = {
        2: calc_inclin_l2_off,
        3: calc_inclin_l3_off,
        4: calc_inclin_l4_off,
        5: calc_inclin_l5_off,
        6: calc_inclin_l6_off,
        7: calc_inclin_l7_off
    }

    inclination_functions_infunc = {
        True: inclination_functions_on_infunc,
        False: inclination_functions_off_infunc
    }

    return inclination_functions_infunc[inclination_nonzero][tidal_order_lvl]