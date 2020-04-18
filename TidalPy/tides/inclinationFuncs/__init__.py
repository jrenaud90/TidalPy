from .orderl2 import calc_inclination as calc_inclin_l2
from .orderl2 import calc_inclination_off as calc_inclin_l2_off
from .orderl3 import calc_inclination as calc_inclin_l3
from .orderl3 import calc_inclination_off as calc_inclin_l3_off
from .orderl4 import calc_inclination as calc_inclin_l4
from .orderl4 import calc_inclination_off as calc_inclin_l4_off
from .orderl5 import calc_inclination as calc_inclin_l5
from .orderl5 import calc_inclination_off as calc_inclin_l5_off
from .orderl5 import calc_inclination as calc_inclin_l6
from .orderl5 import calc_inclination_off as calc_inclin_l6_off
from .orderl5 import calc_inclination as calc_inclin_l7
from .orderl5 import calc_inclination_off as calc_inclin_l7_off

# Build Truncation Tables
#    Two different tables for the case that inclination (obliquity) is used or not. In the case that obliquity is not
#    used (or is always = 0) then it is much more efficient to use the "_off" version.
inclination_functions_on = (
    calc_inclin_l2,
    calc_inclin_l3,
    calc_inclin_l4,
    calc_inclin_l5,
    calc_inclin_l6,
    calc_inclin_l7
)

inclination_functions_off = (
    calc_inclin_l2_off,
    calc_inclin_l3_off,
    calc_inclin_l4_off,
    calc_inclin_l5_off,
    calc_inclin_l6_off,
    calc_inclin_l7_off
)

inclination_functions = {
    True: inclination_functions_on,
    False: inclination_functions_off
}