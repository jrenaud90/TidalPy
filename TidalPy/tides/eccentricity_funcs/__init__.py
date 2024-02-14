from typing import Dict, TYPE_CHECKING

if TYPE_CHECKING:
    from ...utilities.types import FloatArray

EccenOutput = Dict[int, Dict[int, 'FloatArray']]

from .orderl2 import eccentricity_funcs_trunc2 as eccentricity_funcs_l2_trunc2
from .orderl2 import eccentricity_funcs_trunc4 as eccentricity_funcs_l2_trunc4
from .orderl2 import eccentricity_funcs_trunc6 as eccentricity_funcs_l2_trunc6
from .orderl2 import eccentricity_funcs_trunc8 as eccentricity_funcs_l2_trunc8
from .orderl2 import eccentricity_funcs_trunc10 as eccentricity_funcs_l2_trunc10
from .orderl2 import eccentricity_funcs_trunc12 as eccentricity_funcs_l2_trunc12
from .orderl2 import eccentricity_funcs_trunc14 as eccentricity_funcs_l2_trunc14
from .orderl2 import eccentricity_funcs_trunc16 as eccentricity_funcs_l2_trunc16
from .orderl2 import eccentricity_funcs_trunc18 as eccentricity_funcs_l2_trunc18
from .orderl2 import eccentricity_funcs_trunc20 as eccentricity_funcs_l2_trunc20
from .orderl2 import eccentricity_funcs_trunc22 as eccentricity_funcs_l2_trunc22

from .orderl3 import eccentricity_funcs_trunc2 as eccentricity_funcs_l3_trunc2
from .orderl3 import eccentricity_funcs_trunc4 as eccentricity_funcs_l3_trunc4
from .orderl3 import eccentricity_funcs_trunc6 as eccentricity_funcs_l3_trunc6
from .orderl3 import eccentricity_funcs_trunc8 as eccentricity_funcs_l3_trunc8
from .orderl3 import eccentricity_funcs_trunc10 as eccentricity_funcs_l3_trunc10
from .orderl3 import eccentricity_funcs_trunc12 as eccentricity_funcs_l3_trunc12
from .orderl3 import eccentricity_funcs_trunc14 as eccentricity_funcs_l3_trunc14
from .orderl3 import eccentricity_funcs_trunc16 as eccentricity_funcs_l3_trunc16
from .orderl3 import eccentricity_funcs_trunc18 as eccentricity_funcs_l3_trunc18
from .orderl3 import eccentricity_funcs_trunc20 as eccentricity_funcs_l3_trunc20

from .orderl4 import eccentricity_funcs_trunc2 as eccentricity_funcs_l4_trunc2
from .orderl4 import eccentricity_funcs_trunc4 as eccentricity_funcs_l4_trunc4
from .orderl4 import eccentricity_funcs_trunc6 as eccentricity_funcs_l4_trunc6
from .orderl4 import eccentricity_funcs_trunc8 as eccentricity_funcs_l4_trunc8
from .orderl4 import eccentricity_funcs_trunc10 as eccentricity_funcs_l4_trunc10
from .orderl4 import eccentricity_funcs_trunc12 as eccentricity_funcs_l4_trunc12
from .orderl4 import eccentricity_funcs_trunc14 as eccentricity_funcs_l4_trunc14
from .orderl4 import eccentricity_funcs_trunc16 as eccentricity_funcs_l4_trunc16
from .orderl4 import eccentricity_funcs_trunc18 as eccentricity_funcs_l4_trunc18
from .orderl4 import eccentricity_funcs_trunc20 as eccentricity_funcs_l4_trunc20

from .orderl5 import eccentricity_funcs_trunc2 as eccentricity_funcs_l5_trunc2
from .orderl5 import eccentricity_funcs_trunc4 as eccentricity_funcs_l5_trunc4
from .orderl5 import eccentricity_funcs_trunc6 as eccentricity_funcs_l5_trunc6
from .orderl5 import eccentricity_funcs_trunc8 as eccentricity_funcs_l5_trunc8
from .orderl5 import eccentricity_funcs_trunc10 as eccentricity_funcs_l5_trunc10
from .orderl5 import eccentricity_funcs_trunc12 as eccentricity_funcs_l5_trunc12
from .orderl5 import eccentricity_funcs_trunc14 as eccentricity_funcs_l5_trunc14
from .orderl5 import eccentricity_funcs_trunc16 as eccentricity_funcs_l5_trunc16
from .orderl5 import eccentricity_funcs_trunc18 as eccentricity_funcs_l5_trunc18
from .orderl5 import eccentricity_funcs_trunc20 as eccentricity_funcs_l5_trunc20

from .orderl6 import eccentricity_funcs_trunc2 as eccentricity_funcs_l6_trunc2
from .orderl6 import eccentricity_funcs_trunc4 as eccentricity_funcs_l6_trunc4
from .orderl6 import eccentricity_funcs_trunc6 as eccentricity_funcs_l6_trunc6
from .orderl6 import eccentricity_funcs_trunc8 as eccentricity_funcs_l6_trunc8
from .orderl6 import eccentricity_funcs_trunc10 as eccentricity_funcs_l6_trunc10
from .orderl6 import eccentricity_funcs_trunc12 as eccentricity_funcs_l6_trunc12
from .orderl6 import eccentricity_funcs_trunc14 as eccentricity_funcs_l6_trunc14
from .orderl6 import eccentricity_funcs_trunc16 as eccentricity_funcs_l6_trunc16
from .orderl6 import eccentricity_funcs_trunc18 as eccentricity_funcs_l6_trunc18
from .orderl6 import eccentricity_funcs_trunc20 as eccentricity_funcs_l6_trunc20

from .orderl7 import eccentricity_funcs_trunc2 as eccentricity_funcs_l7_trunc2
from .orderl7 import eccentricity_funcs_trunc4 as eccentricity_funcs_l7_trunc4
from .orderl7 import eccentricity_funcs_trunc6 as eccentricity_funcs_l7_trunc6
from .orderl7 import eccentricity_funcs_trunc8 as eccentricity_funcs_l7_trunc8
from .orderl7 import eccentricity_funcs_trunc10 as eccentricity_funcs_l7_trunc10
from .orderl7 import eccentricity_funcs_trunc12 as eccentricity_funcs_l7_trunc12
from .orderl7 import eccentricity_funcs_trunc14 as eccentricity_funcs_l7_trunc14
from .orderl7 import eccentricity_funcs_trunc16 as eccentricity_funcs_l7_trunc16
from .orderl7 import eccentricity_funcs_trunc18 as eccentricity_funcs_l7_trunc18
from .orderl7 import eccentricity_funcs_trunc20 as eccentricity_funcs_l7_trunc20

eccentricity_truncations = {
    2 :  # Truncation Level 2
        {
            2: eccentricity_funcs_l2_trunc2,
            3: eccentricity_funcs_l3_trunc2,
            4: eccentricity_funcs_l4_trunc2,
            5: eccentricity_funcs_l5_trunc2,
            6: eccentricity_funcs_l6_trunc2,
            7: eccentricity_funcs_l7_trunc2
            },
    4 :  # Truncation Level 4
        {
            2: eccentricity_funcs_l2_trunc4,
            3: eccentricity_funcs_l3_trunc4,
            4: eccentricity_funcs_l4_trunc4,
            5: eccentricity_funcs_l5_trunc4,
            6: eccentricity_funcs_l6_trunc4,
            7: eccentricity_funcs_l7_trunc4
            },
    6 :  # Truncation Level 6
        {
            2: eccentricity_funcs_l2_trunc6,
            3: eccentricity_funcs_l3_trunc6,
            4: eccentricity_funcs_l4_trunc6,
            5: eccentricity_funcs_l5_trunc6,
            6: eccentricity_funcs_l6_trunc6,
            7: eccentricity_funcs_l7_trunc6
            },
    8 :  # Truncation Level 8
        {
            2: eccentricity_funcs_l2_trunc8,
            3: eccentricity_funcs_l3_trunc8,
            4: eccentricity_funcs_l4_trunc8,
            5: eccentricity_funcs_l5_trunc8,
            6: eccentricity_funcs_l6_trunc8,
            7: eccentricity_funcs_l7_trunc8
            },
    10:  # Truncation Level 10
        {
            2: eccentricity_funcs_l2_trunc10,
            3: eccentricity_funcs_l3_trunc10,
            4: eccentricity_funcs_l4_trunc10,
            5: eccentricity_funcs_l5_trunc10,
            6: eccentricity_funcs_l6_trunc10,
            7: eccentricity_funcs_l7_trunc10
            },
    12:  # Truncation Level 12
        {
            2: eccentricity_funcs_l2_trunc12,
            3: eccentricity_funcs_l3_trunc12,
            4: eccentricity_funcs_l4_trunc12,
            5: eccentricity_funcs_l5_trunc12,
            6: eccentricity_funcs_l6_trunc12,
            7: eccentricity_funcs_l7_trunc12
            },
    14:  # Truncation Level 14
        {
            2: eccentricity_funcs_l2_trunc14,
            3: eccentricity_funcs_l3_trunc14,
            4: eccentricity_funcs_l4_trunc14,
            5: eccentricity_funcs_l5_trunc14,
            6: eccentricity_funcs_l6_trunc14,
            7: eccentricity_funcs_l7_trunc14
            },
    16:  # Truncation Level 16
        {
            2: eccentricity_funcs_l2_trunc16,
            3: eccentricity_funcs_l3_trunc16,
            4: eccentricity_funcs_l4_trunc16,
            5: eccentricity_funcs_l5_trunc16,
            6: eccentricity_funcs_l6_trunc16,
            7: eccentricity_funcs_l7_trunc16
            },
    18:  # Truncation Level 18
        {
            2: eccentricity_funcs_l2_trunc18,
            3: eccentricity_funcs_l3_trunc18,
            4: eccentricity_funcs_l4_trunc18,
            5: eccentricity_funcs_l5_trunc18,
            6: eccentricity_funcs_l6_trunc18,
            7: eccentricity_funcs_l7_trunc18
            },
    20:  # Truncation Level 20
        {
            2: eccentricity_funcs_l2_trunc20,
            3: eccentricity_funcs_l3_trunc20,
            4: eccentricity_funcs_l4_trunc20,
            5: eccentricity_funcs_l5_trunc20,
            6: eccentricity_funcs_l6_trunc20,
            7: eccentricity_funcs_l7_trunc20
            },
    22:  # Truncation Level 22
    # FIXME: This is not right! Only the l=2 has been implemented so far but njit won't compile the dict unless all
    #    these dicts have the same signature. Putting them all for l=2 for now...
        {
            2: eccentricity_funcs_l2_trunc22,
            3: eccentricity_funcs_l2_trunc22,
            4: eccentricity_funcs_l2_trunc22,
            5: eccentricity_funcs_l2_trunc22,
            6: eccentricity_funcs_l2_trunc22,
            7: eccentricity_funcs_l2_trunc22
            }
    }
