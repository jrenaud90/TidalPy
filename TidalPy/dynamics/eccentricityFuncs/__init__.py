from .orderl2 import eccentricity_funcs_trunc2 as eccentricity_funcs_l2_trunc2
from .orderl2 import eccentricity_funcs_trunc6 as eccentricity_funcs_l2_trunc6
from .orderl2 import eccentricity_funcs_trunc10 as eccentricity_funcs_l2_trunc10
from .orderl2 import eccentricity_funcs_trunc18 as eccentricity_funcs_l2_trunc12
from .orderl2 import eccentricity_funcs_trunc18 as eccentricity_funcs_l2_trunc14
from .orderl2 import eccentricity_funcs_trunc18 as eccentricity_funcs_l2_trunc16
from .orderl2 import eccentricity_funcs_trunc18 as eccentricity_funcs_l2_trunc18
from .orderl2 import eccentricity_funcs_trunc20 as eccentricity_funcs_l2_trunc20

from .orderl3 import eccentricity_funcs_trunc2 as eccentricity_funcs_l3_trunc2
from .orderl3 import eccentricity_funcs_trunc10 as eccentricity_funcs_l3_trunc10
from .orderl3 import eccentricity_funcs_trunc12 as eccentricity_funcs_l3_trunc12
from .orderl3 import eccentricity_funcs_trunc14 as eccentricity_funcs_l3_trunc14

from .orderl4 import eccentricity_funcs_trunc2 as eccentricity_funcs_l4_trunc2
from .orderl4 import eccentricity_funcs_trunc10 as eccentricity_funcs_l4_trunc10
from .orderl4 import eccentricity_funcs_trunc12 as eccentricity_funcs_l4_trunc12
from .orderl4 import eccentricity_funcs_trunc14 as eccentricity_funcs_l4_trunc14

from .orderl5 import eccentricity_funcs_trunc2 as eccentricity_funcs_l5_trunc2
from .orderl5 import eccentricity_funcs_trunc10 as eccentricity_funcs_l5_trunc10
from .orderl5 import eccentricity_funcs_trunc12 as eccentricity_funcs_l5_trunc12
from .orderl5 import eccentricity_funcs_trunc14 as eccentricity_funcs_l5_trunc14


# Change this from a dict to a tuple.
eccentricity_truncations = {
    # Truncation Level 2
    2: (
        eccentricity_funcs_l2_trunc2,
        eccentricity_funcs_l3_trunc2,
        eccentricity_funcs_l4_trunc2,
        eccentricity_funcs_l5_trunc2
    ),
    # # Truncation Level 4
    # (
    #     eccentricity_funcs_l2_trunc4,
    #     eccentricity_funcs_l3_trunc4,
    # ),
    # Truncation Level 6
    6: (
        eccentricity_funcs_l2_trunc6,
        # eccentricity_funcs_l3_trunc6,
    ),
    # # Truncation Level 8
    # (
    #     eccentricity_funcs_l2_trunc8,
    #     eccentricity_funcs_l3_trunc8,
    # ),
    # Truncation Level 10
    10: (
        eccentricity_funcs_l2_trunc10,
        eccentricity_funcs_l3_trunc10,
        eccentricity_funcs_l4_trunc10,
        eccentricity_funcs_l5_trunc10
    ),
    # Truncation Level 12
    12 : (
        eccentricity_funcs_l2_trunc12,
        eccentricity_funcs_l3_trunc12,
        eccentricity_funcs_l4_trunc12,
        eccentricity_funcs_l5_trunc12
    ),
    # Truncation Level 14
    14 : (
        eccentricity_funcs_l2_trunc14,
        eccentricity_funcs_l3_trunc14,
        eccentricity_funcs_l4_trunc14,
        eccentricity_funcs_l5_trunc14
    ),
    # # Truncation Level 16
    # (
    #     eccentricity_funcs_l2_trunc16,
    #     eccentricity_funcs_l3_trunc16,
    # ),
    # Truncation Level 18
    18: (
        eccentricity_funcs_l2_trunc18,
        # eccentricity_funcs_l3_trunc18,
    ),
    # Truncation Level 20
    20: (
        eccentricity_funcs_l2_trunc20,
        # eccentricity_funcs_l3_trunc20,
    ),
}