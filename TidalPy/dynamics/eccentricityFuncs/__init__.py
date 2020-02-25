from .orderl2 import eccentricity_funcs_trunc6 as eccentricity_funcs_l2_trunc6
from .orderl2 import eccentricity_funcs_trunc10 as eccentricity_funcs_l2_trunc10
from .orderl2 import eccentricity_funcs_trunc18 as eccentricity_funcs_l2_trunc18

# Change this from a dict to a tuple.
eccentricity_truncations = {
    # # Truncation Level 2
    # (
    #     eccentricity_funcs_l2_trunc2,
    #     eccentricity_funcs_l3_trunc2,
    # ),
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
        # eccentricity_funcs_l3_trunc10,
    ),
    # # Truncation Level 12
    # (
    #     eccentricity_funcs_l2_trunc12,
    #     eccentricity_funcs_l3_trunc12,
    # ),
    # # Truncation Level 14
    # (
    #     eccentricity_funcs_l2_trunc14,
    #     eccentricity_funcs_l3_trunc14,
    # ),
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
    # # Truncation Level 20
    # (
    #     eccentricity_funcs_l2_trunc20,
    #     eccentricity_funcs_l3_trunc20,
    # ),
}