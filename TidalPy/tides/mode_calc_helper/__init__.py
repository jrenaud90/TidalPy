from . import eccen_calc_orderl2, eccen_calc_orderl3, eccen_calc_orderl4, eccen_calc_orderl5, eccen_calc_orderl6, \
    eccen_calc_orderl7, inclin_calc_orderl2, inclin_calc_orderl3, inclin_calc_orderl4, inclin_calc_orderl5, \
    inclin_calc_orderl6, inclin_calc_orderl7

eccentricity_functions_lookup = {
    # Eccentricity functions are stored by desired eccentricity truncation level and tidal order l.
    2: {
        2: eccen_calc_orderl2.eccentricity_truncation_2_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_2_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_2_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_2_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_2_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_2_maxl_7
    },
    4: {
        2: eccen_calc_orderl2.eccentricity_truncation_4_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_4_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_4_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_4_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_4_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_4_maxl_7
    },
    6: {
        2: eccen_calc_orderl2.eccentricity_truncation_6_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_6_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_6_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_6_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_6_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_6_maxl_7
    },
    8: {
        2: eccen_calc_orderl2.eccentricity_truncation_8_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_8_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_8_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_8_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_8_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_8_maxl_7
    },
    10: {
        2: eccen_calc_orderl2.eccentricity_truncation_10_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_10_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_10_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_10_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_10_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_10_maxl_7
    },
    12: {
        2: eccen_calc_orderl2.eccentricity_truncation_12_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_12_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_12_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_12_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_12_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_12_maxl_7
    },
    14: {
        2: eccen_calc_orderl2.eccentricity_truncation_14_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_14_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_14_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_14_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_14_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_14_maxl_7
    },
    16: {
        2: eccen_calc_orderl2.eccentricity_truncation_16_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_16_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_16_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_16_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_16_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_16_maxl_7
    },
    18: {
        2: eccen_calc_orderl2.eccentricity_truncation_18_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_18_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_18_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_18_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_18_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_18_maxl_7
    },
    20: {
        2: eccen_calc_orderl2.eccentricity_truncation_20_maxl_2,
        3: eccen_calc_orderl3.eccentricity_truncation_20_maxl_3,
        4: eccen_calc_orderl4.eccentricity_truncation_20_maxl_4,
        5: eccen_calc_orderl5.eccentricity_truncation_20_maxl_5,
        6: eccen_calc_orderl6.eccentricity_truncation_20_maxl_6,
        7: eccen_calc_orderl7.eccentricity_truncation_20_maxl_7
    },
    22: {
        2: eccen_calc_orderl2.eccentricity_truncation_22_maxl_2
    }
}

inclination_functions_lookup = {
    # Obliquity functions are stored by if the obliquity is always zero or not and tidal order l.
    True: {
        2: inclin_calc_orderl2.inclination_on_maxl_2,
        3: inclin_calc_orderl3.inclination_on_maxl_3,
        4: inclin_calc_orderl4.inclination_on_maxl_4,
        5: inclin_calc_orderl5.inclination_on_maxl_5,
        6: inclin_calc_orderl6.inclination_on_maxl_6,
        7: inclin_calc_orderl7.inclination_on_maxl_7
    },
    False: {
        2: inclin_calc_orderl2.inclination_off_maxl_2,
        3: inclin_calc_orderl3.inclination_off_maxl_3,
        4: inclin_calc_orderl4.inclination_off_maxl_4,
        5: inclin_calc_orderl5.inclination_off_maxl_5,
        6: inclin_calc_orderl6.inclination_off_maxl_6,
        7: inclin_calc_orderl7.inclination_off_maxl_7
    }
}