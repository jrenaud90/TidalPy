""" Eccentricity functions (squared) for various truncations of e at tidal order-l = 7
"""

from typing import Tuple, Dict

import numpy as np

from ...types import FloatArray
from ...performance import njit


@njit
def eccentricity_funcs_trunc2(eccentricity: FloatArray) -> Dict[int, Dict[int, FloatArray]]:
    """ Calculates the eccentricity functions (by mode) truncated to e^2 for order-l = 7
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : Dict[int, Dict[int, FloatArray]]
    """
    # Eccentricity functions calculated at truncation level 2.
    #     and order-l = 7.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -1: 9.0*e2,
            0: 1.0 - 70.0*e2,
            1: 121.0*e2,
        },
        1: {
            -1: e2,
            0: 1.0 - 22.0*e2,
            1: 81.0*e2,
        },
        2: {
            -1: e2,
            0: 10.0*e2 + 1.0,
            1: 49.0*e2,
        },
        3: {
            -1: -0.140625*e2*(5.0*e**4 + 20.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            0: 26.0*e2 + 1.0,
            1: 25.0*e2,
        },
        4: {
            -1: 25.0*e2,
        },
        5: {
            -1: 49.0*e2,
        },
        6: {
            -1: 81.0*e2,
        },
        7: {
            -1: 121.0*e2,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[3][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[7][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[7][1] = eccentricity_results_bymode[0][-1]

    return eccentricity_results_bymode

@njit
def eccentricity_funcs_trunc4(eccentricity: FloatArray) -> Dict[int, Dict[int, FloatArray]]:
    """ Calculates the eccentricity functions (by mode) truncated to e^4 for order-l = 7
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : Dict[int, Dict[int, FloatArray]]
    """
    # Eccentricity functions calculated at truncation level 4.
    #     and order-l = 7.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -2: 9.765625*e4,
            -1: -238.5*e4 + 9.0*e2,
            0: 1785.21875*e4 - 70.0*e2 + 1.0,
            1: -5016.0*e4 + 121.0*e2,
            2: 4607.015625*e4,
        },
        1: {
            -2: 0.140625*e4,
            -1: -9.0*e4 + e2,
            0: 183.96875*e4 - 22.0*e2 + 1.0,
            1: -1201.5*e4 + 81.0*e2,
            2: 2173.890625*e4,
        },
        2: {
            -3: -0.09765625*e**6*(3.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            -2: 2.640625*e4,
            -1: 16.5*e4 + e2,
            0: 76.46875*e4 + 10.0*e2 + 1.0,
            1: 189.0*e4 + 49.0*e2,
            2: 862.890625*e4,
        },
        3: {
            -2: 47.265625*e4,
            -1: -0.140625*e2*(5.0*e4 + 20.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            0: 310.71875*e4 + 26.0*e2 + 1.0,
            1: 367.5*e4 + 25.0*e2,
            2: 260.015625*e4,
        },
        4: {
            -2: 260.015625*e4,
        },
        5: {
            -2: 862.890625*e4,
        },
        6: {
            -2: 2173.890625*e4,
        },
        7: {
            -2: 4607.015625*e4,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[3][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[3][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[7][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[7][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[7][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[7][2] = eccentricity_results_bymode[0][-2]

    return eccentricity_results_bymode

@njit
def eccentricity_funcs_trunc6(eccentricity: FloatArray) -> Dict[int, Dict[int, FloatArray]]:
    """ Calculates the eccentricity functions (by mode) truncated to e^6 for order-l = 7
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : Dict[int, Dict[int, FloatArray]]
    """
    # Eccentricity functions calculated at truncation level 6.
    #     and order-l = 7.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e
    e6 = e*e*e*e*e*e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -3: 1.77777777777778*e6,
            -2: -120.442708333333*e6 + 9.765625*e4,
            -1: 2496.9375*e6 - 238.5*e4 + 9.0*e2,
            0: -21364.4618055556*e6 + 1785.21875*e4 - 70.0*e2 + 1.0,
            1: 83775.8333333333*e6 - 5016.0*e4 + 121.0*e2,
            2: -148518.984375*e6 + 4607.015625*e4,
            3: 95841.8402777778*e6,
        },
        1: {
            -3: 0.00694444444444444*e6,
            -2: 0.140625*e6 + 0.140625*e4,
            -1: 24.5833333333333*e6 - 9.0*e4 + e2,
            0: -741.545138888889*e6 + 183.96875*e4 - 22.0*e2 + 1.0,
            1: 7388.4375*e6 - 1201.5*e4 + 81.0*e2,
            2: -27743.8177083333*e6 + 2173.890625*e4,
            3: 33184.6944444444*e6,
        },
        2: {
            -3: -0.09765625*e6*(3.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            -2: 38.9322916666667*e6 + 2.640625*e4,
            -1: 139.270833333333*e6 + 16.5*e4 + e2,
            0: 424.84375*e6 + 76.46875*e4 + 10.0*e2 + 1.0,
            1: 1080.58333333333*e6 + 189.0*e4 + 49.0*e2,
            2: 1254.55729166667*e6 + 862.890625*e4,
            3: 9168.0625*e6,
        },
        3: {
            -3: 193.673611111111*e6,
            -2: 688.932291666667*e6 + 47.265625*e4,
            -1: -0.140625*e2*(5.0*e4 + 20.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            0: 2353.37152777778*e6 + 310.71875*e4 + 26.0*e2 + 1.0,
            1: 2885.77083333333*e6 + 367.5*e4 + 25.0*e2,
            2: 2727.140625*e6 + 260.015625*e4,
            3: 1792.11111111111*e6,
        },
        4: {
            -3: 1792.11111111111*e6,
        },
        5: {
            -3: 9168.0625*e6,
        },
        6: {
            -3: 33184.6944444444*e6,
        },
        7: {
            -3: 95841.8402777778*e6,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[3][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[3][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[3][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[7][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[7][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[7][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[7][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[7][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[7][3] = eccentricity_results_bymode[0][-3]

    return eccentricity_results_bymode

@njit
def eccentricity_funcs_trunc8(eccentricity: FloatArray) -> Dict[int, Dict[int, FloatArray]]:
    """ Calculates the eccentricity functions (by mode) truncated to e^8 for order-l = 7
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : Dict[int, Dict[int, FloatArray]]
    """
    # Eccentricity functions calculated at truncation level 8.
    #     and order-l = 7.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e
    e6 = e*e*e*e*e*e
    e8 = e**8

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -4: 0.04449462890625*e8,
            -3: -9.77777777777778*e8 + 1.77777777777778*e6,
            -2: 591.854519314236*e8 - 120.442708333333*e6 + 9.765625*e4,
            -1: -13602.09375*e8 + 2496.9375*e6 - 238.5*e4 + 9.0*e2,
            0: 142589.589938694*e8 - 21364.4618055556*e6 + 1785.21875*e4 - 70.0*e2 + 1.0,
            1: -748898.888888889*e8 + 83775.8333333333*e6 - 5016.0*e4 + 121.0*e2,
            2: 2024000.86157227*e8 - 148518.984375*e6 + 4607.015625*e4,
            3: -2669833.76736111*e8 + 95841.8402777778*e6,
            4: 1352599.28662788*e8,
        },
        1: {
            -5: -0.03515625*e**10/(e2 - 1.0)**13,
            -4: 0.0190497504340278*e8,
            -3: 0.118055555555556*e8 + 0.00694444444444444*e6,
            -2: 1.742431640625*e8 + 0.140625*e6 + 0.140625*e4,
            -1: -32.8888888888889*e8 + 24.5833333333333*e6 - 9.0*e4 + e2,
            0: 1576.83668348524*e8 - 741.545138888889*e6 + 183.96875*e4 - 22.0*e2 + 1.0,
            1: -24466.5*e8 + 7388.4375*e6 - 1201.5*e4 + 81.0*e2,
            2: 152617.184597439*e8 - 27743.8177083333*e6 + 2173.890625*e4,
            3: -397130.923611111*e8 + 33184.6944444444*e6,
            4: 354871.521057129*e8,
        },
        2: {
            -4: 14.396491156684*e8,
            -3: -0.09765625*e6*(3.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            -2: 303.058295355903*e8 + 38.9322916666667*e6 + 2.640625*e4,
            -1: 811.631944444444*e8 + 139.270833333333*e6 + 16.5*e4 + e2,
            0: 1945.65222167969*e8 + 424.84375*e6 + 76.46875*e4 + 10.0*e2 + 1.0,
            1: 4115.31944444444*e8 + 1080.58333333333*e6 + 189.0*e4 + 49.0*e2,
            2: 8894.1413031684*e8 + 1254.55729166667*e6 + 862.890625*e4,
            3: -95.75*e8 + 9168.0625*e6,
            4: 71362.7221747504*e8,
        },
        3: {
            -4: 685.376037597656*e8,
            -3: 2413.96180555556*e8 + 193.673611111111*e6,
            -2: 5301.66408962674*e8 + 688.932291666667*e6 + 47.265625*e4,
            -1: -0.140625*e2*(5.0*e4 + 20.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            0: 13103.2553032769*e8 + 2353.37152777778*e6 + 310.71875*e4 + 26.0*e2 + 1.0,
            1: 16012.0798611111*e8 + 2885.77083333333*e6 + 367.5*e4 + 25.0*e2,
            2: 16688.9113769531*e8 + 2727.140625*e6 + 260.015625*e4,
            3: 14474.4722222222*e8 + 1792.11111111111*e6,
            4: 9661.76367865668*e8,
        },
        4: {
            -4: 9661.76367865668*e8,
        },
        5: {
            -4: 71362.7221747504*e8,
        },
        6: {
            -4: 354871.521057129*e8,
        },
        7: {
            -4: 1352599.28662788*e8,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[3][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[3][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[3][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[3][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[7][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[7][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[7][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[7][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[7][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[7][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[7][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[7][4] = eccentricity_results_bymode[0][-4]

    return eccentricity_results_bymode

@njit
def eccentricity_funcs_trunc10(eccentricity: FloatArray) -> Dict[int, Dict[int, FloatArray]]:
    """ Calculates the eccentricity functions (by mode) truncated to e^10 for order-l = 7
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : Dict[int, Dict[int, FloatArray]]
    """
    # Eccentricity functions calculated at truncation level 10.
    #     and order-l = 7.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e
    e6 = e*e*e*e*e*e
    e8 = e**8
    e10 = e**10

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -5: 6.94444444444444e-5*e10,
            -4: -0.07119140625*e10 + 0.04449462890625*e8,
            -3: 20.1555555555556*e10 - 9.77777777777778*e8 + 1.77777777777778*e6,
            -2: -1516.85926649306*e10 + 591.854519314236*e8 - 120.442708333333*e6 + 9.765625*e4,
            -1: 43713.97734375*e10 - 13602.09375*e8 + 2496.9375*e6 - 238.5*e4 + 9.0*e2,
            0: -586640.268715278*e10 + 142589.589938694*e8 - 21364.4618055556*e6 + 1785.21875*e4 - 70.0*e2 + 1.0,
            1: 4085875.57048611*e10 - 748898.888888889*e8 + 83775.8333333333*e6 - 5016.0*e4 + 121.0*e2,
            2: -15548508.0931641*e10 + 2024000.86157227*e8 - 148518.984375*e6 + 4607.015625*e4,
            3: 32361560.6011285*e10 - 2669833.76736111*e8 + 95841.8402777778*e6,
            4: -34373096.9709961*e10 + 1352599.28662788*e8,
            5: 14501244.8025*e10,
        },
        1: {
            -5: -0.03515625*e10/(e2 - 1.0)**13,
            -4: 0.257063802083333*e10 + 0.0190497504340278*e8,
            -3: 0.971701388888889*e10 + 0.118055555555556*e8 + 0.00694444444444444*e6,
            -2: 5.8693359375*e10 + 1.742431640625*e8 + 0.140625*e6 + 0.140625*e4,
            -1: 38.0814236111111*e10 - 32.8888888888889*e8 + 24.5833333333333*e6 - 9.0*e4 + e2,
            0: -2021.14040798611*e10 + 1576.83668348524*e8 - 741.545138888889*e6 + 183.96875*e4 - 22.0*e2 + 1.0,
            1: 48573.09140625*e10 - 24466.5*e8 + 7388.4375*e6 - 1201.5*e4 + 81.0*e2,
            2: -475122.020203993*e10 + 152617.184597439*e8 - 27743.8177083333*e6 + 2173.890625*e4,
            3: 2097863.73862847*e10 - 397130.923611111*e8 + 33184.6944444444*e6,
            4: -4146304.49912109*e10 + 354871.521057129*e8,
            5: 2959618.45876736*e10,
        },
        2: {
            -5: 32.4425173611111*e10,
            -4: 185.11298828125*e10 + 14.396491156684*e8,
            -3: -0.09765625*e6*(3.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            -2: 1659.23604600694*e10 + 303.058295355903*e8 + 38.9322916666667*e6 + 2.640625*e4,
            -1: 3694.66293402778*e10 + 811.631944444444*e8 + 139.270833333333*e6 + 16.5*e4 + e2,
            0: 7565.199296875*e10 + 1945.65222167969*e8 + 424.84375*e6 + 76.46875*e4 + 10.0*e2 + 1.0,
            1: 14360.1105902778*e10 + 4115.31944444444*e8 + 1080.58333333333*e6 + 189.0*e4 + 49.0*e2,
            2: 24376.8419053819*e10 + 8894.1413031684*e8 + 1254.55729166667*e6 + 862.890625*e4,
            3: 59750.0453125*e10 - 95.75*e8 + 9168.0625*e6,
            4: -78294.8146809896*e10 + 71362.7221747504*e8,
            5: 448693.440434028*e10,
        },
        3: {
            -5: 2198.82840277778*e10,
            -4: 7485.75439453125*e10 + 685.376037597656*e8,
            -3: 16387.2709201389*e10 + 2413.96180555556*e8 + 193.673611111111*e6,
            -2: 28704.1544053819*e10 + 5301.66408962674*e8 + 688.932291666667*e6 + 47.265625*e4,
            -1: -0.140625*e2*(5.0*e4 + 20.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            0: 58347.2753993056*e10 + 13103.2553032769*e8 + 2353.37152777778*e6 + 310.71875*e4 + 26.0*e2 + 1.0,
            1: 70311.6384548611*e10 + 16012.0798611111*e8 + 2885.77083333333*e6 + 367.5*e4 + 25.0*e2,
            2: 76178.1462890625*e10 + 16688.9113769531*e8 + 2727.140625*e6 + 260.015625*e4,
            3: 73701.4461805556*e10 + 14474.4722222222*e8 + 1792.11111111111*e6,
            4: 61808.6661783854*e10 + 9661.76367865668*e8,
            5: 44163.0225*e10,
        },
        4: {
            -5: 44163.0225*e10,
        },
        5: {
            -5: 448693.440434028*e10,
        },
        6: {
            -5: 2959618.45876736*e10,
        },
        7: {
            -5: 14501244.8025*e10,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[3][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[3][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[3][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[3][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[3][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[7][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[7][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[7][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[7][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[7][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[7][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[7][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[7][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[7][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[7][5] = eccentricity_results_bymode[0][-5]

    return eccentricity_results_bymode

@njit
def eccentricity_funcs_trunc12(eccentricity: FloatArray) -> Dict[int, Dict[int, FloatArray]]:
    """ Calculates the eccentricity functions (by mode) truncated to e^12 for order-l = 7
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : Dict[int, Dict[int, FloatArray]]
    """
    # Eccentricity functions calculated at truncation level 12.
    #     and order-l = 7.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e
    e6 = e*e*e*e*e*e
    e8 = e**8
    e10 = e**10
    e12 = e**12

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -6: 4.7095027970679e-10*e12,
            -5: 5.78703703703704e-5*e12 + 6.94444444444444e-5*e10,
            -4: 0.0320917510986328*e12 - 0.07119140625*e10 + 0.04449462890625*e8,
            -3: -20.1716049382716*e12 + 20.1555555555556*e10 - 9.77777777777778*e8 + 1.77777777777778*e6,
            -2: 2269.74851555294*e12 - 1516.85926649306*e10 + 591.854519314236*e8 - 120.442708333333*e6 + 9.765625*e4,
            -1: -89131.943671875*e12 + 43713.97734375*e10 - 13602.09375*e8 + 2496.9375*e6 - 238.5*e4 + 9.0*e2,
            0: 1586997.38291086*e12 - 586640.268715278*e10 + 142589.589938694*e8 - 21364.4618055556*e6 + 1785.21875*e4 - 70.0*e2 + 1.0,
            1: -14695176.0341667*e12 + 4085875.57048611*e10 - 748898.888888889*e8 + 83775.8333333333*e6 - 5016.0*e4 + 121.0*e2,
            2: 76176310.2957281*e12 - 15548508.0931641*e10 + 2024000.86157227*e8 - 148518.984375*e6 + 4607.015625*e4,
            3: -227752424.952377*e12 + 32361560.6011285*e10 - 2669833.76736111*e8 + 95841.8402777778*e6,
            4: 387648466.517466*e12 - 34373096.9709961*e10 + 1352599.28662788*e8,
            5: -347139933.975*e12 + 14501244.8025*e10,
            6: 126441284.408644*e12,
        },
        1: {
            -6: 0.0649431247475706*e12,
            -5: -0.03515625*e10/(e2 - 1.0)**13,
            -4: 1.8597327875208*e12 + 0.257063802083333*e10 + 0.0190497504340278*e8,
            -3: 5.39066358024691*e12 + 0.971701388888889*e10 + 0.118055555555556*e8 + 0.00694444444444444*e6,
            -2: 20.3473989486694*e12 + 5.8693359375*e10 + 1.742431640625*e8 + 0.140625*e6 + 0.140625*e4,
            -1: 30.2398958333333*e12 + 38.0814236111111*e10 - 32.8888888888889*e8 + 24.5833333333333*e6 - 9.0*e4 + e2,
            0: 1725.84582929258*e12 - 2021.14040798611*e10 + 1576.83668348524*e8 - 741.545138888889*e6 + 183.96875*e4 - 22.0*e2 + 1.0,
            1: -62900.057109375*e12 + 48573.09140625*e10 - 24466.5*e8 + 7388.4375*e6 - 1201.5*e4 + 81.0*e2,
            2: 938547.352011341*e12 - 475122.020203993*e10 + 152617.184597439*e8 - 27743.8177083333*e6 + 2173.890625*e4,
            3: -6464235.02122878*e12 + 2097863.73862847*e10 - 397130.923611111*e8 + 33184.6944444444*e6,
            4: 21747975.9750501*e12 - 4146304.49912109*e10 + 354871.521057129*e8,
            5: -34500980.0248119*e12 + 2959618.45876736*e10,
            6: 20518448.8791875*e12,
        },
        2: {
            -6: 71.8156937026978*e12,
            -5: 390.619458912037*e12 + 32.4425173611111*e10,
            -4: 1283.45398924368*e12 + 185.11298828125*e10 + 14.396491156684*e8,
            -3: -0.09765625*e6*(3.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            -2: 7178.43844695621*e12 + 1659.23604600694*e10 + 303.058295355903*e8 + 38.9322916666667*e6 + 2.640625*e4,
            -1: 14041.6768142361*e12 + 3694.66293402778*e10 + 811.631944444444*e8 + 139.270833333333*e6 + 16.5*e4 + e2,
            0: 25731.1676902771*e12 + 7565.199296875*e10 + 1945.65222167969*e8 + 424.84375*e6 + 76.46875*e4 + 10.0*e2 + 1.0,
            1: 44644.5701736111*e12 + 14360.1105902778*e10 + 4115.31944444444*e8 + 1080.58333333333*e6 + 189.0*e4 + 49.0*e2,
            2: 73906.7868524128*e12 + 24376.8419053819*e10 + 8894.1413031684*e8 + 1254.55729166667*e6 + 862.890625*e4,
            3: 97963.396875*e12 + 59750.0453125*e10 - 95.75*e8 + 9168.0625*e6,
            4: 392310.672480862*e12 - 78294.8146809896*e10 + 71362.7221747504*e8,
            5: -894602.368952546*e12 + 448693.440434028*e10,
            6: 2412439.34078983*e12,
        },
        3: {
            -6: 6572.37797547564*e12,
            -5: 21292.203287037*e12 + 2198.82840277778*e10,
            -4: 45732.4679340363*e12 + 7485.75439453125*e10 + 685.376037597656*e8,
            -3: 80064.9051456404*e12 + 16387.2709201389*e10 + 2413.96180555556*e8 + 193.673611111111*e6,
            -2: 122862.972033925*e12 + 28704.1544053819*e10 + 5301.66408962674*e8 + 688.932291666667*e6 + 47.265625*e4,
            -1: -0.140625*e2*(5.0*e4 + 20.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            0: 219108.287543295*e12 + 58347.2753993056*e10 + 13103.2553032769*e8 + 2353.37152777778*e6 + 310.71875*e4 + 26.0*e2 + 1.0,
            1: 259842.755598958*e12 + 70311.6384548611*e10 + 16012.0798611111*e8 + 2885.77083333333*e6 + 367.5*e4 + 25.0*e2,
            2: 285788.216771507*e12 + 76178.1462890625*e10 + 16688.9113769531*e8 + 2727.140625*e6 + 260.015625*e4,
            3: 290946.807098765*e12 + 73701.4461805556*e10 + 14474.4722222222*e8 + 1792.11111111111*e6,
            4: 271561.904476236*e12 + 61808.6661783854*e10 + 9661.76367865668*e8,
            5: 224786.9475*e12 + 44163.0225*e10,
            6: 179171.885103873*e12,
        },
        4: {
            -6: 179171.885103873*e12,
        },
        5: {
            -6: 2412439.34078983*e12,
        },
        6: {
            -6: 20518448.8791875*e12,
        },
        7: {
            -6: 126441284.408644*e12,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[3][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[3][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[3][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[3][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[3][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[3][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[3][-6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[6][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[6][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[7][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[7][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[7][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[7][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[7][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[7][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[7][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[7][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[7][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[7][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[7][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[7][6] = eccentricity_results_bymode[0][-6]

    return eccentricity_results_bymode

@njit
def eccentricity_funcs_trunc14(eccentricity: FloatArray) -> Dict[int, Dict[int, FloatArray]]:
    """ Calculates the eccentricity functions (by mode) truncated to e^14 for order-l = 7
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : Dict[int, Dict[int, FloatArray]]
    """
    # Eccentricity functions calculated at truncation level 14.
    #     and order-l = 7.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e
    e6 = e*e*e*e*e*e
    e8 = e**8
    e10 = e**10
    e12 = e**12
    e14 = e**14

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -6: 1.14373639357363e-9*e14 + 4.7095027970679e-10*e12,
            -5: 6.28995811287478e-5*e14 + 5.78703703703704e-5*e12 + 6.94444444444444e-5*e10,
            -4: -0.0112746211460658*e14 + 0.0320917510986328*e12 - 0.07119140625*e10 + 0.04449462890625*e8,
            -3: 11.5056128747795*e14 - 20.1716049382716*e12 + 20.1555555555556*e10 - 9.77777777777778*e8 + 1.77777777777778*e6,
            -2: -2132.90065789769*e14 + 2269.74851555294*e12 - 1516.85926649306*e10 + 591.854519314236*e8 - 120.442708333333*e6 + 9.765625*e4,
            -1: 121302.196799665*e14 - 89131.943671875*e12 + 43713.97734375*e10 - 13602.09375*e8 + 2496.9375*e6 - 238.5*e4 + 9.0*e2,
            0: -2955172.8970563*e14 + 1586997.38291086*e12 - 586640.268715278*e10 + 142589.589938694*e8 - 21364.4618055556*e6 + 1785.21875*e4 - 70.0*e2 + 1.0,
            1: 36681269.0568395*e14 - 14695176.0341667*e12 + 4085875.57048611*e10 - 748898.888888889*e8 + 83775.8333333333*e6 - 5016.0*e4 + 121.0*e2,
            2: -255474118.187811*e14 + 76176310.2957281*e12 - 15548508.0931641*e10 + 2024000.86157227*e8 - 148518.984375*e6 + 4607.015625*e4,
            3: 1050342781.37476*e14 - 227752424.952377*e12 + 32361560.6011285*e10 - 2669833.76736111*e8 + 95841.8402777778*e6,
            4: -2589179614.74784*e14 + 387648466.517466*e12 - 34373096.9709961*e10 + 1352599.28662788*e8,
            5: 3741770950.29402*e14 - 347139933.975*e12 + 14501244.8025*e10,
            6: -2910686193.59795*e14 + 126441284.408644*e12,
            7: 937784248.676216*e14,
        },
        1: {
            -7: 0.119978389648684*e14,
            -6: 0.811787874265533*e14 + 0.0649431247475706*e12,
            -5: -0.03515625*e10/(e2 - 1.0)**13,
            -4: 9.57492866166264*e14 + 1.8597327875208*e12 + 0.257063802083333*e10 + 0.0190497504340278*e8,
            -3: 23.2613752652392*e14 + 5.39066358024691*e12 + 0.971701388888889*e10 + 0.118055555555556*e8 + 0.00694444444444444*e6,
            -2: 64.667712789263*e14 + 20.3473989486694*e12 + 5.8693359375*e10 + 1.742431640625*e8 + 0.140625*e6 + 0.140625*e4,
            -1: 139.310944113757*e14 + 30.2398958333333*e12 + 38.0814236111111*e10 - 32.8888888888889*e8 + 24.5833333333333*e6 - 9.0*e4 + e2,
            0: -752.495595552844*e14 + 1725.84582929258*e12 - 2021.14040798611*e10 + 1576.83668348524*e8 - 741.545138888889*e6 + 183.96875*e4 - 22.0*e2 + 1.0,
            1: 56545.2157073103*e14 - 62900.057109375*e12 + 48573.09140625*e10 - 24466.5*e8 + 7388.4375*e6 - 1201.5*e4 + 81.0*e2,
            2: -1266309.79616931*e14 + 938547.352011341*e12 - 475122.020203993*e10 + 152617.184597439*e8 - 27743.8177083333*e6 + 2173.890625*e4,
            3: 13062709.1406972*e14 - 6464235.02122878*e12 + 2097863.73862847*e10 - 397130.923611111*e8 + 33184.6944444444*e6,
            4: -67933871.4747141*e14 + 21747975.9750501*e12 - 4146304.49912109*e10 + 354871.521057129*e8,
            5: 182888425.081913*e14 - 34500980.0248119*e12 + 2959618.45876736*e10,
            6: -241659566.56067*e14 + 20518448.8791875*e12,
            7: 123255910.704605*e14,
        },
        2: {
            -7: 156.64211077255*e14,
            -6: 809.854400596619*e14 + 71.8156937026978*e12,
            -5: 2563.02166482377*e14 + 390.619458912037*e12 + 32.4425173611111*e10,
            -4: 6363.25240061133*e14 + 1283.45398924368*e12 + 185.11298828125*e10 + 14.396491156684*e8,
            -3: -0.09765625*e6*(3.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            -2: 26135.9927566851*e14 + 7178.43844695621*e12 + 1659.23604600694*e10 + 303.058295355903*e8 + 38.9322916666667*e6 + 2.640625*e4,
            -1: 46429.6170135272*e14 + 14041.6768142361*e12 + 3694.66293402778*e10 + 811.631944444444*e8 + 139.270833333333*e6 + 16.5*e4 + e2,
            0: 78327.7874431159*e14 + 25731.1676902771*e12 + 7565.199296875*e10 + 1945.65222167969*e8 + 424.84375*e6 + 76.46875*e4 + 10.0*e2 + 1.0,
            1: 126743.549716435*e14 + 44644.5701736111*e12 + 14360.1105902778*e10 + 4115.31944444444*e8 + 1080.58333333333*e6 + 189.0*e4 + 49.0*e2,
            2: 197042.647860248*e14 + 73906.7868524128*e12 + 24376.8419053819*e10 + 8894.1413031684*e8 + 1254.55729166667*e6 + 862.890625*e4,
            3: 306304.572200056*e14 + 97963.396875*e12 + 59750.0453125*e10 - 95.75*e8 + 9168.0625*e6,
            4: 196218.790338876*e14 + 392310.672480862*e12 - 78294.8146809896*e10 + 71362.7221747504*e8,
            5: 2608659.94096175*e14 - 894602.368952546*e12 + 448693.440434028*e10,
            6: -6709222.00565983*e14 + 2412439.34078983*e12,
            7: 11503209.3055423*e14,
        },
        3: {
            -7: 18617.627869898*e14,
            -6: 56701.7701593513*e14 + 6572.37797547564*e12,
            -5: 118393.583558477*e14 + 21292.203287037*e12 + 2198.82840277778*e10,
            -4: 204879.526902172*e14 + 45732.4679340363*e12 + 7485.75439453125*e10 + 685.376037597656*e8,
            -3: 314588.587124842*e14 + 80064.9051456404*e12 + 16387.2709201389*e10 + 2413.96180555556*e8 + 193.673611111111*e6,
            -2: 442814.440752841*e14 + 122862.972033925*e12 + 28704.1544053819*e10 + 5301.66408962674*e8 + 688.932291666667*e6 + 47.265625*e4,
            -1: -0.140625*e2*(5.0*e4 + 20.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            0: 719186.088857198*e14 + 219108.287543295*e12 + 58347.2753993056*e10 + 13103.2553032769*e8 + 2353.37152777778*e6 + 310.71875*e4 + 26.0*e2 + 1.0,
            1: 839864.510332445*e14 + 259842.755598958*e12 + 70311.6384548611*e10 + 16012.0798611111*e8 + 2885.77083333333*e6 + 367.5*e4 + 25.0*e2,
            2: 927769.057548736*e14 + 285788.216771507*e12 + 76178.1462890625*e10 + 16688.9113769531*e8 + 2727.140625*e6 + 260.015625*e4,
            3: 969087.477329971*e14 + 290946.807098765*e12 + 73701.4461805556*e10 + 14474.4722222222*e8 + 1792.11111111111*e6,
            4: 953426.381187765*e14 + 271561.904476236*e12 + 61808.6661783854*e10 + 9661.76367865668*e8,
            5: 876003.884196429*e14 + 224786.9475*e12 + 44163.0225*e10,
            6: 716483.787616746*e14 + 179171.885103873*e12,
            7: 664011.24140688*e14,
        },
        4: {
            -7: 664011.24140688*e14,
        },
        5: {
            -7: 11503209.3055423*e14,
        },
        6: {
            -7: 123255910.704605*e14,
        },
        7: {
            -7: 937784248.676216*e14,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[3][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[3][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[3][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[3][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[3][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[3][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[3][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[3][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[3][-7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[6][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[6][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[6][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[6][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[7][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[7][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[7][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[7][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[7][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[7][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[7][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[7][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[7][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[7][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[7][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[7][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[7][6] = eccentricity_results_bymode[0][-6]

    return eccentricity_results_bymode

@njit
def eccentricity_funcs_trunc16(eccentricity: FloatArray) -> Dict[int, Dict[int, FloatArray]]:
    """ Calculates the eccentricity functions (by mode) truncated to e^16 for order-l = 7
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : Dict[int, Dict[int, FloatArray]]
    """
    # Eccentricity functions calculated at truncation level 16.
    #     and order-l = 7.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e
    e6 = e*e*e*e*e*e
    e8 = e**8
    e10 = e**10
    e12 = e**12
    e14 = e**14
    e16 = e**16

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -8: 9.38596699032984e-15*e16,
            -6: 1.87967129335142e-9*e16 + 1.14373639357363e-9*e14 + 4.7095027970679e-10*e12,
            -5: 6.23484347442681e-5*e16 + 6.28995811287478e-5*e14 + 5.78703703703704e-5*e12 + 6.94444444444444e-5*e10,
            -4: 0.00108947594250952*e16 - 0.0112746211460658*e14 + 0.0320917510986328*e12 - 0.07119140625*e10 + 0.04449462890625*e8,
            -3: -4.40869268077601*e16 + 11.5056128747795*e14 - 20.1716049382716*e12 + 20.1555555555556*e10 - 9.77777777777778*e8 + 1.77777777777778*e6,
            -2: 1348.87160996971*e16 - 2132.90065789769*e14 + 2269.74851555294*e12 - 1516.85926649306*e10 + 591.854519314236*e8 - 120.442708333333*e6 + 9.765625*e4,
            -1: -115035.60841179*e16 + 121302.196799665*e14 - 89131.943671875*e12 + 43713.97734375*e10 - 13602.09375*e8 + 2496.9375*e6 - 238.5*e4 + 9.0*e2,
            0: 3928011.79991695*e16 - 2955172.8970563*e14 + 1586997.38291086*e12 - 586640.268715278*e10 + 142589.589938694*e8 - 21364.4618055556*e6 + 1785.21875*e4 - 70.0*e2 + 1.0,
            1: -66006759.5382433*e16 + 36681269.0568395*e14 - 14695176.0341667*e12 + 4085875.57048611*e10 - 748898.888888889*e8 + 83775.8333333333*e6 - 5016.0*e4 + 121.0*e2,
            2: 614791138.403345*e16 - 255474118.187811*e14 + 76176310.2957281*e12 - 15548508.0931641*e10 + 2024000.86157227*e8 - 148518.984375*e6 + 4607.015625*e4,
            3: -3399179717.5205*e16 + 1050342781.37476*e14 - 227752424.952377*e12 + 32361560.6011285*e10 - 2669833.76736111*e8 + 95841.8402777778*e6,
            4: 11548436061.762*e16 - 2589179614.74784*e14 + 387648466.517466*e12 - 34373096.9709961*e10 + 1352599.28662788*e8,
            5: -24236865800.0458*e16 + 3741770950.29402*e14 - 347139933.975*e12 + 14501244.8025*e10,
            6: 30511860315.4628*e16 - 2910686193.59795*e14 + 126441284.408644*e12,
            7: -21052727705.1045*e16 + 937784248.676216*e14,
            8: 6104436340.31089*e16,
        },
        1: {
            -8: 0.221719882537668*e16,
            -7: 1.43968913129488*e16 + 0.119978389648684*e14,
            -6: 5.49174032711125*e16 + 0.811787874265533*e14 + 0.0649431247475706*e12,
            -5: -0.03515625*e10/(e2 - 1.0)**13,
            -4: 39.3153428329901*e16 + 9.57492866166264*e14 + 1.8597327875208*e12 + 0.257063802083333*e10 + 0.0190497504340278*e8,
            -3: 83.9811453855269*e16 + 23.2613752652392*e14 + 5.39066358024691*e12 + 0.971701388888889*e10 + 0.118055555555556*e8 + 0.00694444444444444*e6,
            -2: 191.185387344744*e16 + 64.667712789263*e14 + 20.3473989486694*e12 + 5.8693359375*e10 + 1.742431640625*e8 + 0.140625*e6 + 0.140625*e4,
            -1: 373.925813039336*e16 + 139.310944113757*e14 + 30.2398958333333*e12 + 38.0814236111111*e10 - 32.8888888888889*e8 + 24.5833333333333*e6 - 9.0*e4 + e2,
            0: 1001.54276161753*e16 - 752.495595552844*e14 + 1725.84582929258*e12 - 2021.14040798611*e10 + 1576.83668348524*e8 - 741.545138888889*e6 + 183.96875*e4 - 22.0*e2 + 1.0,
            1: -35946.0589629504*e16 + 56545.2157073103*e14 - 62900.057109375*e12 + 48573.09140625*e10 - 24466.5*e8 + 7388.4375*e6 - 1201.5*e4 + 81.0*e2,
            2: 1231228.06259342*e16 - 1266309.79616931*e14 + 938547.352011341*e12 - 475122.020203993*e10 + 152617.184597439*e8 - 27743.8177083333*e6 + 2173.890625*e4,
            3: -18573722.9224831*e16 + 13062709.1406972*e14 - 6464235.02122878*e12 + 2097863.73862847*e10 - 397130.923611111*e8 + 33184.6944444444*e6,
            4: 142257563.766493*e16 - 67933871.4747141*e14 + 21747975.9750501*e12 - 4146304.49912109*e10 + 354871.521057129*e8,
            5: -586239684.517029*e16 + 182888425.081913*e14 - 34500980.0248119*e12 + 2959618.45876736*e10,
            6: 1307624112.7669*e16 - 241659566.56067*e14 + 20518448.8791875*e12,
            7: -1477892171.3568*e16 + 123255910.704605*e14,
            8: 660185101.545129*e16,
        },
        2: {
            -8: 337.448907736202*e16,
            -7: 1653.23773598415*e16 + 156.64211077255*e14,
            -6: 5032.83500185247*e16 + 809.854400596619*e14 + 71.8156937026978*e12,
            -5: 12123.5319398803*e16 + 2563.02166482377*e14 + 390.619458912037*e12 + 32.4425173611111*e10,
            -4: 25261.633121265*e16 + 6363.25240061133*e14 + 1283.45398924368*e12 + 185.11298828125*e10 + 14.396491156684*e8,
            -3: -0.09765625*e6*(3.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            -2: 83279.1471429178*e16 + 26135.9927566851*e14 + 7178.43844695621*e12 + 1659.23604600694*e10 + 303.058295355903*e8 + 38.9322916666667*e6 + 2.640625*e4,
            -1: 137331.061940094*e16 + 46429.6170135272*e14 + 14041.6768142361*e12 + 3694.66293402778*e10 + 811.631944444444*e8 + 139.270833333333*e6 + 16.5*e4 + e2,
            0: 217273.143247915*e16 + 78327.7874431159*e14 + 25731.1676902771*e12 + 7565.199296875*e10 + 1945.65222167969*e8 + 424.84375*e6 + 76.46875*e4 + 10.0*e2 + 1.0,
            1: 332514.430101101*e16 + 126743.549716435*e14 + 44644.5701736111*e12 + 14360.1105902778*e10 + 4115.31944444444*e8 + 1080.58333333333*e6 + 189.0*e4 + 49.0*e2,
            2: 493808.611203993*e16 + 197042.647860248*e14 + 73906.7868524128*e12 + 24376.8419053819*e10 + 8894.1413031684*e8 + 1254.55729166667*e6 + 862.890625*e4,
            3: 706998.515643834*e16 + 306304.572200056*e14 + 97963.396875*e12 + 59750.0453125*e10 - 95.75*e8 + 9168.0625*e6,
            4: 1172588.42289078*e16 + 196218.790338876*e14 + 392310.672480862*e12 - 78294.8146809896*e10 + 71362.7221747504*e8,
            5: -1015685.91263426*e16 + 2608659.94096175*e14 - 894602.368952546*e12 + 448693.440434028*e10,
            6: 16725535.5538483*e16 - 6709222.00565983*e14 + 2412439.34078983*e12,
            7: -40255987.811503*e16 + 11503209.3055423*e14,
            8: 49868524.3845965*e16,
        },
        3: {
            -8: 50553.5356753209*e16,
            -7: 143177.836475606*e16 + 18617.627869898*e14,
            -6: 288970.153736016*e16 + 56701.7701593513*e14 + 6572.37797547564*e12,
            -5: 491118.834556809*e16 + 118393.583558477*e14 + 21292.203287037*e12 + 2198.82840277778*e10,
            -4: 748464.032461628*e16 + 204879.526902172*e14 + 45732.4679340363*e12 + 7485.75439453125*e10 + 685.376037597656*e8,
            -3: 1054606.98825424*e16 + 314588.587124842*e14 + 80064.9051456404*e12 + 16387.2709201389*e10 + 2413.96180555556*e8 + 193.673611111111*e6,
            -2: 1397522.06146574*e16 + 442814.440752841*e14 + 122862.972033925*e12 + 28704.1544053819*e10 + 5301.66408962674*e8 + 688.932291666667*e6 + 47.265625*e4,
            -1: -0.140625*e2*(5.0*e4 + 20.0*e2 + 8.0)**2/(e2 - 1.0)**13,
            0: 2116390.96404087*e16 + 719186.088857198*e14 + 219108.287543295*e12 + 58347.2753993056*e10 + 13103.2553032769*e8 + 2353.37152777778*e6 + 310.71875*e4 + 26.0*e2 + 1.0,
            1: 2437240.94540722*e16 + 839864.510332445*e14 + 259842.755598958*e12 + 70311.6384548611*e10 + 16012.0798611111*e8 + 2885.77083333333*e6 + 367.5*e4 + 25.0*e2,
            2: 2690460.52736734*e16 + 927769.057548736*e14 + 285788.216771507*e12 + 76178.1462890625*e10 + 16688.9113769531*e8 + 2727.140625*e6 + 260.015625*e4,
            3: 2847310.73636333*e16 + 969087.477329971*e14 + 290946.807098765*e12 + 73701.4461805556*e10 + 14474.4722222222*e8 + 1792.11111111111*e6,
            4: 2884135.51696895*e16 + 953426.381187765*e14 + 271561.904476236*e12 + 61808.6661783854*e10 + 9661.76367865668*e8,
            5: 2784035.39506975*e16 + 876003.884196429*e14 + 224786.9475*e12 + 44163.0225*e10,
            6: 2551140.62734588*e16 + 716483.787616746*e14 + 179171.885103873*e12,
            7: 2022000.82930082*e16 + 664011.24140688*e14,
            8: 2291716.29166404*e16,
        },
        4: {
            -8: 2291716.29166404*e16,
        },
        5: {
            -8: 49868524.3845965*e16,
        },
        6: {
            -8: 660185101.545129*e16,
        },
        7: {
            -8: 6104436340.31089*e16,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[3][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[3][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[3][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[3][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[3][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[3][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[3][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[3][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[3][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[3][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[3][-8]
    eccentricity_results_bymode[5][-7] = eccentricity_results_bymode[2][7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[5][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[6][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[6][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[6][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[6][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[6][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[6][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[7][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[7][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[7][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[7][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[7][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[7][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[7][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[7][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[7][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[7][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[7][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[7][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[7][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[7][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[7][8] = eccentricity_results_bymode[0][-8]

    return eccentricity_results_bymode