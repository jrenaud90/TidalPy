""" Eccentricity functions (squared) for various truncations of e at tidal order-l = 4
"""

from typing import Dict, TYPE_CHECKING

from . import EccenOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def eccentricity_funcs_trunc2(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^2 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 2.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -1: 2.25 * e2,
            0 : 1.0 - 22.0 * e2,
            1 : 42.25 * e2,
            },
        1: {
            -2: -0.5625 * e**4 / (e2 - 1.0)**7,
            -1: 0.25 * e2,
            0 : 2.0 * e2 + 1.0,
            1 : 20.25 * e2,
            },
        2: {
            -1: 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -1: 20.25 * e2,
            },
        4: {
            -1: 42.25 * e2,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc4(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^4 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 4.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -2: 0.25 * e4,
            -1: -14.0625 * e4 + 2.25 * e2,
            0 : 170.75 * e4 - 22.0 * e2 + 1.0,
            1 : -621.5625 * e4 + 42.25 * e2,
            2 : 650.25 * e4,
            },
        1: {
            -2: -0.5625 * e4 / (e2 - 1.0)**7,
            -1: 2.0625 * e4 + 0.25 * e2,
            0 : 9.125 * e4 + 2.0 * e2 + 1.0,
            1 : -1.6875 * e4 + 20.25 * e2,
            2 : 175.5625 * e4,
            },
        2: {
            -2: 25.0 * e4,
            -1: 42.1875 * e4 + 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -2: 175.5625 * e4,
            },
        4: {
            -2: 650.25 * e4,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc6(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^6 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 6.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -3: 0.000434027777777778 * e6,
            -2: -0.333333333333333 * e6 + 0.25 * e4,
            -1: 31.18359375 * e6 - 14.0625 * e4 + 2.25 * e2,
            0 : -583.638888888889 * e6 + 170.75 * e4 - 22.0 * e2 + 1.0,
            1 : 3569.95442708333 * e6 - 621.5625 * e4 + 42.25 * e2,
            2 : -8185.5 * e6 + 650.25 * e4,
            3 : 6106.77126736111 * e6,
            },
        1: {
            -3: 1.04210069444444 * e6,
            -2: -0.5625 * e4 / (e2 - 1.0)**7,
            -1: 9.11067708333333 * e6 + 2.0625 * e4 + 0.25 * e2,
            0 : 23.5694444444444 * e6 + 9.125 * e4 + 2.0 * e2 + 1.0,
            1 : 67.74609375 * e6 - 1.6875 * e4 + 20.25 * e2,
            2 : -197.645833333333 * e6 + 175.5625 * e4,
            3 : 1030.67751736111 * e6,
            },
        2: {
            -3: 82.12890625 * e6,
            -2: 129.166666666667 * e6 + 25.0 * e4,
            -1: 166.048177083333 * e6 + 42.1875 * e4 + 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -3: 1030.67751736111 * e6,
            },
        4: {
            -3: 6106.77126736111 * e6,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[0][-3]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc8(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^8 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 8.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -3: 0.000379774305555556 * e8 + 0.000434027777777778 * e6,
            -2: 0.111111111111111 * e8 - 0.333333333333333 * e6 + 0.25 * e4,
            -1: -30.61669921875 * e8 + 31.18359375 * e6 - 14.0625 * e4 + 2.25 * e2,
            0 : 1030.61024305556 * e8 - 583.638888888889 * e6 + 170.75 * e4 - 22.0 * e2 + 1.0,
            1 : -10497.2230360243 * e8 + 3569.95442708333 * e6 - 621.5625 * e4 + 42.25 * e2,
            2 : 42418.125 * e8 - 8185.5 * e6 + 650.25 * e4,
            3 : -71820.7014431424 * e8 + 6106.77126736111 * e6,
            4 : 42418.8350694444 * e8,
            },
        1: {
            -4: 1.94835069444444 * e8,
            -3: 6.76567925347222 * e8 + 1.04210069444444 * e6,
            -2: -0.5625 * e4 / (e2 - 1.0)**7,
            -1: 29.2154405381944 * e8 + 9.11067708333333 * e6 + 2.0625 * e4 + 0.25 * e2,
            0 : 58.2599826388889 * e8 + 23.5694444444444 * e6 + 9.125 * e4 + 2.0 * e2 + 1.0,
            1 : 97.09716796875 * e8 + 67.74609375 * e6 - 1.6875 * e4 + 20.25 * e2,
            2 : 460.165798611111 * e8 - 197.645833333333 * e6 + 175.5625 * e4,
            3 : -1962.95241970486 * e8 + 1030.67751736111 * e6,
            4 : 4812.890625 * e8,
            },
        2: {
            -4: 240.896267361111 * e8,
            -3: 333.82568359375 * e8 + 82.12890625 * e6,
            -2: 427.777777777778 * e8 + 129.166666666667 * e6 + 25.0 * e4,
            -1: 494.570583767361 * e8 + 166.048177083333 * e6 + 42.1875 * e4 + 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -4: 4812.890625 * e8,
            },
        4: {
            -4: 42418.8350694444 * e8,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[0][-3]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc10(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^10 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 10.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8
    e10 = e**10

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -5: 6.78168402777778e-8 * e10,
            -3: 0.000366550021701389 * e10 + 0.000379774305555556 * e8 + 0.000434027777777778 * e6,
            -2: -0.0305555555555556 * e10 + 0.111111111111111 * e8 - 0.333333333333333 * e6 + 0.25 * e4,
            -1: 15.6179992675781 * e10 - 30.61669921875 * e8 + 31.18359375 * e6 - 14.0625 * e4 + 2.25 * e2,
            0 : -1034.92486111111 * e10 + 1030.61024305556 * e8 - 583.638888888889 * e6 + 170.75 * e4 - 22.0 * e2 + 1.0,
            1 : 17934.7833421495 * e10 - 10497.2230360243 * e8 + 3569.95442708333 * e6 - 621.5625 * e4 + 42.25 * e2,
            2 : -120013.70625 * e10 + 42418.125 * e8 - 8185.5 * e6 + 650.25 * e4,
            3 : 359793.938875665 * e10 - 71820.7014431424 * e8 + 6106.77126736111 * e6,
            4 : -485876.304166667 * e10 + 42418.8350694444 * e8,
            5 : 239949.195562134 * e10,
            },
        1: {
            -5: 3.66662658691406 * e10,
            -4: 11.6203125 * e10 + 1.94835069444444 * e8,
            -3: 25.5508222791884 * e10 + 6.76567925347222 * e8 + 1.04210069444444 * e6,
            -2: -0.5625 * e4 / (e2 - 1.0)**7,
            -1: 76.7125413682726 * e10 + 29.2154405381944 * e8 + 9.11067708333333 * e6 + 2.0625 * e4 + 0.25 * e2,
            0 : 129.900034722222 * e10 + 58.2599826388889 * e8 + 23.5694444444444 * e6 + 9.125 * e4 + 2.0 * e2 + 1.0,
            1 : 213.248858642578 * e10 + 97.09716796875 * e8 + 67.74609375 * e6 - 1.6875 * e4 + 20.25 * e2,
            2 : 147.488194444444 * e10 + 460.165798611111 * e8 - 197.645833333333 * e6 + 175.5625 * e4,
            3 : 3201.53101433648 * e10 - 1962.95241970486 * e8 + 1030.67751736111 * e6,
            4 : -12399.046875 * e10 + 4812.890625 * e8,
            5 : 19310.6487826199 * e10,
            },
        2: {
            -5: 655.906780666775 * e10,
            -4: 764.401041666667 * e10 + 240.896267361111 * e8,
            -3: 966.631546020508 * e10 + 333.82568359375 * e8 + 82.12890625 * e6,
            -2: 1123.84548611111 * e10 + 427.777777777778 * e8 + 129.166666666667 * e6 + 25.0 * e4,
            -1: 1233.00030178494 * e10 + 494.570583767361 * e8 + 166.048177083333 * e6 + 42.1875 * e4 + 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -5: 19310.6487826199 * e10,
            },
        4: {
            -5: 239949.195562134 * e10,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[0][-5]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc12(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^12 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 12.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8
    e10 = e**10
    e12 = e**12

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -6: 1.92901234567901e-6 * e12,
            -5: 1.75193504050926e-7 * e12 + 6.78168402777778e-8 * e10,
            -3: 0.000341919322072724 * e12 + 0.000366550021701389 * e10 + 0.000379774305555556 * e8 + 0.000434027777777778 * e6,
            -2: -0.000998263888888889 * e12 - 0.0305555555555556 * e10 + 0.111111111111111 * e8 - 0.333333333333333 * e6 + 0.25 * e4,
            -1: -5.07659271240234 * e12 + 15.6179992675781 * e10 - 30.61669921875 * e8 + 31.18359375 * e6 - 14.0625 * e4 + 2.25 * e2,
            0 : 646.549230324074 * e12 - 1034.92486111111 * e10 + 1030.61024305556 * e8 - 583.638888888889 * e6 + 170.75 * e4 - 22.0 * e2 + 1.0,
            1 : -19313.339785258 * e12 + 17934.7833421495 * e10 - 10497.2230360243 * e8 + 3569.95442708333 * e6 - 621.5625 * e4 + 42.25 * e2,
            2 : 209706.022265625 * e12 - 120013.70625 * e10 + 42418.125 * e8 - 8185.5 * e6 + 650.25 * e4,
            3 : -1022671.96280935 * e12 + 359793.938875665 * e10 - 71820.7014431424 * e8 + 6106.77126736111 * e6,
            4 : 2428954.85556713 * e12 - 485876.304166667 * e10 + 42418.8350694444 * e8,
            5 : -2737716.23672287 * e12 + 239949.195562134 * e10,
            6 : 1168816.25004823 * e12,
            },
        1: {
            -6: 6.92530394000772 * e12,
            -5: 19.8938053894043 * e12 + 3.66662658691406 * e10,
            -4: 41.270521556713 * e12 + 11.6203125 * e10 + 1.94835069444444 * e8,
            -3: 73.2350860124753 * e12 + 25.5508222791884 * e10 + 6.76567925347222 * e8 + 1.04210069444444 * e6,
            -2: -0.5625 * e4 / (e2 - 1.0)**7,
            -1: 175.357567647298 * e12 + 76.7125413682726 * e10 + 29.2154405381944 * e8 + 9.11067708333333 * e6 + 2.0625 * e4 + 0.25 * e2,
            0 : 266.809137008102 * e12 + 129.900034722222 * e10 + 58.2599826388889 * e8 + 23.5694444444444 * e6 + 9.125 * e4 + 2.0 * e2 + 1.0,
            1 : 397.592113952637 * e12 + 213.248858642578 * e10 + 97.09716796875 * e8 + 67.74609375 * e6 - 1.6875 * e4 + 20.25 * e2,
            2 : 644.418093532986 * e12 + 147.488194444444 * e10 + 460.165798611111 * e8 - 197.645833333333 * e6 + 175.5625 * e4,
            3 : -1290.54145761184 * e12 + 3201.53101433648 * e10 - 1962.95241970486 * e8 + 1030.67751736111 * e6,
            4 : 20109.588046875 * e12 - 12399.046875 * e10 + 4812.890625 * e8,
            5 : -61669.7520816323 * e12 + 19310.6487826199 * e10,
            6 : 69524.1394102045 * e12,
            },
        2: {
            -6: 1693.8369140625 * e12,
            -5: 1570.4741746408 * e12 + 655.906780666775 * e10,
            -4: 1995.2357494213 * e12 + 764.401041666667 * e10 + 240.896267361111 * e8,
            -3: 2301.10216140747 * e12 + 966.631546020508 * e10 + 333.82568359375 * e8 + 82.12890625 * e6,
            -2: 2545.29351128472 * e12 + 1123.84548611111 * e10 + 427.777777777778 * e8 + 129.166666666667 * e6 + 25.0 * e4,
            -1: 2711.33776262071 * e12 + 1233.00030178494 * e10 + 494.570583767361 * e8 + 166.048177083333 * e6 + 42.1875 * e4 + 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -6: 69524.1394102045 * e12,
            },
        4: {
            -6: 1168816.25004823 * e12,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[0][-6]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc14(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^14 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 14.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8
    e10 = e**10
    e12 = e**12
    e14 = e**14

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -7: 1.14925540223414e-5 * e14,
            -6: 5.51146384479718e-6 * e14 + 1.92901234567901e-6 * e12,
            -5: 2.88440226667562e-7 * e14 + 1.75193504050926e-7 * e12 + 6.78168402777778e-8 * e10,
            -3: 0.000316369310678418 * e14 + 0.000341919322072724 * e12 + 0.000366550021701389 * e10 + 0.000379774305555556 * e8 + 0.000434027777777778 * e6,
            -2: -0.00308807319223986 * e14 - 0.000998263888888889 * e12 - 0.0305555555555556 * e10 + 0.111111111111111 * e8 - 0.333333333333333 * e6 + 0.25 * e4,
            -1: 1.05155106272016 * e14 - 5.07659271240234 * e12 + 15.6179992675781 * e10 - 30.61669921875 * e8 + 31.18359375 * e6 - 14.0625 * e4 + 2.25 * e2,
            0 : -274.421647415911 * e14 + 646.549230324074 * e12 - 1034.92486111111 * e10 + 1030.61024305556 * e8 - 583.638888888889 * e6 + 170.75 * e4 - 22.0 * e2 + 1.0,
            1 : 14027.1271612972 * e14 - 19313.339785258 * e12 + 17934.7833421495 * e10 - 10497.2230360243 * e8 + 3569.95442708333 * e6 - 621.5625 * e4 + 42.25 * e2,
            2 : -244265.454441964 * e14 + 209706.022265625 * e12 - 120013.70625 * e10 + 42418.125 * e8 - 8185.5 * e6 + 650.25 * e4,
            3 : 1864648.11796509 * e14 - 1022671.96280935 * e12 + 359793.938875665 * e10 - 71820.7014431424 * e8 + 6106.77126736111 * e6,
            4 : -7073073.58282352 * e14 + 2428954.85556713 * e12 - 485876.304166667 * e10 + 42418.8350694444 * e8,
            5 : 13885275.9401448 * e14 - 2737716.23672287 * e12 + 239949.195562134 * e10,
            6 : -13459442.5123663 * e14 + 1168816.25004823 * e12,
            7 : 5080086.80850045 * e14,
            },
        1: {
            -7: 13.1045189610835 * e14,
            -6: 33.824509118028 * e14 + 6.92530394000772 * e12,
            -5: 66.2469145815713 * e14 + 19.8938053894043 * e12 + 3.66662658691406 * e10,
            -4: 112.829789220955 * e14 + 41.270521556713 * e12 + 11.6203125 * e10 + 1.94835069444444 * e8,
            -3: 176.401245009929 * e14 + 73.2350860124753 * e12 + 25.5508222791884 * e10 + 6.76567925347222 * e8 + 1.04210069444444 * e6,
            -2: -0.5625 * e4 / (e2 - 1.0)**7,
            -1: 361.787971839249 * e14 + 175.357567647298 * e12 + 76.7125413682726 * e10 + 29.2154405381944 * e8 + 9.11067708333333 * e6 + 2.0625 * e4 + 0.25 * e2,
            0 : 510.991944606836 * e14 + 266.809137008102 * e12 + 129.900034722222 * e10 + 58.2599826388889 * e8 + 23.5694444444444 * e6 + 9.125 * e4 + 2.0 * e2 + 1.0,
            1 : 715.117373326165 * e14 + 397.592113952637 * e12 + 213.248858642578 * e10 + 97.09716796875 * e8 + 67.74609375 * e6 - 1.6875 * e4 + 20.25 * e2,
            2 : 958.088995001791 * e14 + 644.418093532986 * e12 + 147.488194444444 * e10 + 460.165798611111 * e8 - 197.645833333333 * e6 + 175.5625 * e4,
            3 : 2405.01878383467 * e14 - 1290.54145761184 * e12 + 3201.53101433648 * e10 - 1962.95241970486 * e8 + 1030.67751736111 * e6,
            4 : -15810.425859375 * e14 + 20109.588046875 * e12 - 12399.046875 * e10 + 4812.890625 * e8,
            5 : 109198.707912504 * e14 - 61669.7520816323 * e12 + 19310.6487826199 * e10,
            6 : -262758.011185585 * e14 + 69524.1394102045 * e12,
            7 : 230739.513376018 * e14,
            },
        2: {
            -7: 4203.9590960638 * e14,
            -6: 2851.90764508929 * e14 + 1693.8369140625 * e12,
            -5: 3856.25669750377 * e14 + 1570.4741746408 * e12 + 655.906780666775 * e10,
            -4: 4371.62514812059 * e14 + 1995.2357494213 * e12 + 764.401041666667 * e10 + 240.896267361111 * e8,
            -3: 4829.9199915954 * e14 + 2301.10216140747 * e12 + 966.631546020508 * e10 + 333.82568359375 * e8 + 82.12890625 * e6,
            -2: 5187.27980668541 * e14 + 2545.29351128472 * e12 + 1123.84548611111 * e10 + 427.777777777778 * e8 + 129.166666666667 * e6 + 25.0 * e4,
            -1: 5426.83725947299 * e14 + 2711.33776262071 * e12 + 1233.00030178494 * e10 + 494.570583767361 * e8 + 166.048177083333 * e6 + 42.1875 * e4 + 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -7: 230739.513376018 * e14,
            },
        4: {
            -7: 5080086.80850045 * e14,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[0][-7]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc16(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^16 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 16.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8
    e10 = e**10
    e12 = e**12
    e14 = e**14
    e16 = e**16

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -8: 4.03124212648022e-5 * e16,
            -7: 3.37593774406277e-5 * e16 + 1.14925540223414e-5 * e14,
            -6: 9.89603017132779e-6 * e16 + 5.51146384479718e-6 * e14 + 1.92901234567901e-6 * e12,
            -5: 3.89251969688994e-7 * e16 + 2.88440226667562e-7 * e14 + 1.75193504050926e-7 * e12 + 6.78168402777778e-8 * e10,
            -3: 0.000291948672745383 * e16 + 0.000316369310678418 * e14 + 0.000341919322072724 * e12 + 0.000366550021701389 * e10 + 0.000379774305555556 * e8 + 0.000434027777777778 * e6,
            -2: -0.00191220238095238 * e16 - 0.00308807319223986 * e14 - 0.000998263888888889 * e12 - 0.0305555555555556 * e10 + 0.111111111111111 * e8 - 0.333333333333333 * e6 + 0.25 * e4,
            -1: -0.210739373321436 * e16 + 1.05155106272016 * e14 - 5.07659271240234 * e12 + 15.6179992675781 * e10 - 30.61669921875 * e8 + 31.18359375 * e6 - 14.0625 * e4 + 2.25 * e2,
            0 : 83.7341605691484 * e16 - 274.421647415911 * e14 + 646.549230324074 * e12 - 1034.92486111111 * e10 + 1030.61024305556 * e8 - 583.638888888889 * e6 + 170.75 * e4 - 22.0 * e2 + 1.0,
            1 : -7287.22223439182 * e16 + 14027.1271612972 * e14 - 19313.339785258 * e12 + 17934.7833421495 * e10 - 10497.2230360243 * e8 + 3569.95442708333 * e6 - 621.5625 * e4 + 42.25 * e2,
            2 : 200995.432721819 * e16 - 244265.454441964 * e14 + 209706.022265625 * e12 - 120013.70625 * e10 + 42418.125 * e8 - 8185.5 * e6 + 650.25 * e4,
            3 : -2347634.97832888 * e16 + 1864648.11796509 * e14 - 1022671.96280935 * e12 + 359793.938875665 * e10 - 71820.7014431424 * e8 + 6106.77126736111 * e6,
            4 : 13560151.1710948 * e16 - 7073073.58282352 * e14 + 2428954.85556713 * e12 - 485876.304166667 * e10 + 42418.8350694444 * e8,
            5 : -41820478.6210939 * e16 + 13885275.9401448 * e14 - 2737716.23672287 * e12 + 239949.195562134 * e10,
            6 : 69895407.2647045 * e16 - 13459442.5123663 * e14 + 1168816.25004823 * e12,
            7 : -59518599.4607124 * e16 + 5080086.80850045 * e14,
            8 : 20181795.1394302 * e16,
            },
        1: {
            -8: 24.8184113047074 * e16,
            -7: 56.9153849110289 * e16 + 13.1045189610835 * e14,
            -6: 105.429229452451 * e16 + 33.824509118028 * e14 + 6.92530394000772 * e12,
            -5: 172.58139504305 * e16 + 66.2469145815713 * e14 + 19.8938053894043 * e12 + 3.66662658691406 * e10,
            -4: 261.681921161739 * e16 + 112.829789220955 * e14 + 41.270521556713 * e12 + 11.6203125 * e10 + 1.94835069444444 * e8,
            -3: 376.234101765403 * e16 + 176.401245009929 * e14 + 73.2350860124753 * e12 + 25.5508222791884 * e10 + 6.76567925347222 * e8 + 1.04210069444444 * e6,
            -2: -0.5625 * e4 / (e2 - 1.0)**7,
            -1: 689.540891238697 * e16 + 361.787971839249 * e14 + 175.357567647298 * e12 + 76.7125413682726 * e10 + 29.2154405381944 * e8 + 9.11067708333333 * e6 + 2.0625 * e4 + 0.25 * e2,
            0 : 922.540749251985 * e16 + 510.991944606836 * e14 + 266.809137008102 * e12 + 129.900034722222 * e10 + 58.2599826388889 * e8 + 23.5694444444444 * e6 + 9.125 * e4 + 2.0 * e2 + 1.0,
            1 : 1228.54824858416 * e16 + 715.117373326165 * e14 + 397.592113952637 * e12 + 213.248858642578 * e10 + 97.09716796875 * e8 + 67.74609375 * e6 - 1.6875 * e4 + 20.25 * e2,
            2 : 1619.29401058459 * e16 + 958.088995001791 * e14 + 644.418093532986 * e12 + 147.488194444444 * e10 + 460.165798611111 * e8 - 197.645833333333 * e6 + 175.5625 * e4,
            3 : 1659.24271223582 * e16 + 2405.01878383467 * e14 - 1290.54145761184 * e12 + 3201.53101433648 * e10 - 1962.95241970486 * e8 + 1030.67751736111 * e6,
            4 : 13946.7877608817 * e16 - 15810.425859375 * e14 + 20109.588046875 * e12 - 12399.046875 * e10 + 4812.890625 * e8,
            5 : -112211.596750438 * e16 + 109198.707912504 * e14 - 61669.7520816323 * e12 + 19310.6487826199 * e10,
            6 : 518979.357460918 * e16 - 262758.011185585 * e14 + 69524.1394102045 * e12,
            7 : -1002607.24567437 * e16 + 230739.513376018 * e14,
            8 : 718688.367614519 * e16,
            },
        2: {
            -8: 10115.7471973691 * e16,
            -7: 4278.53728709011 * e16 + 4203.9590960638 * e14,
            -6: 7150.50067263233 * e16 + 2851.90764508929 * e14 + 1693.8369140625 * e12,
            -5: 7828.84497193997 * e16 + 3856.25669750377 * e14 + 1570.4741746408 * e12 + 655.906780666775 * e10,
            -4: 8613.101790117 * e16 + 4371.62514812059 * e14 + 1995.2357494213 * e12 + 764.401041666667 * e10 + 240.896267361111 * e8,
            -3: 9260.66661050703 * e16 + 4829.9199915954 * e14 + 2301.10216140747 * e12 + 966.631546020508 * e10 + 333.82568359375 * e8 + 82.12890625 * e6,
            -2: 9760.96896830564 * e16 + 5187.27980668541 * e14 + 2545.29351128472 * e12 + 1123.84548611111 * e10 + 427.777777777778 * e8 + 129.166666666667 * e6 + 25.0 * e4,
            -1: 10092.7562355593 * e16 + 5426.83725947299 * e14 + 2711.33776262071 * e12 + 1233.00030178494 * e10 + 494.570583767361 * e8 + 166.048177083333 * e6 + 42.1875 * e4 + 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -8: 718688.367614519 * e16,
            },
        4: {
            -8: 20181795.1394302 * e16,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[0][-8]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc18(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^18 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 18.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8
    e10 = e**10
    e12 = e**12
    e14 = e**14
    e16 = e**16
    e18 = e**18

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -9: 0.000110508094485363 * e18,
            -8: 0.000116458105876095 * e18 + 4.03124212648022e-5 * e16,
            -7: 6.29509135852663e-5 * e18 + 3.37593774406277e-5 * e16 + 1.14925540223414e-5 * e14,
            -6: 1.43602429826526e-5 * e18 + 9.89603017132779e-6 * e16 + 5.51146384479718e-6 * e14 + 1.92901234567901e-6 * e12,
            -5: 4.71147415257258e-7 * e18 + 3.89251969688994e-7 * e16 + 2.88440226667562e-7 * e14 + 1.75193504050926e-7 * e12 + 6.78168402777778e-8 * e10,
            -3: 0.000269459533033123 * e18 + 0.000291948672745383 * e16 + 0.000316369310678418 * e14 + 0.000341919322072724 * e12 + 0.000366550021701389 * e10 + 0.000379774305555556 * e8 + 0.000434027777777778 * e6,
            -2: -0.00132808918642045 * e18 - 0.00191220238095238 * e16 - 0.00308807319223986 * e14 - 0.000998263888888889 * e12 - 0.0305555555555556 * e10 + 0.111111111111111 * e8 - 0.333333333333333 * e6 + 0.25 * e4,
            -1: -0.00624679957763577 * e18 - 0.210739373321436 * e16 + 1.05155106272016 * e14 - 5.07659271240234 * e12 + 15.6179992675781 * e10 - 30.61669921875 * e8 + 31.18359375 * e6 - 14.0625 * e4 + 2.25 * e2,
            0 : -19.6095487951493 * e18 + 83.7341605691484 * e16 - 274.421647415911 * e14 + 646.549230324074 * e12 - 1034.92486111111 * e10 + 1030.61024305556 * e8 - 583.638888888889 * e6 + 170.75 * e4 - 22.0 * e2 + 1.0,
            1 : 2835.36492864121 * e18 - 7287.22223439182 * e16 + 14027.1271612972 * e14 - 19313.339785258 * e12 + 17934.7833421495 * e10 - 10497.2230360243 * e8 + 3569.95442708333 * e6 - 621.5625 * e4 + 42.25 * e2,
            2 : -122475.1553694 * e18 + 200995.432721819 * e16 - 244265.454441964 * e14 + 209706.022265625 * e12 - 120013.70625 * e10 + 42418.125 * e8 - 8185.5 * e6 + 650.25 * e4,
            3 : 2151920.00018377 * e18 - 2347634.97832888 * e16 + 1864648.11796509 * e14 - 1022671.96280935 * e12 + 359793.938875665 * e10 - 71820.7014431424 * e8 + 6106.77126736111 * e6,
            4 : -18391093.9448307 * e18 + 13560151.1710948 * e16 - 7073073.58282352 * e14 + 2428954.85556713 * e12 - 485876.304166667 * e10 + 42418.8350694444 * e8,
            5 : 84511835.2637969 * e18 - 41820478.6210939 * e16 + 13885275.9401448 * e14 - 2737716.23672287 * e12 + 239949.195562134 * e10,
            6 : -218755505.070142 * e18 + 69895407.2647045 * e16 - 13459442.5123663 * e14 + 1168816.25004823 * e12,
            7 : 318124431.806029 * e18 - 59518599.4607124 * e16 + 5080086.80850045 * e14,
            8 : -241783123.600429 * e18 + 20181795.1394302 * e16,
            9 : 74550951.0043993 * e18,
            },
        1: {
            -9: 47.0155082797035 * e18,
            -8: 94.4097103196747 * e18 + 24.8184113047074 * e16,
            -7: 166.038771294539 * e18 + 56.9153849110289 * e16 + 13.1045189610835 * e14,
            -6: 261.673075206366 * e18 + 105.429229452451 * e16 + 33.824509118028 * e14 + 6.92530394000772 * e12,
            -5: 385.330895126831 * e18 + 172.58139504305 * e16 + 66.2469145815713 * e14 + 19.8938053894043 * e12 + 3.66662658691406 * e10,
            -4: 541.013072546844 * e18 + 261.681921161739 * e16 + 112.829789220955 * e14 + 41.270521556713 * e12 + 11.6203125 * e10 + 1.94835069444444 * e8,
            -3: 732.945076815878 * e18 + 376.234101765403 * e16 + 176.401245009929 * e14 + 73.2350860124753 * e12 + 25.5508222791884 * e10 + 6.76567925347222 * e8 + 1.04210069444444 * e6,
            -2: -0.5625 * e4 / (e2 - 1.0)**7,
            -1: 1233.62829513752 * e18 + 689.540891238697 * e16 + 361.787971839249 * e14 + 175.357567647298 * e12 + 76.7125413682726 * e10 + 29.2154405381944 * e8 + 9.11067708333333 * e6 + 2.0625 * e4 + 0.25 * e2,
            0 : 1584.24247010098 * e18 + 922.540749251985 * e16 + 510.991944606836 * e14 + 266.809137008102 * e12 + 129.900034722222 * e10 + 58.2599826388889 * e8 + 23.5694444444444 * e6 + 9.125 * e4 + 2.0 * e2 + 1.0,
            1 : 2029.47967558932 * e18 + 1228.54824858416 * e16 + 715.117373326165 * e14 + 397.592113952637 * e12 + 213.248858642578 * e10 + 97.09716796875 * e8 + 67.74609375 * e6 - 1.6875 * e4 + 20.25 * e2,
            2 : 2577.36986845289 * e18 + 1619.29401058459 * e16 + 958.088995001791 * e14 + 644.418093532986 * e12 + 147.488194444444 * e10 + 460.165798611111 * e8 - 197.645833333333 * e6 + 175.5625 * e4,
            3 : 3365.37682166583 * e18 + 1659.24271223582 * e16 + 2405.01878383467 * e14 - 1290.54145761184 * e12 + 3201.53101433648 * e10 - 1962.95241970486 * e8 + 1030.67751736111 * e6,
            4 : -1443.24249389648 * e18 + 13946.7877608817 * e16 - 15810.425859375 * e14 + 20109.588046875 * e12 - 12399.046875 * e10 + 4812.890625 * e8,
            5 : 93297.3472560036 * e18 - 112211.596750438 * e16 + 109198.707912504 * e14 - 61669.7520816323 * e12 + 19310.6487826199 * e10,
            6 : -632971.254111344 * e18 + 518979.357460918 * e16 - 262758.011185585 * e14 + 69524.1394102045 * e12,
            7 : 2209279.34743102 * e18 - 1002607.24567437 * e16 + 230739.513376018 * e14,
            8 : -3519755.67639182 * e18 + 718688.367614519 * e16,
            9 : 2127365.75708672 * e18,
            },
        2: {
            -9: 23743.0948315122 * e18,
            -8: 3873.16651570201 * e18 + 10115.7471973691 * e16,
            -7: 13182.0511517127 * e18 + 4278.53728709011 * e16 + 4203.9590960638 * e14,
            -6: 13313.5535776566 * e18 + 7150.50067263233 * e16 + 2851.90764508929 * e14 + 1693.8369140625 * e12,
            -5: 14625.2196112 * e18 + 7828.84497193997 * e16 + 3856.25669750377 * e14 + 1570.4741746408 * e12 + 655.906780666775 * e10,
            -4: 15689.441801609 * e18 + 8613.101790117 * e16 + 4371.62514812059 * e14 + 1995.2357494213 * e12 + 764.401041666667 * e10 + 240.896267361111 * e8,
            -3: 16572.114646754 * e18 + 9260.66661050703 * e16 + 4829.9199915954 * e14 + 2301.10216140747 * e12 + 966.631546020508 * e10 + 333.82568359375 * e8 + 82.12890625 * e6,
            -2: 17248.5291050075 * e18 + 9760.96896830564 * e16 + 5187.27980668541 * e14 + 2545.29351128472 * e12 + 1123.84548611111 * e10 + 427.777777777778 * e8 + 129.166666666667 * e6 + 25.0 * e4,
            -1: 17693.354767452 * e18 + 10092.7562355593 * e16 + 5426.83725947299 * e14 + 2711.33776262071 * e12 + 1233.00030178494 * e10 + 494.570583767361 * e8 + 166.048177083333 * e6 + 42.1875 * e4 + 6.25 * e2,
            0 : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -9: 2127365.75708672 * e18,
            },
        4: {
            -9: 74550951.0043993 * e18,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[2][9] = eccentricity_results_bymode[2][-9]
    eccentricity_results_bymode[3][-8] = eccentricity_results_bymode[1][8]
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[3][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[4][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[4][9] = eccentricity_results_bymode[0][-9]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc20(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^20 for order-l = 4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20.
    #     and order-l = 4.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8
    e10 = e**10
    e12 = e**12
    e14 = e**14
    e16 = e**16
    e18 = e**18
    e20 = e**20

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -10: 0.000264788444674745 * e20,
            -9 : 0.000303897259834749 * e20 + 0.000110508094485363 * e18,
            -8 : 0.000219155243258712 * e20 + 0.000116458105876095 * e18 + 4.03124212648022e-5 * e16,
            -7 : 9.48803487081765e-5 * e20 + 6.29509135852663e-5 * e18 + 3.37593774406277e-5 * e16 + 1.14925540223414e-5 * e14,
            -6 : 1.84668086761668e-5 * e20 + 1.43602429826526e-5 * e18 + 9.89603017132779e-6 * e16 + 5.51146384479718e-6 * e14 + 1.92901234567901e-6 * e12,
            -5 : 5.33776651940077e-7 * e20 + 4.71147415257258e-7 * e18 + 3.89251969688994e-7 * e16 + 2.88440226667562e-7 * e14 + 1.75193504050926e-7 * e12 + 6.78168402777778e-8 * e10,
            -3 : 0.000249089377844448 * e20 + 0.000269459533033123 * e18 + 0.000291948672745383 * e16 + 0.000316369310678418 * e14 + 0.000341919322072724 * e12 + 0.000366550021701389 * e10 + 0.000379774305555556 * e8 + 0.000434027777777778 * e6,
            -2 : -0.000900413512258068 * e20 - 0.00132808918642045 * e18 - 0.00191220238095238 * e16 - 0.00308807319223986 * e14 - 0.000998263888888889 * e12 - 0.0305555555555556 * e10 + 0.111111111111111 * e8 - 0.333333333333333 * e6 + 0.25 * e4,
            -1 : -0.0243703577957705 * e20 - 0.00624679957763577 * e18 - 0.210739373321436 * e16 + 1.05155106272016 * e14 - 5.07659271240234 * e12 + 15.6179992675781 * e10 - 30.61669921875 * e8 + 31.18359375 * e6 - 14.0625 * e4 + 2.25 * e2,
            0  : 3.41340665092355 * e20 - 19.6095487951493 * e18 + 83.7341605691484 * e16 - 274.421647415911 * e14 + 646.549230324074 * e12 - 1034.92486111111 * e10 + 1030.61024305556 * e8 - 583.638888888889 * e6 + 170.75 * e4 - 22.0 * e2 + 1.0,
            1  : -858.491204550473 * e20 + 2835.36492864121 * e18 - 7287.22223439182 * e16 + 14027.1271612972 * e14 - 19313.339785258 * e12 + 17934.7833421495 * e10 - 10497.2230360243 * e8 + 3569.95442708333 * e6 - 621.5625 * e4 + 42.25 * e2,
            2  : 57396.1703826711 * e20 - 122475.1553694 * e18 + 200995.432721819 * e16 - 244265.454441964 * e14 + 209706.022265625 * e12 - 120013.70625 * e10 + 42418.125 * e8 - 8185.5 * e6 + 650.25 * e4,
            3  : -1496508.19846058 * e20 + 2151920.00018377 * e18 - 2347634.97832888 * e16 + 1864648.11796509 * e14 - 1022671.96280935 * e12 + 359793.938875665 * e10 - 71820.7014431424 * e8 + 6106.77126736111 * e6,
            4  : 18548423.7749046 * e20 - 18391093.9448307 * e18 + 13560151.1710948 * e16 - 7073073.58282352 * e14 + 2428954.85556713 * e12 - 485876.304166667 * e10 + 42418.8350694444 * e8,
            5  : -122970261.610087 * e20 + 84511835.2637969 * e18 - 41820478.6210939 * e16 + 13885275.9401448 * e14 - 2737716.23672287 * e12 + 239949.195562134 * e10,
            6  : 466038806.296261 * e20 - 218755505.070142 * e18 + 69895407.2647045 * e16 - 13459442.5123663 * e14 + 1168816.25004823 * e12,
            7  : -1036833911.39453 * e20 + 318124431.806029 * e18 - 59518599.4607124 * e16 + 5080086.80850045 * e14,
            8  : 1334060842.90354 * e20 - 241783123.600429 * e18 + 20181795.1394302 * e16,
            9  : -916177658.262001 * e20 + 74550951.0043993 * e18,
            10 : 259321364.360223 * e20,
            },
        1: {
            -10: 89.0575497114019 * e20,
            -9 : 153.608232871895 * e20 + 47.0155082797035 * e18,
            -8 : 258.427641819857 * e20 + 94.4097103196747 * e18 + 24.8184113047074 * e16,
            -7 : 392.81916766517 * e20 + 166.038771294539 * e18 + 56.9153849110289 * e16 + 13.1045189610835 * e14,
            -6 : 562.648882297901 * e20 + 261.673075206366 * e18 + 105.429229452451 * e16 + 33.824509118028 * e14 + 6.92530394000772 * e12,
            -5 : 772.35536512039 * e20 + 385.330895126831 * e18 + 172.58139504305 * e16 + 66.2469145815713 * e14 + 19.8938053894043 * e12 + 3.66662658691406 * e10,
            -4 : 1026.73231107527 * e20 + 541.013072546844 * e18 + 261.681921161739 * e16 + 112.829789220955 * e14 + 41.270521556713 * e12 + 11.6203125 * e10 + 1.94835069444444 * e8,
            -3 : 1330.76366549821 * e20 + 732.945076815878 * e18 + 376.234101765403 * e16 + 176.401245009929 * e14 + 73.2350860124753 * e12 + 25.5508222791884 * e10 + 6.76567925347222 * e8 + 1.04210069444444 * e6,
            -2 : -0.5625 * e4 / (e2 - 1.0)**7,
            -1 : 2095.67032556917 * e20 + 1233.62829513752 * e18 + 689.540891238697 * e16 + 361.787971839249 * e14 + 175.357567647298 * e12 + 76.7125413682726 * e10 + 29.2154405381944 * e8 + 9.11067708333333 * e6 + 2.0625 * e4 + 0.25 * e2,
            0  : 2606.78747425335 * e20 + 1584.24247010098 * e18 + 922.540749251985 * e16 + 510.991944606836 * e14 + 266.809137008102 * e12 + 129.900034722222 * e10 + 58.2599826388889 * e8 + 23.5694444444444 * e6 + 9.125 * e4 + 2.0 * e2 + 1.0,
            1  : 3237.83712084324 * e20 + 2029.47967558932 * e18 + 1228.54824858416 * e16 + 715.117373326165 * e14 + 397.592113952637 * e12 + 213.248858642578 * e10 + 97.09716796875 * e8 + 67.74609375 * e6 - 1.6875 * e4 + 20.25 * e2,
            2  : 3999.99210257598 * e20 + 2577.36986845289 * e18 + 1619.29401058459 * e16 + 958.088995001791 * e14 + 644.418093532986 * e12 + 147.488194444444 * e10 + 460.165798611111 * e8 - 197.645833333333 * e6 + 175.5625 * e4,
            3  : 4870.85267815359 * e20 + 3365.37682166583 * e18 + 1659.24271223582 * e16 + 2405.01878383467 * e14 - 1290.54145761184 * e12 + 3201.53101433648 * e10 - 1962.95241970486 * e8 + 1030.67751736111 * e6,
            4  : 8007.57559164865 * e20 - 1443.24249389648 * e18 + 13946.7877608817 * e16 - 15810.425859375 * e14 + 20109.588046875 * e12 - 12399.046875 * e10 + 4812.890625 * e8,
            5  : -43938.7933704723 * e20 + 93297.3472560036 * e18 - 112211.596750438 * e16 + 109198.707912504 * e14 - 61669.7520816323 * e12 + 19310.6487826199 * e10,
            6  : 574669.555148506 * e20 - 632971.254111344 * e18 + 518979.357460918 * e16 - 262758.011185585 * e14 + 69524.1394102045 * e12,
            7  : -3081356.18575181 * e20 + 2209279.34743102 * e18 - 1002607.24567437 * e16 + 230739.513376018 * e14,
            8  : 8600912.63852483 * e20 - 3519755.67639182 * e18 + 718688.367614519 * e16,
            9  : -11573173.6099864 * e20 + 2127365.75708672 * e18,
            10 : 6039431.01813882 * e20,
            },
        2: {
            -10: 54600.977190567 * e20,
            -9 : -5445.46991862148 * e20 + 23743.0948315122 * e18,
            -8 : 25519.9999508399 * e20 + 3873.16651570201 * e18 + 10115.7471973691 * e16,
            -7 : 21405.5811719681 * e20 + 13182.0511517127 * e18 + 4278.53728709011 * e16 + 4203.9590960638 * e14,
            -6 : 23867.8924682306 * e20 + 13313.5535776566 * e18 + 7150.50067263233 * e16 + 2851.90764508929 * e14 + 1693.8369140625 * e12,
            -5 : 25488.0510220735 * e20 + 14625.2196112 * e18 + 7828.84497193997 * e16 + 3856.25669750377 * e14 + 1570.4741746408 * e12 + 655.906780666775 * e10,
            -4 : 26907.7728517542 * e20 + 15689.441801609 * e18 + 8613.101790117 * e16 + 4371.62514812059 * e14 + 1995.2357494213 * e12 + 764.401041666667 * e10 + 240.896267361111 * e8,
            -3 : 28075.2628323042 * e20 + 16572.114646754 * e18 + 9260.66661050703 * e16 + 4829.9199915954 * e14 + 2301.10216140747 * e12 + 966.631546020508 * e10 + 333.82568359375 * e8 + 82.12890625 * e6,
            -2 : 28964.3814115806 * e20 + 17248.5291050075 * e18 + 9760.96896830564 * e16 + 5187.27980668541 * e14 + 2545.29351128472 * e12 + 1123.84548611111 * e10 + 427.777777777778 * e8 + 129.166666666667 * e6 + 25.0 * e4,
            -1 : 29545.1451755686 * e20 + 17693.354767452 * e18 + 10092.7562355593 * e16 + 5426.83725947299 * e14 + 2711.33776262071 * e12 + 1233.00030178494 * e10 + 494.570583767361 * e8 + 166.048177083333 * e6 + 42.1875 * e4 + 6.25 * e2,
            0  : -0.25 * (3.0 * e2 + 2.0)**2 / (e2 - 1.0)**7,
            },
        3: {
            -10: 6039431.01813882 * e20,
            },
        4: {
            -10: 259321364.360223 * e20,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[2][9] = eccentricity_results_bymode[2][-9]
    eccentricity_results_bymode[2][10] = eccentricity_results_bymode[2][-10]
    eccentricity_results_bymode[3][-9] = eccentricity_results_bymode[1][9]
    eccentricity_results_bymode[3][-8] = eccentricity_results_bymode[1][8]
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[3][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[3][10] = eccentricity_results_bymode[1][-10]
    eccentricity_results_bymode[4][-9] = eccentricity_results_bymode[0][9]
    eccentricity_results_bymode[4][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[4][9] = eccentricity_results_bymode[0][-9]
    eccentricity_results_bymode[4][10] = eccentricity_results_bymode[0][-10]

    return eccentricity_results_bymode
