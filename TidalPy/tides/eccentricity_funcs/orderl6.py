""" Eccentricity functions (squared) for various truncations of e at tidal order-l = 6
"""

from typing import Dict, TYPE_CHECKING

from . import EccenOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def eccentricity_funcs_trunc2(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^2 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 2.
    #     and order-l = 6.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -1: 6.25 * e2,
            0 : 1.0 - 51.0 * e2,
            1 : 90.25 * e2,
            },
        1: {
            -1: 0.25 * e2,
            0 : 1.0 - 11.0 * e2,
            1 : 56.25 * e2,
            },
        2: {
            -2: -1.5625 * e**4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 2.25 * e2,
            0 : 13.0 * e2 + 1.0,
            1 : 30.25 * e2,
            },
        3: {
            -1: 12.25 * e2,
            0 : -0.015625 * (15.0 * e**4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -1: 30.25 * e2,
            },
        5: {
            -1: 56.25 * e2,
            },
        6: {
            -1: 90.25 * e2,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc4(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^4 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 4.
    #     and order-l = 6.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -2: 4.0 * e4,
            -1: -114.0625 * e4 + 6.25 * e2,
            0 : 941.625 * e4 - 51.0 * e2 + 1.0,
            1 : -2801.3125 * e4 + 90.25 * e2,
            2 : 2652.25 * e4,
            },
        1: {
            -2: 0.0625 * e4,
            -1: -1.9375 * e4 + 0.25 * e2,
            0 : 52.25 * e4 - 11.0 * e2 + 1.0,
            1 : -473.4375 * e4 + 56.25 * e2,
            2 : 1105.5625 * e4,
            },
        2: {
            -2: -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 30.1875 * e4 + 2.25 * e2,
            0 : 94.625 * e4 + 13.0 * e2 + 1.0,
            1 : 190.4375 * e4 + 30.25 * e2,
            2 : 361.0 * e4,
            },
        3: {
            -2: 76.5625 * e4,
            -1: 162.3125 * e4 + 12.25 * e2,
            0 : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -2: 361.0 * e4,
            },
        5: {
            -2: 1105.5625 * e4,
            },
        6: {
            -2: 2652.25 * e4,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc6(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^6 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 6.
    #     and order-l = 6.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -3: 0.31640625 * e6,
            -2: -30.6666666666667 * e6 + 4.0 * e4,
            -1: 808.756510416667 * e6 - 114.0625 * e4 + 6.25 * e2,
            0 : -8060.0 * e6 + 941.625 * e4 - 51.0 * e2 + 1.0,
            1 : 34797.1002604167 * e6 - 2801.3125 * e4 + 90.25 * e2,
            2 : -65473.6666666667 * e6 + 2652.25 * e4,
            3 : 43759.41015625 * e6,
            },
        1: {
            -4: -0.09765625 * e**8 / (e2 - 1.0)**11,
            -3: 0.0525173611111111 * e6,
            -2: 0.5 * e6 + 0.0625 * e4,
            -1: 1.37109375 * e6 - 1.9375 * e4 + 0.25 * e2,
            0 : -117.222222222222 * e6 + 52.25 * e4 - 11.0 * e2 + 1.0,
            1 : 1769.43359375 * e6 - 473.4375 * e4 + 56.25 * e2,
            2 : -8711.5 * e6 + 1105.5625 * e4,
            3 : 12858.6150173611 * e6,
            },
        2: {
            -3: 16.1671006944444 * e6,
            -2: -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 212.49609375 * e6 + 30.1875 * e4 + 2.25 * e2,
            0 : 491.111111111111 * e6 + 94.625 * e4 + 13.0 * e2 + 1.0,
            1 : 909.27734375 * e6 + 190.4375 * e4 + 30.25 * e2,
            2 : 1330.0 * e6 + 361.0 * e4,
            3 : 2767.19835069444 * e6,
            },
        3: {
            -3: 353.91015625 * e6,
            -2: 775.833333333333 * e6 + 76.5625 * e4,
            -1: 1130.38151041667 * e6 + 162.3125 * e4 + 12.25 * e2,
            0 : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -3: 2767.19835069444 * e6,
            },
        5: {
            -3: 12858.6150173611 * e6,
            },
        6: {
            -3: 43759.41015625 * e6,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[0][-3]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc8(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^8 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 8.
    #     and order-l = 6.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -4: 0.00173611111111111 * e8,
            -3: -0.83056640625 * e8 + 0.31640625 * e6,
            -2: 89.6111111111111 * e8 - 30.6666666666667 * e6 + 4.0 * e4,
            -1: -2910.81678602431 * e8 + 808.756510416667 * e6 - 114.0625 * e4 + 6.25 * e2,
            0 : 37913.58984375 * e8 - 8060.0 * e6 + 941.625 * e4 - 51.0 * e2 + 1.0,
            1 : -229156.124728733 * e8 + 34797.1002604167 * e6 - 2801.3125 * e4 + 90.25 * e2,
            2 : 679077.819444444 * e8 - 65473.6666666667 * e6 + 2652.25 * e4,
            3 : -951271.984863281 * e8 + 43759.41015625 * e6,
            4 : 500703.656684028 * e8,
            },
        1: {
            -4: -0.09765625 * e8 / (e2 - 1.0)**11,
            -3: 0.605740017361111 * e8 + 0.0525173611111111 * e6,
            -2: 2.76302083333333 * e8 + 0.5 * e6 + 0.0625 * e4,
            -1: 1.62646484375 * e8 + 1.37109375 * e6 - 1.9375 * e4 + 0.25 * e2,
            0 : 130.514756944444 * e8 - 117.222222222222 * e6 + 52.25 * e4 - 11.0 * e2 + 1.0,
            1 : -3489.35139973958 * e8 + 1769.43359375 * e6 - 473.4375 * e4 + 56.25 * e2,
            2 : 30777.9140625 * e8 - 8711.5 * e6 + 1105.5625 * e4,
            3 : -101014.721082899 * e8 + 12858.6150173611 * e6,
            4 : 108035.47265625 * e8,
            },
        2: {
            -4: 39.84765625 * e8,
            -3: 175.900987413194 * e8 + 16.1671006944444 * e6,
            -2: -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 1049.76285807292 * e8 + 212.49609375 * e6 + 30.1875 * e4 + 2.25 * e2,
            0 : 2020.14366319444 * e8 + 491.111111111111 * e6 + 94.625 * e4 + 13.0 * e2 + 1.0,
            1 : 3412.37548828125 * e8 + 909.27734375 * e6 + 190.4375 * e4 + 30.25 * e2,
            2 : 5180.95833333333 * e8 + 1330.0 * e6 + 361.0 * e4,
            3 : 5826.87038845486 * e8 + 2767.19835069444 * e6,
            4 : 16256.25 * e8,
            },
        3: {
            -4: 1372.08506944444 * e8,
            -3: 2916.67236328125 * e8 + 353.91015625 * e6,
            -2: 4417.99652777778 * e8 + 775.833333333333 * e6 + 76.5625 * e4,
            -1: 5531.57416449653 * e8 + 1130.38151041667 * e6 + 162.3125 * e4 + 12.25 * e2,
            0 : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -4: 16256.25 * e8,
            },
        5: {
            -4: 108035.47265625 * e8,
            },
        6: {
            -4: 500703.656684028 * e8,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[0][-4]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc10(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^10 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 10.
    #     and order-l = 6.

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
            -4: 0.000347222222222222 * e10 + 0.00173611111111111 * e8,
            -3: 0.696340942382813 * e10 - 0.83056640625 * e8 + 0.31640625 * e6,
            -2: -128.8 * e10 + 89.6111111111111 * e8 - 30.6666666666667 * e6 + 4.0 * e4,
            -1: 6002.7229309082 * e10 - 2910.81678602431 * e8 + 808.756510416667 * e6 - 114.0625 * e4 + 6.25 * e2,
            0 : -108050.40984375 * e10 + 37913.58984375 * e8 - 8060.0 * e6 + 941.625 * e4 - 51.0 * e2 + 1.0,
            1 : 910799.197898356 * e10 - 229156.124728733 * e8 + 34797.1002604167 * e6 - 2801.3125 * e4 + 90.25 * e2,
            2 : -3940116.19166667 * e10 + 679077.819444444 * e8 - 65473.6666666667 * e6 + 2652.25 * e4,
            3 : 8950487.88320618 * e10 - 951271.984863281 * e8 + 43759.41015625 * e6,
            4 : -10092278.1358507 * e10 + 500703.656684028 * e8,
            5 : 4432418.57569994 * e10,
            },
        1: {
            -5: 0.180403713650174 * e10,
            -4: -0.09765625 * e8 / (e2 - 1.0)**11,
            -3: 3.78169793023003 * e10 + 0.605740017361111 * e8 + 0.0525173611111111 * e6,
            -2: 11.6824652777778 * e10 + 2.76302083333333 * e8 + 0.5 * e6 + 0.0625 * e4,
            -1: 18.1928527832031 * e10 + 1.62646484375 * e8 + 1.37109375 * e6 - 1.9375 * e4 + 0.25 * e2,
            0 : -68.4535069444444 * e10 + 130.514756944444 * e8 - 117.222222222222 * e6 + 52.25 * e4 - 11.0 * e2 + 1.0,
            1 : 4141.78814358181 * e10 - 3489.35139973958 * e8 + 1769.43359375 * e6 - 473.4375 * e4 + 56.25 * e2,
            2 : -61175.3546875 * e10 + 30777.9140625 * e8 - 8711.5 * e6 + 1105.5625 * e4,
            3 : 359143.306014336 * e10 - 101014.721082899 * e8 + 12858.6150173611 * e6,
            4 : -868382.78828125 * e10 + 108035.47265625 * e8,
            5 : 725442.129741821 * e10,
            },
        2: {
            -5: 94.7732672119141 * e10,
            -4: 395.504427083333 * e10 + 39.84765625 * e8,
            -3: 1046.12310418023 * e10 + 175.900987413194 * e8 + 16.1671006944444 * e6,
            -2: -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 4098.99631415473 * e10 + 1049.76285807292 * e8 + 212.49609375 * e6 + 30.1875 * e4 + 2.25 * e2,
            0 : 6983.08883680556 * e10 + 2020.14366319444 * e8 + 491.111111111111 * e6 + 94.625 * e4 + 13.0 * e2 + 1.0,
            1 : 10944.0567199707 * e10 + 3412.37548828125 * e8 + 909.27734375 * e6 + 190.4375 * e4 + 30.25 * e2,
            2 : 15846.3694444444 * e10 + 5180.95833333333 * e8 + 1330.0 * e6 + 361.0 * e4,
            3 : 22338.8989342584 * e10 + 5826.87038845486 * e8 + 2767.19835069444 * e6,
            4 : 15357.375 * e10 + 16256.25 * e8,
            5 : 79813.1655143907 * e10,
            },
        3: {
            -5: 4734.26420254178 * e10,
            -4: 9414.44826388889 * e10 + 1372.08506944444 * e8,
            -3: 14262.7660675049 * e10 + 2916.67236328125 * e8 + 353.91015625 * e6,
            -2: 18484.9947916667 * e10 + 4417.99652777778 * e8 + 775.833333333333 * e6 + 76.5625 * e4,
            -1: 21415.4602681478 * e10 + 5531.57416449653 * e8 + 1130.38151041667 * e6 + 162.3125 * e4 + 12.25 * e2,
            0 : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -5: 79813.1655143907 * e10,
            },
        5: {
            -5: 725442.129741821 * e10,
            },
        6: {
            -5: 4432418.57569994 * e10,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[0][-5]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc12(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^12 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 12.
    #     and order-l = 6.

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
            -5: 1.29982277199074e-7 * e12 + 6.78168402777778e-8 * e10,
            -4: 0.000552662037037037 * e12 + 0.000347222222222222 * e10 + 0.00173611111111111 * e8,
            -3: -0.268883514404297 * e12 + 0.696340942382813 * e10 - 0.83056640625 * e8 + 0.31640625 * e6,
            -2: 102.849652777778 * e12 - 128.8 * e10 + 89.6111111111111 * e8 - 30.6666666666667 * e6 + 4.0 * e4,
            -1: -7628.348520067 * e12 + 6002.7229309082 * e10 - 2910.81678602431 * e8 + 808.756510416667 * e6 - 114.0625 * e4 + 6.25 * e2,
            0 : 198843.244707031 * e12 - 108050.40984375 * e10 + 37913.58984375 * e8 - 8060.0 * e6 + 941.625 * e4 - 51.0 * e2 + 1.0,
            1 : -2357802.65863441 * e12 + 910799.197898356 * e10 - 229156.124728733 * e8 + 34797.1002604167 * e6 - 2801.3125 * e4 + 90.25 * e2,
            2 : 14457660.8707465 * e12 - 3940116.19166667 * e10 + 679077.819444444 * e8 - 65473.6666666667 * e6 + 2652.25 * e4,
            3 : -48600122.2348045 * e12 + 8950487.88320618 * e10 - 951271.984863281 * e8 + 43759.41015625 * e6,
            4 : 89871728.1598126 * e12 - 10092278.1358507 * e10 + 500703.656684028 * e8,
            5 : -85358456.808223 * e12 + 4432418.57569994 * e10,
            6 : 32405694.76 * e12,
            },
        1: {
            -6: 0.333426046489198 * e12,
            -5: 1.8942113410102 * e12 + 0.180403713650174 * e10,
            -4: -0.09765625 * e8 / (e2 - 1.0)**11,
            -3: 16.9546108198755 * e12 + 3.78169793023003 * e10 + 0.605740017361111 * e8 + 0.0525173611111111 * e6,
            -2: 41.1977593315972 * e12 + 11.6824652777778 * e10 + 2.76302083333333 * e8 + 0.5 * e6 + 0.0625 * e4,
            -1: 68.1135433959961 * e12 + 18.1928527832031 * e10 + 1.62646484375 * e8 + 1.37109375 * e6 - 1.9375 * e4 + 0.25 * e2,
            0 : 130.753796296296 * e12 - 68.4535069444444 * e10 + 130.514756944444 * e8 - 117.222222222222 * e6 + 52.25 * e4 - 11.0 * e2 + 1.0,
            1 : -3053.98570166694 * e12 + 4141.78814358181 * e10 - 3489.35139973958 * e8 + 1769.43359375 * e6 - 473.4375 * e4 + 56.25 * e2,
            2 : 77679.9017089844 * e12 - 61175.3546875 * e10 + 30777.9140625 * e8 - 8711.5 * e6 + 1105.5625 * e4,
            3 : -745768.239982143 * e12 + 359143.306014336 * e10 - 101014.721082899 * e8 + 12858.6150173611 * e6,
            4 : 3188086.33699653 * e12 - 868382.78828125 * e10 + 108035.47265625 * e8,
            5 : -6038281.67580887 * e12 + 725442.129741821 * e10,
            6 : 4131091.53653369 * e12,
            },
        2: {
            -6: 219.327872299383 * e12,
            -5: 860.110561676025 * e12 + 94.7732672119141 * e10,
            -4: 2185.01003255208 * e12 + 395.504427083333 * e10 + 39.84765625 * e8,
            -3: 4500.02279358852 * e12 + 1046.12310418023 * e10 + 175.900987413194 * e8 + 16.1671006944444 * e6,
            -2: -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 13491.608223504 * e12 + 4098.99631415473 * e10 + 1049.76285807292 * e8 + 212.49609375 * e6 + 30.1875 * e4 + 2.25 * e2,
            0 : 21076.3362550637 * e12 + 6983.08883680556 * e10 + 2020.14366319444 * e8 + 491.111111111111 * e6 + 94.625 * e4 + 13.0 * e2 + 1.0,
            1 : 31090.9299649048 * e12 + 10944.0567199707 * e10 + 3412.37548828125 * e8 + 909.27734375 * e6 + 190.4375 * e4 + 30.25 * e2,
            2 : 43423.2953993056 * e12 + 15846.3694444444 * e10 + 5180.95833333333 * e8 + 1330.0 * e6 + 361.0 * e4,
            3 : 57148.0307461303 * e12 + 22338.8989342584 * e10 + 5826.87038845486 * e8 + 2767.19835069444 * e6,
            4 : 84145.69125 * e12 + 15357.375 * e10 + 16256.25 * e8,
            5 : -407.969009518094 * e12 + 79813.1655143907 * e10,
            6 : 344140.297316744 * e12,
            },
        3: {
            -6: 15027.69515625 * e12,
            -5: 27211.1487738998 * e12 + 4734.26420254178 * e10,
            -4: 40533.1550231481 * e12 + 9414.44826388889 * e10 + 1372.08506944444 * e8,
            -3: 52945.9238456726 * e12 + 14262.7660675049 * e10 + 2916.67236328125 * e8 + 353.91015625 * e6,
            -2: 63143.2810980903 * e12 + 18484.9947916667 * e10 + 4417.99652777778 * e8 + 775.833333333333 * e6 + 76.5625 * e4,
            -1: 69942.9386810642 * e12 + 21415.4602681478 * e10 + 5531.57416449653 * e8 + 1130.38151041667 * e6 + 162.3125 * e4 + 12.25 * e2,
            0 : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -6: 344140.297316744 * e12,
            },
        5: {
            -6: 4131091.53653369 * e12,
            },
        6: {
            -6: 32405694.76 * e12,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[3][-6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[6][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[0][-5]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc14(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^14 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 14.
    #     and order-l = 6.

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
            -7: 2.40280754952444e-12 * e14,
            -5: 1.82274577899374e-7 * e14 + 1.29982277199074e-7 * e12 + 6.78168402777778e-8 * e10,
            -4: 0.000423831569664903 * e14 + 0.000552662037037037 * e12 + 0.000347222222222222 * e10 + 0.00173611111111111 * e8,
            -3: 0.0818452875954764 * e14 - 0.268883514404297 * e12 + 0.696340942382813 * e10 - 0.83056640625 * e8 + 0.31640625 * e6,
            -2: -51.4749018959436 * e14 + 102.849652777778 * e12 - 128.8 * e10 + 89.6111111111111 * e8 - 30.6666666666667 * e6 + 4.0 * e4,
            -1: 6344.85087845691 * e14 - 7628.348520067 * e12 + 6002.7229309082 * e10 - 2910.81678602431 * e8 + 808.756510416667 * e6 - 114.0625 * e4 + 6.25 * e2,
            0 : -247941.106069037 * e14 + 198843.244707031 * e12 - 108050.40984375 * e10 + 37913.58984375 * e8 - 8060.0 * e6 + 941.625 * e4 - 51.0 * e2 + 1.0,
            1 : 4187599.83676463 * e14 - 2357802.65863441 * e12 + 910799.197898356 * e10 - 229156.124728733 * e8 + 34797.1002604167 * e6 - 2801.3125 * e4 + 90.25 * e2,
            2 : -36008756.8366826 * e14 + 14457660.8707465 * e12 - 3940116.19166667 * e10 + 679077.819444444 * e8 - 65473.6666666667 * e6 + 2652.25 * e4,
            3 : 171836856.188733 * e14 - 48600122.2348045 * e12 + 8950487.88320618 * e10 - 951271.984863281 * e8 + 43759.41015625 * e6,
            4 : -471706639.341561 * e14 + 89871728.1598126 * e12 - 10092278.1358507 * e10 + 500703.656684028 * e8,
            5 : 737419222.764519 * e14 - 85358456.808223 * e12 + 4432418.57569994 * e10,
            6 : -607649410.257857 * e14 + 32405694.76 * e12,
            7 : 204185799.366173 * e14,
            },
        1: {
            -7: 0.616886645141913 * e14,
            -6: 3.33357304756393 * e14 + 0.333426046489198 * e12,
            -5: 10.9254498859688 * e14 + 1.8942113410102 * e12 + 0.180403713650174 * e10,
            -4: -0.09765625 * e8 / (e2 - 1.0)**11,
            -3: 61.124635589825 * e14 + 16.9546108198755 * e12 + 3.78169793023003 * e10 + 0.605740017361111 * e8 + 0.0525173611111111 * e6,
            -2: 126.394356708829 * e14 + 41.1977593315972 * e12 + 11.6824652777778 * e10 + 2.76302083333333 * e8 + 0.5 * e6 + 0.0625 * e4,
            -1: 208.170076889311 * e14 + 68.1135433959961 * e12 + 18.1928527832031 * e10 + 1.62646484375 * e8 + 1.37109375 * e6 - 1.9375 * e4 + 0.25 * e2,
            0 : 281.988204601285 * e14 + 130.753796296296 * e12 - 68.4535069444444 * e10 + 130.514756944444 * e8 - 117.222222222222 * e6 + 52.25 * e4 - 11.0 * e2 + 1.0,
            1 : 2109.77532191882 * e14 - 3053.98570166694 * e12 + 4141.78814358181 * e10 - 3489.35139973958 * e8 + 1769.43359375 * e6 - 473.4375 * e4 + 56.25 * e2,
            2 : -66859.2437695313 * e14 + 77679.9017089844 * e12 - 61175.3546875 * e10 + 30777.9140625 * e8 - 8711.5 * e6 + 1105.5625 * e4,
            3 : 1022138.69151758 * e14 - 745768.239982143 * e12 + 359143.306014336 * e10 - 101014.721082899 * e8 + 12858.6150173611 * e6,
            4 : -7000489.52094411 * e14 + 3188086.33699653 * e12 - 868382.78828125 * e10 + 108035.47265625 * e8,
            5 : 23151140.4241943 * e14 - 6038281.67580887 * e12 + 725442.129741821 * e10,
            6 : -35811174.726429 * e14 + 4131091.53653369 * e12,
            7 : 20725438.1895036 * e14,
            },
        2: {
            -7: 496.705916584755 * e14,
            -6: 1819.56243496473 * e14 + 219.327872299383 * e12,
            -5: 4428.17964847701 * e14 + 860.110561676025 * e12 + 94.7732672119141 * e10,
            -4: 8845.59142629795 * e14 + 2185.01003255208 * e12 + 395.504427083333 * e10 + 39.84765625 * e8,
            -3: 15651.9454827029 * e14 + 4500.02279358852 * e12 + 1046.12310418023 * e10 + 175.900987413194 * e8 + 16.1671006944444 * e6,
            -2: -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 38957.2608385778 * e14 + 13491.608223504 * e12 + 4098.99631415473 * e10 + 1049.76285807292 * e8 + 212.49609375 * e6 + 30.1875 * e4 + 2.25 * e2,
            0 : 57049.3071968596 * e14 + 21076.3362550637 * e12 + 6983.08883680556 * e10 + 2020.14366319444 * e8 + 491.111111111111 * e6 + 94.625 * e4 + 13.0 * e2 + 1.0,
            1 : 80147.647066352 * e14 + 31090.9299649048 * e12 + 10944.0567199707 * e10 + 3412.37548828125 * e8 + 909.27734375 * e6 + 190.4375 * e4 + 30.25 * e2,
            2 : 108181.088184524 * e14 + 43423.2953993056 * e12 + 15846.3694444444 * e10 + 5180.95833333333 * e8 + 1330.0 * e6 + 361.0 * e4,
            3 : 140797.275858066 * e14 + 57148.0307461303 * e12 + 22338.8989342584 * e10 + 5826.87038845486 * e8 + 2767.19835069444 * e6,
            4 : 169084.547075893 * e14 + 84145.69125 * e12 + 15357.375 * e10 + 16256.25 * e8,
            5 : 314179.867007791 * e14 - 407.969009518094 * e12 + 79813.1655143907 * e10,
            6 : -286610.403578593 * e14 + 344140.297316744 * e12,
            7 : 1344523.62270713 * e14,
            },
        3: {
            -7: 44792.9712759683 * e14,
            -6: 72003.30046875 * e14 + 15027.69515625 * e12,
            -5: 104830.957339307 * e14 + 27211.1487738998 * e12 + 4734.26420254178 * e10,
            -4: 136234.820096451 * e14 + 40533.1550231481 * e12 + 9414.44826388889 * e10 + 1372.08506944444 * e8,
            -3: 164107.183336115 * e14 + 52945.9238456726 * e12 + 14262.7660675049 * e10 + 2916.67236328125 * e8 + 353.91015625 * e6,
            -2: 186188.84331115 * e14 + 63143.2810980903 * e12 + 18484.9947916667 * e10 + 4417.99652777778 * e8 + 775.833333333333 * e6 + 76.5625 * e4,
            -1: 200534.725930072 * e14 + 69942.9386810642 * e12 + 21415.4602681478 * e10 + 5531.57416449653 * e8 + 1130.38151041667 * e6 + 162.3125 * e4 + 12.25 * e2,
            0 : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -7: 1344523.62270713 * e14,
            },
        5: {
            -7: 20725438.1895036 * e14,
            },
        6: {
            -7: 204185799.366173 * e14,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[3][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[3][-7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[6][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[6][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[6][7] = eccentricity_results_bymode[0][-7]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc16(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^16 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 16.
    #     and order-l = 6.

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
            -8: 6.15118732678257e-10 * e16,
            -7: 8.86035283887137e-12 * e16 + 2.40280754952444e-12 * e14,
            -5: 2.23275635398044e-7 * e16 + 1.82274577899374e-7 * e14 + 1.29982277199074e-7 * e12 + 6.78168402777778e-8 * e10,
            -4: 0.000364889047343474 * e16 + 0.000423831569664903 * e14 + 0.000552662037037037 * e12 + 0.000347222222222222 * e10 + 0.00173611111111111 * e8,
            -3: -0.00410451778343746 * e16 + 0.0818452875954764 * e14 - 0.268883514404297 * e12 + 0.696340942382813 * e10 - 0.83056640625 * e8 + 0.31640625 * e6,
            -2: 17.8718484071869 * e16 - 51.4749018959436 * e14 + 102.849652777778 * e12 - 128.8 * e10 + 89.6111111111111 * e8 - 30.6666666666667 * e6 + 4.0 * e4,
            -1: -3669.99875953992 * e16 + 6344.85087845691 * e14 - 7628.348520067 * e12 + 6002.7229309082 * e10 - 2910.81678602431 * e8 + 808.756510416667 * e6 - 114.0625 * e4 + 6.25 * e2,
            0 : 218602.008076882 * e16 - 247941.106069037 * e14 + 198843.244707031 * e12 - 108050.40984375 * e10 + 37913.58984375 * e8 - 8060.0 * e6 + 941.625 * e4 - 51.0 * e2 + 1.0,
            1 : -5312856.42342945 * e16 + 4187599.83676463 * e14 - 2357802.65863441 * e12 + 910799.197898356 * e10 - 229156.124728733 * e8 + 34797.1002604167 * e6 - 2801.3125 * e4 + 90.25 * e2,
            2 : 63871121.2968106 * e16 - 36008756.8366826 * e14 + 14457660.8707465 * e12 - 3940116.19166667 * e10 + 679077.819444444 * e8 - 65473.6666666667 * e6 + 2652.25 * e4,
            3 : -423808849.698558 * e16 + 171836856.188733 * e14 - 48600122.2348045 * e12 + 8950487.88320618 * e10 - 951271.984863281 * e8 + 43759.41015625 * e6,
            4 : 1645491640.43578 * e16 - 471706639.341561 * e14 + 89871728.1598126 * e12 - 10092278.1358507 * e10 + 500703.656684028 * e8,
            5 : -3813865880.14918 * e16 + 737419222.764519 * e14 - 85358456.808223 * e12 + 4432418.57569994 * e10,
            6 : 5172892203.72458 * e16 - 607649410.257857 * e14 + 32405694.76 * e12,
            7 : -3775523947.20958 * e16 + 204185799.366173 * e14,
            8 : 1142015861.10956 * e16,
            },
        1: {
            -8: 1.14273024524661 * e16,
            -7: 5.85592993927245 * e16 + 0.616886645141913 * e14,
            -6: 18.4574599376245 * e16 + 3.33357304756393 * e14 + 0.333426046489198 * e12,
            -5: 45.7630135075178 * e16 + 10.9254498859688 * e14 + 1.8942113410102 * e12 + 0.180403713650174 * e10,
            -4: -0.09765625 * e8 / (e2 - 1.0)**11,
            -3: 188.217183966076 * e16 + 61.124635589825 * e14 + 16.9546108198755 * e12 + 3.78169793023003 * e10 + 0.605740017361111 * e8 + 0.0525173611111111 * e6,
            -2: 346.866849481492 * e16 + 126.394356708829 * e14 + 41.1977593315972 * e12 + 11.6824652777778 * e10 + 2.76302083333333 * e8 + 0.5 * e6 + 0.0625 * e4,
            -1: 553.680295894864 * e16 + 208.170076889311 * e14 + 68.1135433959961 * e12 + 18.1928527832031 * e10 + 1.62646484375 * e8 + 1.37109375 * e6 - 1.9375 * e4 + 0.25 * e2,
            0 : 792.65418115502 * e16 + 281.988204601285 * e14 + 130.753796296296 * e12 - 68.4535069444444 * e10 + 130.514756944444 * e8 - 117.222222222222 * e6 + 52.25 * e4 - 11.0 * e2 + 1.0,
            1 : 369.822816342269 * e16 + 2109.77532191882 * e14 - 3053.98570166694 * e12 + 4141.78814358181 * e10 - 3489.35139973958 * e8 + 1769.43359375 * e6 - 473.4375 * e4 + 56.25 * e2,
            2 : 44089.1411303711 * e16 - 66859.2437695313 * e14 + 77679.9017089844 * e12 - 61175.3546875 * e10 + 30777.9140625 * e8 - 8711.5 * e6 + 1105.5625 * e4,
            3 : -988542.074944686 * e16 + 1022138.69151758 * e14 - 745768.239982143 * e12 + 359143.306014336 * e10 - 101014.721082899 * e8 + 12858.6150173611 * e6,
            4 : 10369096.3572048 * e16 - 7000489.52094411 * e14 + 3188086.33699653 * e12 - 868382.78828125 * e10 + 108035.47265625 * e8,
            5 : -53965266.192221 * e16 + 23151140.4241943 * e14 - 6038281.67580887 * e12 + 725442.129741821 * e10,
            6 : 144077742.421739 * e16 - 35811174.726429 * e14 + 4131091.53653369 * e12,
            7 : -187589612.523247 * e16 + 20725438.1895036 * e14,
            8 : 94024315.6193818 * e16,
            },
        2: {
            -8: 1105.2582240638 * e16,
            -7: 3758.40170018442 * e16 + 496.705916584755 * e14,
            -6: 8745.56848149821 * e16 + 1819.56243496473 * e14 + 219.327872299383 * e12,
            -5: 16923.8945209788 * e16 + 4428.17964847701 * e14 + 860.110561676025 * e12 + 94.7732672119141 * e10,
            -4: 29238.0507502426 * e16 + 8845.59142629795 * e14 + 2185.01003255208 * e12 + 395.504427083333 * e10 + 39.84765625 * e8,
            -3: 46702.2873911873 * e16 + 15651.9454827029 * e14 + 4500.02279358852 * e12 + 1046.12310418023 * e10 + 175.900987413194 * e8 + 16.1671006944444 * e6,
            -2: -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 101373.176843053 * e16 + 38957.2608385778 * e14 + 13491.608223504 * e12 + 4098.99631415473 * e10 + 1049.76285807292 * e8 + 212.49609375 * e6 + 30.1875 * e4 + 2.25 * e2,
            0 : 141220.816219917 * e16 + 57049.3071968596 * e14 + 21076.3362550637 * e12 + 6983.08883680556 * e10 + 2020.14366319444 * e8 + 491.111111111111 * e6 + 94.625 * e4 + 13.0 * e2 + 1.0,
            1 : 190675.584203681 * e16 + 80147.647066352 * e14 + 31090.9299649048 * e12 + 10944.0567199707 * e10 + 3412.37548828125 * e8 + 909.27734375 * e6 + 190.4375 * e4 + 30.25 * e2,
            2 : 249801.02931992 * e16 + 108181.088184524 * e14 + 43423.2953993056 * e12 + 15846.3694444444 * e10 + 5180.95833333333 * e8 + 1330.0 * e6 + 361.0 * e4,
            3 : 317908.286921933 * e16 + 140797.275858066 * e14 + 57148.0307461303 * e12 + 22338.8989342584 * e10 + 5826.87038845486 * e8 + 2767.19835069444 * e6,
            4 : 397402.852876674 * e16 + 169084.547075893 * e14 + 84145.69125 * e12 + 15357.375 * e10 + 16256.25 * e8,
            5 : 396921.039788735 * e16 + 314179.867007791 * e14 - 407.969009518094 * e12 + 79813.1655143907 * e10,
            6 : 1270114.88737564 * e16 - 286610.403578593 * e14 + 344140.297316744 * e12,
            7 : -2128818.75773916 * e16 + 1344523.62270713 * e14,
            8 : 4860787.77403088 * e16,
            },
        3: {
            -8: 127095.901675536 * e16,
            -7: 176327.753474399 * e16 + 44792.9712759683 * e14,
            -6: 251722.407875977 * e16 + 72003.30046875 * e14 + 15027.69515625 * e12,
            -5: 323221.433356364 * e16 + 104830.957339307 * e14 + 27211.1487738998 * e12 + 4734.26420254178 * e10,
            -4: 389453.523706838 * e16 + 136234.820096451 * e14 + 40533.1550231481 * e12 + 9414.44826388889 * e10 + 1372.08506944444 * e8,
            -3: 446353.172189103 * e16 + 164107.183336115 * e14 + 52945.9238456726 * e12 + 14262.7660675049 * e10 + 2916.67236328125 * e8 + 353.91015625 * e6,
            -2: 490358.387215742 * e16 + 186188.84331115 * e14 + 63143.2810980903 * e12 + 18484.9947916667 * e10 + 4417.99652777778 * e8 + 775.833333333333 * e6 + 76.5625 * e4,
            -1: 518442.206898254 * e16 + 200534.725930072 * e14 + 69942.9386810642 * e12 + 21415.4602681478 * e10 + 5531.57416449653 * e8 + 1130.38151041667 * e6 + 162.3125 * e4 + 12.25 * e2,
            0 : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -8: 4860787.77403088 * e16,
            },
        5: {
            -8: 94024315.6193818 * e16,
            },
        6: {
            -8: 1142015861.10956 * e16,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[3][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[3][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[3][-8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[2][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[5][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[5][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[6][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[6][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[6][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[6][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[6][8] = eccentricity_results_bymode[0][-8]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc18(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^18 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 18.
    #     and order-l = 6.

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
            -9: 1.12231972874427e-8 * e18,
            -8: 2.52882145656617e-9 * e18 + 6.15118732678257e-10 * e16,
            -7: 1.94834424227877e-11 * e18 + 8.86035283887137e-12 * e16 + 2.40280754952444e-12 * e14,
            -5: 2.54098263498964e-7 * e18 + 2.23275635398044e-7 * e16 + 1.82274577899374e-7 * e14 + 1.29982277199074e-7 * e12 + 6.78168402777778e-8 * e10,
            -4: 0.000314514650542573 * e18 + 0.000364889047343474 * e16 + 0.000423831569664903 * e14 + 0.000552662037037037 * e12 + 0.000347222222222222 * e10 + 0.00173611111111111 * e8,
            -3: 0.00951790638800178 * e18 - 0.00410451778343746 * e16 + 0.0818452875954764 * e14 - 0.268883514404297 * e12 + 0.696340942382813 * e10 - 0.83056640625 * e8 + 0.31640625 * e6,
            -2: -4.47036443461129 * e18 + 17.8718484071869 * e16 - 51.4749018959436 * e14 + 102.849652777778 * e12 - 128.8 * e10 + 89.6111111111111 * e8 - 30.6666666666667 * e6 + 4.0 * e4,
            -1: 1556.52804495489 * e18 - 3669.99875953992 * e16 + 6344.85087845691 * e14 - 7628.348520067 * e12 + 6002.7229309082 * e10 - 2910.81678602431 * e8 + 808.756510416667 * e6 - 114.0625 * e4 + 6.25 * e2,
            0 : -141903.055105815 * e18 + 218602.008076882 * e16 - 247941.106069037 * e14 + 198843.244707031 * e12 - 108050.40984375 * e10 + 37913.58984375 * e8 - 8060.0 * e6 + 941.625 * e4 - 51.0 * e2 + 1.0,
            1 : 4985882.10295459 * e18 - 5312856.42342945 * e16 + 4187599.83676463 * e14 - 2357802.65863441 * e12 + 910799.197898356 * e10 - 229156.124728733 * e8 + 34797.1002604167 * e6 - 2801.3125 * e4 + 90.25 * e2,
            2 : -83708417.4232456 * e18 + 63871121.2968106 * e16 - 36008756.8366826 * e14 + 14457660.8707465 * e12 - 3940116.19166667 * e10 + 679077.819444444 * e8 - 65473.6666666667 * e6 + 2652.25 * e4,
            3 : 763624009.156845 * e18 - 423808849.698558 * e16 + 171836856.188733 * e14 - 48600122.2348045 * e12 + 8950487.88320618 * e10 - 951271.984863281 * e8 + 43759.41015625 * e6,
            4 : -4083236400.87082 * e18 + 1645491640.43578 * e16 - 471706639.341561 * e14 + 89871728.1598126 * e12 - 10092278.1358507 * e10 + 500703.656684028 * e8,
            5 : 13311984041.219 * e18 - 3813865880.14918 * e16 + 737419222.764519 * e14 - 85358456.808223 * e12 + 4432418.57569994 * e10,
            6 : -26684504470.8569 * e18 + 5172892203.72458 * e16 - 607649410.257857 * e14 + 32405694.76 * e12,
            7 : 32001453962.3975 * e18 - 3775523947.20958 * e16 + 204185799.366173 * e14,
            8 : -21006631484.0241 * e18 + 1142015861.10956 * e16,
            9 : 5792215177.98438 * e18,
            },
        1: {
            -9: 2.11929181104972 * e18,
            -8: 10.2664729798742 * e18 + 1.14273024524661 * e16,
            -7: 31.076514676145 * e18 + 5.85592993927245 * e16 + 0.616886645141913 * e14,
            -6: 74.6715785950487 * e18 + 18.4574599376245 * e16 + 3.33357304756393 * e14 + 0.333426046489198 * e12,
            -5: 155.504597129133 * e18 + 45.7630135075178 * e16 + 10.9254498859688 * e14 + 1.8942113410102 * e12 + 0.180403713650174 * e10,
            -4: -0.09765625 * e8 / (e2 - 1.0)**11,
            -3: 513.674032637722 * e18 + 188.217183966076 * e16 + 61.124635589825 * e14 + 16.9546108198755 * e12 + 3.78169793023003 * e10 + 0.605740017361111 * e8 + 0.0525173611111111 * e6,
            -2: 868.670156063517 * e18 + 346.866849481492 * e16 + 126.394356708829 * e14 + 41.1977593315972 * e12 + 11.6824652777778 * e10 + 2.76302083333333 * e8 + 0.5 * e6 + 0.0625 * e4,
            -1: 1335.2030446067 * e18 + 553.680295894864 * e16 + 208.170076889311 * e14 + 68.1135433959961 * e12 + 18.1928527832031 * e10 + 1.62646484375 * e8 + 1.37109375 * e6 - 1.9375 * e4 + 0.25 * e2,
            0 : 1885.7486840863 * e18 + 792.65418115502 * e16 + 281.988204601285 * e14 + 130.753796296296 * e12 - 68.4535069444444 * e10 + 130.514756944444 * e8 - 117.222222222222 * e6 + 52.25 * e4 - 11.0 * e2 + 1.0,
            1 : 2734.02074166065 * e18 + 369.822816342269 * e16 + 2109.77532191882 * e14 - 3053.98570166694 * e12 + 4141.78814358181 * e10 - 3489.35139973958 * e8 + 1769.43359375 * e6 - 473.4375 * e4 + 56.25 * e2,
            2 : -17290.9875345285 * e18 + 44089.1411303711 * e16 - 66859.2437695313 * e14 + 77679.9017089844 * e12 - 61175.3546875 * e10 + 30777.9140625 * e8 - 8711.5 * e6 + 1105.5625 * e4,
            3 : 720960.855457223 * e18 - 988542.074944686 * e16 + 1022138.69151758 * e14 - 745768.239982143 * e12 + 359143.306014336 * e10 - 101014.721082899 * e8 + 12858.6150173611 * e6,
            4 : -11092205.6537999 * e18 + 10369096.3572048 * e16 - 7000489.52094411 * e14 + 3188086.33699653 * e12 - 868382.78828125 * e10 + 108035.47265625 * e8,
            5 : 86189199.8191689 * e18 - 53965266.192221 * e16 + 23151140.4241943 * e14 - 6038281.67580887 * e12 + 725442.129741821 * e10,
            6 : -356665472.701488 * e18 + 144077742.421739 * e16 - 35811174.726429 * e14 + 4131091.53653369 * e12,
            7 : 793257308.259555 * e18 - 187589612.523247 * e16 + 20725438.1895036 * e14,
            8 : -889350074.882826 * e18 + 94024315.6193818 * e16,
            9 : 392991619.772678 * e18,
            },
        2: {
            -9: 2423.74006354441 * e18,
            -8: 7597.54040902274 * e18 + 1105.2582240638 * e16,
            -7: 16881.7174939228 * e18 + 3758.40170018442 * e16 + 496.705916584755 * e14,
            -6: 31624.7449636542 * e18 + 8745.56848149821 * e16 + 1819.56243496473 * e14 + 219.327872299383 * e12,
            -5: 53311.4774134068 * e18 + 16923.8945209788 * e16 + 4428.17964847701 * e14 + 860.110561676025 * e12 + 94.7732672119141 * e10,
            -4: 83529.8813267346 * e18 + 29238.0507502426 * e16 + 8845.59142629795 * e14 + 2185.01003255208 * e12 + 395.504427083333 * e10 + 39.84765625 * e8,
            -3: 123948.326111417 * e18 + 46702.2873911873 * e16 + 15651.9454827029 * e14 + 4500.02279358852 * e12 + 1046.12310418023 * e10 + 175.900987413194 * e8 + 16.1671006944444 * e6,
            -2: -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1: 242310.069368053 * e18 + 101373.176843053 * e16 + 38957.2608385778 * e14 + 13491.608223504 * e12 + 4098.99631415473 * e10 + 1049.76285807292 * e8 + 212.49609375 * e6 + 30.1875 * e4 + 2.25 * e2,
            0 : 324435.666552975 * e18 + 141220.816219917 * e16 + 57049.3071968596 * e14 + 21076.3362550637 * e12 + 6983.08883680556 * e10 + 2020.14366319444 * e8 + 491.111111111111 * e6 + 94.625 * e4 + 13.0 * e2 + 1.0,
            1 : 423968.616774776 * e18 + 190675.584203681 * e16 + 80147.647066352 * e14 + 31090.9299649048 * e12 + 10944.0567199707 * e10 + 3412.37548828125 * e8 + 909.27734375 * e6 + 190.4375 * e4 + 30.25 * e2,
            2 : 541239.052023032 * e18 + 249801.02931992 * e16 + 108181.088184524 * e14 + 43423.2953993056 * e12 + 15846.3694444444 * e10 + 5180.95833333333 * e8 + 1330.0 * e6 + 361.0 * e4,
            3 : 675610.059154617 * e18 + 317908.286921933 * e16 + 140797.275858066 * e14 + 57148.0307461303 * e12 + 22338.8989342584 * e10 + 5826.87038845486 * e8 + 2767.19835069444 * e6,
            4 : 824036.876297084 * e18 + 397402.852876674 * e16 + 169084.547075893 * e14 + 84145.69125 * e12 + 15357.375 * e10 + 16256.25 * e8,
            5 : 1033970.02634054 * e18 + 396921.039788735 * e16 + 314179.867007791 * e14 - 407.969009518094 * e12 + 79813.1655143907 * e10,
            6 : 508663.889735776 * e18 + 1270114.88737564 * e16 - 286610.403578593 * e14 + 344140.297316744 * e12,
            7 : 5558336.7958343 * e18 - 2128818.75773916 * e16 + 1344523.62270713 * e14,
            8 : -11090993.4337707 * e18 + 4860787.77403088 * e16,
            9 : 16504761.4644517 * e18,
            },
        3: {
            -9: 346586.72312203 * e18,
            -8: 400508.041841002 * e18 + 127095.901675536 * e16,
            -7: 568948.030465282 * e18 + 176327.753474399 * e16 + 44792.9712759683 * e14,
            -6: 718586.528210449 * e18 + 251722.407875977 * e16 + 72003.30046875 * e14 + 15027.69515625 * e12,
            -5: 861544.089198368 * e18 + 323221.433356364 * e16 + 104830.957339307 * e14 + 27211.1487738998 * e12 + 4734.26420254178 * e10,
            -4: 990005.114006545 * e18 + 389453.523706838 * e16 + 136234.820096451 * e14 + 40533.1550231481 * e12 + 9414.44826388889 * e10 + 1372.08506944444 * e8,
            -3: 1097978.52197659 * e18 + 446353.172189103 * e16 + 164107.183336115 * e14 + 52945.9238456726 * e12 + 14262.7660675049 * e10 + 2916.67236328125 * e8 + 353.91015625 * e6,
            -2: 1180089.36365881 * e18 + 490358.387215742 * e16 + 186188.84331115 * e14 + 63143.2810980903 * e12 + 18484.9947916667 * e10 + 4417.99652777778 * e8 + 775.833333333333 * e6 + 76.5625 * e4,
            -1: 1231830.22560534 * e18 + 518442.206898254 * e16 + 200534.725930072 * e14 + 69942.9386810642 * e12 + 21415.4602681478 * e10 + 5531.57416449653 * e8 + 1130.38151041667 * e6 + 162.3125 * e4 + 12.25 * e2,
            0 : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -9: 16504761.4644517 * e18,
            },
        5: {
            -9: 392991619.772678 * e18,
            },
        6: {
            -9: 5792215177.98438 * e18,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[3][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[3][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[3][-8]
    eccentricity_results_bymode[3][9] = eccentricity_results_bymode[3][-9]
    eccentricity_results_bymode[4][-8] = eccentricity_results_bymode[2][8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[2][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[4][9] = eccentricity_results_bymode[2][-9]
    eccentricity_results_bymode[5][-8] = eccentricity_results_bymode[1][8]
    eccentricity_results_bymode[5][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[5][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[5][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[6][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[6][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[6][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[6][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[6][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[6][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[6][9] = eccentricity_results_bymode[0][-9]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc20(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^20 for order-l = 6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20.
    #     and order-l = 6.

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
            -10: 7.9629474103313e-8 * e20,
            -9 : 4.88209082003759e-8 * e20 + 1.12231972874427e-8 * e18,
            -8 : 6.09669995756076e-9 * e20 + 2.52882145656617e-9 * e18 + 6.15118732678257e-10 * e16,
            -7 : 3.34076504682396e-11 * e20 + 1.94834424227877e-11 * e18 + 8.86035283887137e-12 * e16 + 2.40280754952444e-12 * e14,
            -5 : 2.76446846101869e-7 * e20 + 2.54098263498964e-7 * e18 + 2.23275635398044e-7 * e16 + 1.82274577899374e-7 * e14 + 1.29982277199074e-7 * e12 + 6.78168402777778e-8 * e10,
            -4 : 0.000274386702816147 * e20 + 0.000314514650542573 * e18 + 0.000364889047343474 * e16 + 0.000423831569664903 * e14 + 0.000552662037037037 * e12 + 0.000347222222222222 * e10 + 0.00173611111111111 * e8,
            -3 : 0.00682729907801413 * e20 + 0.00951790638800178 * e18 - 0.00410451778343746 * e16 + 0.0818452875954764 * e14 - 0.268883514404297 * e12 + 0.696340942382813 * e10 - 0.83056640625 * e8 + 0.31640625 * e6,
            -2 : 0.933368535270708 * e20 - 4.47036443461129 * e18 + 17.8718484071869 * e16 - 51.4749018959436 * e14 + 102.849652777778 * e12 - 128.8 * e10 + 89.6111111111111 * e8 - 30.6666666666667 * e6 + 4.0 * e4,
            -1 : -504.418427161525 * e20 + 1556.52804495489 * e18 - 3669.99875953992 * e16 + 6344.85087845691 * e14 - 7628.348520067 * e12 + 6002.7229309082 * e10 - 2910.81678602431 * e8 + 808.756510416667 * e6 - 114.0625 * e4 + 6.25 * e2,
            0  : 70317.5436322466 * e20 - 141903.055105815 * e18 + 218602.008076882 * e16 - 247941.106069037 * e14 + 198843.244707031 * e12 - 108050.40984375 * e10 + 37913.58984375 * e8 - 8060.0 * e6 + 941.625 * e4 - 51.0 * e2 + 1.0,
            1  : -3571873.50712728 * e20 + 4985882.10295459 * e18 - 5312856.42342945 * e16 + 4187599.83676463 * e14 - 2357802.65863441 * e12 + 910799.197898356 * e10 - 229156.124728733 * e8 + 34797.1002604167 * e6 - 2801.3125 * e4 + 90.25 * e2,
            2  : 83589501.397737 * e20 - 83708417.4232456 * e18 + 63871121.2968106 * e16 - 36008756.8366826 * e14 + 14457660.8707465 * e12 - 3940116.19166667 * e10 + 679077.819444444 * e8 - 65473.6666666667 * e6 + 2652.25 * e4,
            3  : -1040749555.24248 * e20 + 763624009.156845 * e18 - 423808849.698558 * e16 + 171836856.188733 * e14 - 48600122.2348045 * e12 + 8950487.88320618 * e10 - 951271.984863281 * e8 + 43759.41015625 * e6,
            4  : 7541060698.92719 * e20 - 4083236400.87082 * e18 + 1645491640.43578 * e16 - 471706639.341561 * e14 + 89871728.1598126 * e12 - 10092278.1358507 * e10 + 500703.656684028 * e8,
            5  : -33546378608.8935 * e20 + 13311984041.219 * e18 - 3813865880.14918 * e16 + 737419222.764519 * e14 - 85358456.808223 * e12 + 4432418.57569994 * e10,
            6  : 94017846081.5128 * e20 - 26684504470.8569 * e18 + 5172892203.72458 * e16 - 607649410.257857 * e14 + 32405694.76 * e12,
            7  : -165975055751.793 * e20 + 32001453962.3975 * e18 - 3775523947.20958 * e16 + 204185799.366173 * e14,
            8  : 178547514320.954 * e20 - 21006631484.0241 * e18 + 1142015861.10956 * e16,
            9  : -106655280236.629 * e20 + 5792215177.98438 * e18,
            10 : 27070739641.9779 * e20,
            },
        1: {
            -10: 3.93444492265351 * e20,
            -9 : 17.9565920499232 * e20 + 2.11929181104972 * e18,
            -8 : 52.1348849140079 * e20 + 10.2664729798742 * e18 + 1.14273024524661 * e16,
            -7 : 121.326731696711 * e20 + 31.076514676145 * e18 + 5.85592993927245 * e16 + 0.616886645141913 * e14,
            -6 : 246.239746841218 * e20 + 74.6715785950487 * e18 + 18.4574599376245 * e16 + 3.33357304756393 * e14 + 0.333426046489198 * e12,
            -5 : 454.540499592287 * e20 + 155.504597129133 * e18 + 45.7630135075178 * e16 + 10.9254498859688 * e14 + 1.8942113410102 * e12 + 0.180403713650174 * e10,
            -4 : -0.09765625 * e8 / (e2 - 1.0)**11,
            -3 : 1273.63610321758 * e20 + 513.674032637722 * e18 + 188.217183966076 * e16 + 61.124635589825 * e14 + 16.9546108198755 * e12 + 3.78169793023003 * e10 + 0.605740017361111 * e8 + 0.0525173611111111 * e6,
            -2 : 2015.30858990235 * e20 + 868.670156063517 * e18 + 346.866849481492 * e16 + 126.394356708829 * e14 + 41.1977593315972 * e12 + 11.6824652777778 * e10 + 2.76302083333333 * e8 + 0.5 * e6 + 0.0625 * e4,
            -1 : 2984.22970048062 * e20 + 1335.2030446067 * e18 + 553.680295894864 * e16 + 208.170076889311 * e14 + 68.1135433959961 * e12 + 18.1928527832031 * e10 + 1.62646484375 * e8 + 1.37109375 * e6 - 1.9375 * e4 + 0.25 * e2,
            0  : 4145.58489923501 * e20 + 1885.7486840863 * e18 + 792.65418115502 * e16 + 281.988204601285 * e14 + 130.753796296296 * e12 - 68.4535069444444 * e10 + 130.514756944444 * e8 - 117.222222222222 * e6 + 52.25 * e4 - 11.0 * e2 + 1.0,
            1  : 5439.42982139739 * e20 + 2734.02074166065 * e18 + 369.822816342269 * e16 + 2109.77532191882 * e14 - 3053.98570166694 * e12 + 4141.78814358181 * e10 - 3489.35139973958 * e8 + 1769.43359375 * e6 - 473.4375 * e4 + 56.25 * e2,
            2  : 14886.0344467632 * e20 - 17290.9875345285 * e18 + 44089.1411303711 * e16 - 66859.2437695313 * e14 + 77679.9017089844 * e12 - 61175.3546875 * e10 + 30777.9140625 * e8 - 8711.5 * e6 + 1105.5625 * e4,
            3  : -393242.00429522 * e20 + 720960.855457223 * e18 - 988542.074944686 * e16 + 1022138.69151758 * e14 - 745768.239982143 * e12 + 359143.306014336 * e10 - 101014.721082899 * e8 + 12858.6150173611 * e6,
            4  : 9029384.94376862 * e20 - 11092205.6537999 * e18 + 10369096.3572048 * e16 - 7000489.52094411 * e14 + 3188086.33699653 * e12 - 868382.78828125 * e10 + 108035.47265625 * e8,
            5  : -100968875.479466 * e20 + 86189199.8191689 * e18 - 53965266.192221 * e16 + 23151140.4241943 * e14 - 6038281.67580887 * e12 + 725442.129741821 * e10,
            6  : 612087503.683281 * e20 - 356665472.701488 * e18 + 144077742.421739 * e16 - 35811174.726429 * e14 + 4131091.53653369 * e12,
            7  : -2083048317.14278 * e20 + 793257308.259555 * e18 - 187589612.523247 * e16 + 20725438.1895036 * e14,
            8  : 3952895967.90602 * e20 - 889350074.882826 * e18 + 94024315.6193818 * e16,
            9  : -3884299971.10362 * e20 + 392991619.772678 * e18,
            10 : 1534449159.58053 * e20,
            },
        2: {
            -10: 5249.91649030485 * e20,
            -9 : 15049.7296110611 * e20 + 2423.74006354441 * e18,
            -8 : 31913.6819331855 * e20 + 7597.54040902274 * e18 + 1105.2582240638 * e16,
            -7 : 57857.0096881669 * e20 + 16881.7174939228 * e18 + 3758.40170018442 * e16 + 496.705916584755 * e14,
            -6 : 95146.0730971624 * e20 + 31624.7449636542 * e18 + 8745.56848149821 * e16 + 1819.56243496473 * e14 + 219.327872299383 * e12,
            -5 : 146190.527291719 * e20 + 53311.4774134068 * e18 + 16923.8945209788 * e16 + 4428.17964847701 * e14 + 860.110561676025 * e12 + 94.7732672119141 * e10,
            -4 : 213518.74084949 * e20 + 83529.8813267346 * e18 + 29238.0507502426 * e16 + 8845.59142629795 * e14 + 2185.01003255208 * e12 + 395.504427083333 * e10 + 39.84765625 * e8,
            -3 : 299747.006601623 * e20 + 123948.326111417 * e18 + 46702.2873911873 * e16 + 15651.9454827029 * e14 + 4500.02279358852 * e12 + 1046.12310418023 * e10 + 175.900987413194 * e8 + 16.1671006944444 * e6,
            -2 : -1.5625 * e4 * (e2 + 2.0)**2 / (e2 - 1.0)**11,
            -1 : 539615.75929479 * e20 + 242310.069368053 * e18 + 101373.176843053 * e16 + 38957.2608385778 * e14 + 13491.608223504 * e12 + 4098.99631415473 * e10 + 1049.76285807292 * e8 + 212.49609375 * e6 + 30.1875 * e4 + 2.25 * e2,
            0  : 699633.201211568 * e20 + 324435.666552975 * e18 + 141220.816219917 * e16 + 57049.3071968596 * e14 + 21076.3362550637 * e12 + 6983.08883680556 * e10 + 2020.14366319444 * e8 + 491.111111111111 * e6 + 94.625 * e4 + 13.0 * e2 + 1.0,
            1  : 889725.20586972 * e20 + 423968.616774776 * e18 + 190675.584203681 * e16 + 80147.647066352 * e14 + 31090.9299649048 * e12 + 10944.0567199707 * e10 + 3412.37548828125 * e8 + 909.27734375 * e6 + 190.4375 * e4 + 30.25 * e2,
            2  : 1110684.66372005 * e20 + 541239.052023032 * e18 + 249801.02931992 * e16 + 108181.088184524 * e14 + 43423.2953993056 * e12 + 15846.3694444444 * e10 + 5180.95833333333 * e8 + 1330.0 * e6 + 361.0 * e4,
            3  : 1361956.27684728 * e20 + 675610.059154617 * e18 + 317908.286921933 * e16 + 140797.275858066 * e14 + 57148.0307461303 * e12 + 22338.8989342584 * e10 + 5826.87038845486 * e8 + 2767.19835069444 * e6,
            4  : 1642041.42656992 * e20 + 824036.876297084 * e18 + 397402.852876674 * e16 + 169084.547075893 * e14 + 84145.69125 * e12 + 15357.375 * e10 + 16256.25 * e8,
            5  : 1925232.9246357 * e20 + 1033970.02634054 * e18 + 396921.039788735 * e16 + 314179.867007791 * e14 - 407.969009518094 * e12 + 79813.1655143907 * e10,
            6  : 2719444.79183074 * e20 + 508663.889735776 * e18 + 1270114.88737564 * e16 - 286610.403578593 * e14 + 344140.297316744 * e12,
            7  : -1841041.26229933 * e20 + 5558336.7958343 * e18 - 2128818.75773916 * e16 + 1344523.62270713 * e14,
            8  : 24828601.6998815 * e20 - 11090993.4337707 * e18 + 4860787.77403088 * e16,
            9  : -48577408.7564983 * e20 + 16504761.4644517 * e18,
            10 : 53213324.1545223 * e20,
            },
        3: {
            -10: 914708.371930524 * e20,
            -9 : 837339.476451071 * e20 + 346586.72312203 * e18,
            -8 : 1224414.59118849 * e20 + 400508.041841002 * e18 + 127095.901675536 * e16,
            -7 : 1513112.10593381 * e20 + 568948.030465282 * e18 + 176327.753474399 * e16 + 44792.9712759683 * e14,
            -6 : 1800410.99445078 * e20 + 718586.528210449 * e18 + 251722.407875977 * e16 + 72003.30046875 * e14 + 15027.69515625 * e12,
            -5 : 2065418.85812685 * e20 + 861544.089198368 * e18 + 323221.433356364 * e16 + 104830.957339307 * e14 + 27211.1487738998 * e12 + 4734.26420254178 * e10,
            -4 : 2298995.60317224 * e20 + 990005.114006545 * e18 + 389453.523706838 * e16 + 136234.820096451 * e14 + 40533.1550231481 * e12 + 9414.44826388889 * e10 + 1372.08506944444 * e8,
            -3 : 2492286.27403881 * e20 + 1097978.52197659 * e18 + 446353.172189103 * e16 + 164107.183336115 * e14 + 52945.9238456726 * e12 + 14262.7660675049 * e10 + 2916.67236328125 * e8 + 353.91015625 * e6,
            -2 : 2637494.71615794 * e20 + 1180089.36365881 * e18 + 490358.387215742 * e16 + 186188.84331115 * e14 + 63143.2810980903 * e12 + 18484.9947916667 * e10 + 4417.99652777778 * e8 + 775.833333333333 * e6 + 76.5625 * e4,
            -1 : 2728144.58162382 * e20 + 1231830.22560534 * e18 + 518442.206898254 * e16 + 200534.725930072 * e14 + 69942.9386810642 * e12 + 21415.4602681478 * e10 + 5531.57416449653 * e8 + 1130.38151041667 * e6 + 162.3125 * e4 + 12.25 * e2,
            0  : -0.015625 * (15.0 * e4 + 40.0 * e2 + 8.0)**2 / (e2 - 1.0)**11,
            },
        4: {
            -10: 53213324.1545223 * e20,
            },
        5: {
            -10: 1534449159.58053 * e20,
            },
        6: {
            -10: 27070739641.9779 * e20,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[3][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[3][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[3][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[3][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[3][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[3][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[3][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[3][-8]
    eccentricity_results_bymode[3][9] = eccentricity_results_bymode[3][-9]
    eccentricity_results_bymode[3][10] = eccentricity_results_bymode[3][-10]
    eccentricity_results_bymode[4][-9] = eccentricity_results_bymode[2][9]
    eccentricity_results_bymode[4][-8] = eccentricity_results_bymode[2][8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[2][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[4][9] = eccentricity_results_bymode[2][-9]
    eccentricity_results_bymode[4][10] = eccentricity_results_bymode[2][-10]
    eccentricity_results_bymode[5][-9] = eccentricity_results_bymode[1][9]
    eccentricity_results_bymode[5][-8] = eccentricity_results_bymode[1][8]
    eccentricity_results_bymode[5][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[5][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[5][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[5][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[5][10] = eccentricity_results_bymode[1][-10]
    eccentricity_results_bymode[6][-9] = eccentricity_results_bymode[0][9]
    eccentricity_results_bymode[6][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[6][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[6][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[6][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[6][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[6][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[6][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[6][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[6][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[6][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[6][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[6][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[6][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[6][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[6][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[6][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[6][9] = eccentricity_results_bymode[0][-9]
    eccentricity_results_bymode[6][10] = eccentricity_results_bymode[0][-10]

    return eccentricity_results_bymode
