""" Eccentricity functions (squared) for various truncations of e at tidal order-l = 3
"""

from typing import Dict, TYPE_CHECKING

from . import EccenOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def eccentricity_funcs_trunc2(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^2 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 2.
    #     and order-l = 3.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -1: e2,
            0 : 1.0 - 12.0 * e2,
            1 : 25.0 * e2,
            },
        1: {
            -1: -e2 / (e2 - 1.0)**5,
            0 : 4.0 * e2 + 1.0,
            1 : 9.0 * e2,
            },
        2: {
            -1: 9.0 * e2,
            },
        3: {
            -1: 25.0 * e2,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc4(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^4 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 4.
    #     and order-l = 3.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -2: 0.015625 * e4,
            -1: -2.5 * e4 + e2,
            0 : 49.21875 * e4 - 12.0 * e2 + 1.0,
            1 : -220.0 * e4 + 25.0 * e2,
            2 : 252.015625 * e4,
            },
        1: {
            -2: 1.890625 * e4,
            -1: -e2 / (e2 - 1.0)**5,
            0 : 11.46875 * e4 + 4.0 * e2 + 1.0,
            1 : 16.5 * e4 + 9.0 * e2,
            2 : 43.890625 * e4,
            },
        2: {
            -2: 43.890625 * e4,
            },
        3: {
            -2: 252.015625 * e4,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc6(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^6 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 6.
    #     and order-l = 3.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -2: 0.00520833333333333 * e6 + 0.015625 * e4,
            -1: 1.85416666666667 * e6 - 2.5 * e4 + e2,
            0 : -83.21875 * e6 + 49.21875 * e4 - 12.0 * e2 + 1.0,
            1 : 736.916666666667 * e6 - 220.0 * e4 + 25.0 * e2,
            2 : -2027.36979166667 * e6 + 252.015625 * e4,
            3 : 1660.5625 * e6,
            },
        1: {
            -3: 3.67361111111111 * e6,
            -2: 8.421875 * e6 + 1.890625 * e4,
            -1: -e2 / (e2 - 1.0)**5,
            0 : 26.4756944444444 * e6 + 11.46875 * e4 + 4.0 * e2 + 1.0,
            1 : 38.1875 * e6 + 16.5 * e4 + 9.0 * e2,
            2 : 32.296875 * e6 + 43.890625 * e4,
            3 : 164.694444444444 * e6,
            },
        2: {
            -3: 164.694444444444 * e6,
            },
        3: {
            -3: 1660.5625 * e6,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc8(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^8 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 8.
    #     and order-l = 3.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e
    e8 = e**8

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -4: 6.78168402777778e-6 * e8,
            -2: 0.00490993923611111 * e8 + 0.00520833333333333 * e6 + 0.015625 * e4,
            -1: -0.524305555555556 * e8 + 1.85416666666667 * e6 - 2.5 * e4 + e2,
            0 : 68.0408935546875 * e8 - 83.21875 * e6 + 49.21875 * e4 - 12.0 * e2 + 1.0,
            1 : -1221.72222222222 * e8 + 736.916666666667 * e6 - 220.0 * e4 + 25.0 * e2,
            2 : 6597.1491156684 * e8 - 2027.36979166667 * e6 + 252.015625 * e4,
            3 : -13126.59375 * e8 + 1660.5625 * e6,
            4 : 8504.77816433377 * e8,
            },
        1: {
            -4: 7.18072509765625 * e8,
            -3: 14.2152777777778 * e8 + 3.67361111111111 * e6,
            -2: 23.4019368489583 * e8 + 8.421875 * e6 + 1.890625 * e4,
            -1: -e2 / (e2 - 1.0)**5,
            0 : 53.2151557074653 * e8 + 26.4756944444444 * e6 + 11.46875 * e4 + 4.0 * e2 + 1.0,
            1 : 71.4791666666667 * e8 + 38.1875 * e6 + 16.5 * e4 + 9.0 * e2,
            2 : 97.048095703125 * e8 + 32.296875 * e6 + 43.890625 * e4,
            3 : -13.3680555555556 * e8 + 164.694444444444 * e6,
            4 : 532.960510253906 * e8,
            },
        2: {
            -4: 532.960510253906 * e8,
            },
        3: {
            -4: 8504.77816433377 * e8,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[0][-4]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc10(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^10 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 10.
    #     and order-l = 3.

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
            -5: 6.94444444444444e-5 * e10,
            -4: 1.35633680555556e-5 * e10 + 6.78168402777778e-6 * e8,
            -2: 0.00393880208333333 * e10 + 0.00490993923611111 * e8 + 0.00520833333333333 * e6 + 0.015625 * e4,
            -1: 0.136197916666667 * e10 - 0.524305555555556 * e8 + 1.85416666666667 * e6 - 2.5 * e4 + e2,
            0 : -31.233955078125 * e10 + 68.0408935546875 * e8 - 83.21875 * e6 + 49.21875 * e4 - 12.0 * e2 + 1.0,
            1 : 1147.2109375 * e10 - 1221.72222222222 * e8 + 736.916666666667 * e6 - 220.0 * e4 + 25.0 * e2,
            2 : -11511.4770507813 * e10 + 6597.1491156684 * e8 - 2027.36979166667 * e6 + 252.015625 * e4,
            3 : 43691.82890625 * e10 - 13126.59375 * e8 + 1660.5625 * e6,
            4 : -68154.558710395 * e10 + 8504.77816433377 * e8,
            5 : 36828.8084027778 * e10,
            },
        1: {
            -5: 14.0312673611111 * e10,
            -4: 23.6063720703125 * e10 + 7.18072509765625 * e8,
            -3: 36.3644097222222 * e10 + 14.2152777777778 * e8 + 3.67361111111111 * e6,
            -2: 51.6582790798611 * e10 + 23.4019368489583 * e8 + 8.421875 * e6 + 1.890625 * e4,
            -1: -e2 / (e2 - 1.0)**5,
            0 : 96.8403244357639 * e10 + 53.2151557074653 * e8 + 26.4756944444444 * e6 + 11.46875 * e4 + 4.0 * e2 + 1.0,
            1 : 124.026996527778 * e10 + 71.4791666666667 * e8 + 38.1875 * e6 + 16.5 * e4 + 9.0 * e2,
            2 : 149.09169921875 * e10 + 97.048095703125 * e8 + 32.296875 * e6 + 43.890625 * e4,
            3 : 254.317795138889 * e10 - 13.3680555555556 * e8 + 164.694444444444 * e6,
            4 : -416.388549804688 * e10 + 532.960510253906 * e8,
            5 : 1567.17015625 * e10,
            },
        2: {
            -5: 1567.17015625 * e10,
            },
        3: {
            -5: 36828.8084027778 * e10,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[0][-5]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc12(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^12 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 12.
    #     and order-l = 3.

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
            -6: 0.000250282287597656 * e12,
            -5: 0.000150462962962963 * e12 + 6.94444444444444e-5 * e10,
            -4: 1.84094464337384e-5 * e12 + 1.35633680555556e-5 * e10 + 6.78168402777778e-6 * e8,
            -2: 0.00329664018419054 * e12 + 0.00393880208333333 * e10 + 0.00490993923611111 * e8 + 0.00520833333333333 * e6 + 0.015625 * e4,
            -1: 0.0147482638888889 * e12 + 0.136197916666667 * e10 - 0.524305555555556 * e8 + 1.85416666666667 * e6 - 2.5 * e4 + e2,
            0 : 9.42518173217773 * e12 - 31.233955078125 * e10 + 68.0408935546875 * e8 - 83.21875 * e6 + 49.21875 * e4 - 12.0 * e2 + 1.0,
            1 : -678.596527777778 * e12 + 1147.2109375 * e10 - 1221.72222222222 * e8 + 736.916666666667 * e6 - 220.0 * e4 + 25.0 * e2,
            2 : 12278.0423667696 * e12 - 11511.4770507813 * e10 + 6597.1491156684 * e8 - 2027.36979166667 * e6 + 252.015625 * e4,
            3 : -81727.431640625 * e12 + 43691.82890625 * e10 - 13126.59375 * e8 + 1660.5625 * e6,
            4 : 236639.213326379 * e12 - 68154.558710395 * e10 + 8504.77816433377 * e8,
            5 : -303761.039259259 * e12 + 36828.8084027778 * e10,
            6 : 141428.145432472 * e12,
            },
        1: {
            -6: 27.3566682721362 * e12,
            -5: 38.0982465277778 * e12 + 14.0312673611111 * e10,
            -4: 55.6978223419189 * e12 + 23.6063720703125 * e10 + 7.18072509765625 * e8,
            -3: 75.7906635802469 * e12 + 36.3644097222222 * e10 + 14.2152777777778 * e8 + 3.67361111111111 * e6,
            -2: 99.0341271930271 * e12 + 51.6582790798611 * e10 + 23.4019368489583 * e8 + 8.421875 * e6 + 1.890625 * e4,
            -1: -e2 / (e2 - 1.0)**5,
            0 : 163.495271295618 * e12 + 96.8403244357639 * e10 + 53.2151557074653 * e8 + 26.4756944444444 * e6 + 11.46875 * e4 + 4.0 * e2 + 1.0,
            1 : 201.607039930556 * e12 + 124.026996527778 * e10 + 71.4791666666667 * e8 + 38.1875 * e6 + 16.5 * e4 + 9.0 * e2,
            2 : 240.380771446228 * e12 + 149.09169921875 * e10 + 97.048095703125 * e8 + 32.296875 * e6 + 43.890625 * e4,
            3 : 246.567255015432 * e12 + 254.317795138889 * e10 - 13.3680555555556 * e8 + 164.694444444444 * e6,
            4 : 831.843743430244 * e12 - 416.388549804688 * e10 + 532.960510253906 * e8,
            5 : -2226.549453125 * e12 + 1567.17015625 * e10,
            6 : 4308.45518784182 * e12,
            },
        2: {
            -6: 4308.45518784182 * e12,
            },
        3: {
            -6: 141428.145432472 * e12,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[0][-6]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc14(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^14 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 14.
    #     and order-l = 3.

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
            -7: 0.000644998740236836 * e14,
            -6: 0.000536319187709263 * e14 + 0.000250282287597656 * e12,
            -5: 0.000220803020282187 * e14 + 0.000150462962962963 * e12 + 6.94444444444444e-5 * e10,
            -4: 2.13616319013861e-5 * e14 + 1.84094464337384e-5 * e12 + 1.35633680555556e-5 * e10 + 6.78168402777778e-6 * e8,
            -2: 0.00281377997558163 * e14 + 0.00329664018419054 * e12 + 0.00393880208333333 * e10 + 0.00490993923611111 * e8 + 0.00520833333333333 * e6 + 0.015625 * e4,
            -1: 0.0241298776455026 * e14 + 0.0147482638888889 * e12 + 0.136197916666667 * e10 - 0.524305555555556 * e8 + 1.85416666666667 * e6 - 2.5 * e4 + e2,
            0 : -1.8579429315061 * e14 + 9.42518173217773 * e12 - 31.233955078125 * e10 + 68.0408935546875 * e8 - 83.21875 * e6 + 49.21875 * e4 - 12.0 * e2 + 1.0,
            1 : 275.77939401455 * e14 - 678.596527777778 * e12 + 1147.2109375 * e10 - 1221.72222222222 * e8 + 736.916666666667 * e6 - 220.0 * e4 + 25.0 * e2,
            2 : -8750.68870295334 * e14 + 12278.0423667696 * e12 - 11511.4770507813 * e10 + 6597.1491156684 * e8 - 2027.36979166667 * e6 + 252.015625 * e4,
            3 : 97701.9738448661 * e14 - 81727.431640625 * e12 + 43691.82890625 * e10 - 13126.59375 * e8 + 1660.5625 * e6,
            4 : -476303.109310042 * e14 + 236639.213326379 * e12 - 68154.558710395 * e10 + 8504.77816433377 * e8,
            5 : 1108883.48523754 * e14 - 303761.039259259 * e12 + 36828.8084027778 * e10,
            6 : -1209861.33836051 * e14 + 141428.145432472 * e12,
            7 : 496174.887629126 * e14,
            },
        1: {
            -7: 53.1922769850128 * e14,
            -6: 58.8906035768082 * e14 + 27.3566682721362 * e12,
            -5: 83.769890666336 * e14 + 38.0982465277778 * e12 + 14.0312673611111 * e10,
            -4: 109.782940444946 * e14 + 55.6978223419189 * e12 + 23.6063720703125 * e10 + 7.18072509765625 * e8,
            -3: 139.202924296599 * e14 + 75.7906635802469 * e12 + 36.3644097222222 * e10 + 14.2152777777778 * e8 + 3.67361111111111 * e6,
            -2: 172.38474787697 * e14 + 99.0341271930271 * e12 + 51.6582790798611 * e10 + 23.4019368489583 * e8 + 8.421875 * e6 + 1.890625 * e4,
            -1: -e2 / (e2 - 1.0)**5,
            0 : 260.314515388356 * e14 + 163.495271295618 * e12 + 96.8403244357639 * e10 + 53.2151557074653 * e8 + 26.4756944444444 * e6 + 11.46875 * e4 + 4.0 * e2 + 1.0,
            1 : 311.5595890687 * e14 + 201.607039930556 * e12 + 124.026996527778 * e10 + 71.4791666666667 * e8 + 38.1875 * e6 + 16.5 * e4 + 9.0 * e2,
            2 : 363.242660745893 * e14 + 240.380771446228 * e12 + 149.09169921875 * e10 + 97.048095703125 * e8 + 32.296875 * e6 + 43.890625 * e4,
            3 : 424.487224633488 * e14 + 246.567255015432 * e12 + 254.317795138889 * e10 - 13.3680555555556 * e8 + 164.694444444444 * e6,
            4 : 183.313497563519 * e14 + 831.843743430244 * e12 - 416.388549804688 * e10 + 532.960510253906 * e8,
            5 : 3203.3710281808 * e14 - 2226.549453125 * e12 + 1567.17015625 * e10,
            6 : -8722.04787448412 * e14 + 4308.45518784182 * e12,
            7 : 11267.6961313067 * e14,
            },
        2: {
            -7: 11267.6961313067 * e14,
            },
        3: {
            -7: 496174.887629126 * e14,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[0][-7]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc16(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^16 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 16.
    #     and order-l = 3.

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
            -8: 0.00143218490453031 * e16,
            -7: 0.00128999748047367 * e16 + 0.000644998740236836 * e14,
            -6: 0.000807711062382679 * e16 + 0.000536319187709263 * e14 + 0.000250282287597656 * e12,
            -5: 0.00027371307319224 * e16 + 0.000220803020282187 * e14 + 0.000150462962962963 * e12 + 6.94444444444444e-5 * e10,
            -4: 2.29065267392147e-5 * e16 + 2.13616319013861e-5 * e14 + 1.84094464337384e-5 * e12 + 1.35633680555556e-5 * e10 + 6.78168402777778e-6 * e8,
            -2: 0.00244013611233844 * e16 + 0.00281377997558163 * e14 + 0.00329664018419054 * e12 + 0.00393880208333333 * e10 + 0.00490993923611111 * e8 + 0.00520833333333333 * e6 + 0.015625 * e4,
            -1: 0.0198851794768833 * e16 + 0.0241298776455026 * e14 + 0.0147482638888889 * e12 + 0.136197916666667 * e10 - 0.524305555555556 * e8 + 1.85416666666667 * e6 - 2.5 * e4 + e2,
            0 : 0.369872448363779 * e16 - 1.8579429315061 * e14 + 9.42518173217773 * e12 - 31.233955078125 * e10 + 68.0408935546875 * e8 - 83.21875 * e6 + 49.21875 * e4 - 12.0 * e2 + 1.0,
            1 : -81.3516651706979 * e16 + 275.77939401455 * e14 - 678.596527777778 * e12 + 1147.2109375 * e10 - 1221.72222222222 * e8 + 736.916666666667 * e6 - 220.0 * e4 + 25.0 * e2,
            2 : 4454.5061560165 * e16 - 8750.68870295334 * e14 + 12278.0423667696 * e12 - 11511.4770507813 * e10 + 6597.1491156684 * e8 - 2027.36979166667 * e6 + 252.015625 * e4,
            3 : -80931.3161132813 * e16 + 97701.9738448661 * e14 - 81727.431640625 * e12 + 43691.82890625 * e10 - 13126.59375 * e8 + 1660.5625 * e6,
            4 : 630833.254275488 * e16 - 476303.109310042 * e14 + 236639.213326379 * e12 - 68154.558710395 * e10 + 8504.77816433377 * e8,
            5 : -2399019.14743937 * e16 + 1108883.48523754 * e14 - 303761.039259259 * e12 + 36828.8084027778 * e10,
            6 : 4656601.00533799 * e16 - 1209861.33836051 * e14 + 141428.145432472 * e12,
            7 : -4419448.48776258 * e16 + 496174.887629126 * e14,
            8 : 1622243.68803112 * e16,
            },
        1: {
            -8: 103.143449445698 * e16,
            -7: 85.0223375717474 * e16 + 53.1922769850128 * e14,
            -6: 123.839627952332 * e16 + 58.8906035768082 * e14 + 27.3566682721362 * e12,
            -5: 156.5581383653 * e16 + 83.769890666336 * e14 + 38.0982465277778 * e12 + 14.0312673611111 * e10,
            -4: 193.491497041434 * e16 + 109.782940444946 * e14 + 55.6978223419189 * e12 + 23.6063720703125 * e10 + 7.18072509765625 * e8,
            -3: 234.328517037313 * e16 + 139.202924296599 * e14 + 75.7906635802469 * e12 + 36.3644097222222 * e10 + 14.2152777777778 * e8 + 3.67361111111111 * e6,
            -2: 279.575049269658 * e16 + 172.38474787697 * e14 + 99.0341271930271 * e12 + 51.6582790798611 * e10 + 23.4019368489583 * e8 + 8.421875 * e6 + 1.890625 * e4,
            -1: -e2 / (e2 - 1.0)**5,
            0 : 395.424782489565 * e16 + 260.314515388356 * e14 + 163.495271295618 * e12 + 96.8403244357639 * e10 + 53.2151557074653 * e8 + 26.4756944444444 * e6 + 11.46875 * e4 + 4.0 * e2 + 1.0,
            1 : 462.1278311504 * e16 + 311.5595890687 * e14 + 201.607039930556 * e12 + 124.026996527778 * e10 + 71.4791666666667 * e8 + 38.1875 * e6 + 16.5 * e4 + 9.0 * e2,
            2 : 529.637537261631 * e16 + 363.242660745893 * e14 + 240.380771446228 * e12 + 149.09169921875 * e10 + 97.048095703125 * e8 + 32.296875 * e6 + 43.890625 * e4,
            3 : 595.336160783179 * e16 + 424.487224633488 * e14 + 246.567255015432 * e12 + 254.317795138889 * e10 - 13.3680555555556 * e8 + 164.694444444444 * e6,
            4 : 777.771468742585 * e16 + 183.313497563519 * e14 + 831.843743430244 * e12 - 416.388549804688 * e10 + 532.960510253906 * e8,
            5 : -1119.99898856027 * e16 + 3203.3710281808 * e14 - 2226.549453125 * e12 + 1567.17015625 * e10,
            6 : 12506.5681938571 * e16 - 8722.04787448412 * e14 + 4308.45518784182 * e12,
            7 : -29352.2401665102 * e16 + 11267.6961313067 * e14,
            8 : 28351.9614085849 * e16,
            },
        2: {
            -8: 28351.9614085849 * e16,
            },
        3: {
            -8: 1622243.68803112 * e16,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[0][-8]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc18(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^18 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 18.
    #     and order-l = 3.

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
            -9: 0.00294209382971939 * e18,
            -8: 0.00254610649694277 * e18 + 0.00143218490453031 * e16,
            -7: 0.00194395453654713 * e18 + 0.00128999748047367 * e16 + 0.000644998740236836 * e14,
            -6: 0.0010325740551462 * e18 + 0.000807711062382679 * e16 + 0.000536319187709263 * e14 + 0.000250282287597656 * e12,
            -5: 0.000309964959009581 * e18 + 0.00027371307319224 * e16 + 0.000220803020282187 * e14 + 0.000150462962962963 * e12 + 6.94444444444444e-5 * e10,
            -4: 2.35019862219653e-5 * e18 + 2.29065267392147e-5 * e16 + 2.13616319013861e-5 * e14 + 1.84094464337384e-5 * e12 + 1.35633680555556e-5 * e10 + 6.78168402777778e-6 * e8,
            -2: 0.00214341433927743 * e18 + 0.00244013611233844 * e16 + 0.00281377997558163 * e14 + 0.00329664018419054 * e12 + 0.00393880208333333 * e10 + 0.00490993923611111 * e8 + 0.00520833333333333 * e6 + 0.015625 * e4,
            -1: 0.0175260220170405 * e18 + 0.0198851794768833 * e16 + 0.0241298776455026 * e14 + 0.0147482638888889 * e12 + 0.136197916666667 * e10 - 0.524305555555556 * e8 + 1.85416666666667 * e6 - 2.5 * e4 + e2,
            0 : 0.0260870858601161 * e18 + 0.369872448363779 * e16 - 1.8579429315061 * e14 + 9.42518173217773 * e12 - 31.233955078125 * e10 + 68.0408935546875 * e8 - 83.21875 * e6 + 49.21875 * e4 - 12.0 * e2 + 1.0,
            1 : 18.5461801783682 * e18 - 81.3516651706979 * e16 + 275.77939401455 * e14 - 678.596527777778 * e12 + 1147.2109375 * e10 - 1221.72222222222 * e8 + 736.916666666667 * e6 - 220.0 * e4 + 25.0 * e2,
            2 : -1700.47068536249 * e18 + 4454.5061560165 * e16 - 8750.68870295334 * e14 + 12278.0423667696 * e12 - 11511.4770507813 * e10 + 6597.1491156684 * e8 - 2027.36979166667 * e6 + 252.015625 * e4,
            3 : 49207.8836325509 * e18 - 80931.3161132813 * e16 + 97701.9738448661 * e14 - 81727.431640625 * e12 + 43691.82890625 * e10 - 13126.59375 * e8 + 1660.5625 * e6,
            4 : -593336.542954374 * e18 + 630833.254275488 * e16 - 476303.109310042 * e14 + 236639.213326379 * e12 - 68154.558710395 * e10 + 8504.77816433377 * e8,
            5 : 3486150.09457409 * e18 - 2399019.14743937 * e16 + 1108883.48523754 * e14 - 303761.039259259 * e12 + 36828.8084027778 * e10,
            6 : -10798182.5623592 * e18 + 4656601.00533799 * e16 - 1209861.33836051 * e14 + 141428.145432472 * e12,
            7 : 17941640.86237 * e18 - 4419448.48776258 * e16 + 496174.887629126 * e14,
            8 : -15071136.1662944 * e18 + 1622243.68803112 * e16,
            9 : 5012301.89930944 * e18,
            },
        1: {
            -9: 199.482423600843 * e18,
            -8: 108.392941165549 * e18 + 103.143449445698 * e16,
            -7: 181.40987406828 * e18 + 85.0223375717474 * e16 + 53.1922769850128 * e14,
            -6: 219.49750673412 * e18 + 123.839627952332 * e16 + 58.8906035768082 * e14 + 27.3566682721362 * e12,
            -5: 265.473238732686 * e18 + 156.5581383653 * e16 + 83.769890666336 * e14 + 38.0982465277778 * e12 + 14.0312673611111 * e10,
            -4: 315.430263128451 * e18 + 193.491497041434 * e16 + 109.782940444946 * e14 + 55.6978223419189 * e12 + 23.6063720703125 * e10 + 7.18072509765625 * e8,
            -3: 369.911178495056 * e18 + 234.328517037313 * e16 + 139.202924296599 * e14 + 75.7906635802469 * e12 + 36.3644097222222 * e10 + 14.2152777777778 * e8 + 3.67361111111111 * e6,
            -2: 429.477882236213 * e18 + 279.575049269658 * e16 + 172.38474787697 * e14 + 99.0341271930271 * e12 + 51.6582790798611 * e10 + 23.4019368489583 * e8 + 8.421875 * e6 + 1.890625 * e4,
            -1: -e2 / (e2 - 1.0)**5,
            0 : 577.946203070403 * e18 + 395.424782489565 * e16 + 260.314515388356 * e14 + 163.495271295618 * e12 + 96.8403244357639 * e10 + 53.2151557074653 * e8 + 26.4756944444444 * e6 + 11.46875 * e4 + 4.0 * e2 + 1.0,
            1 : 662.550597040598 * e18 + 462.1278311504 * e16 + 311.5595890687 * e14 + 201.607039930556 * e12 + 124.026996527778 * e10 + 71.4791666666667 * e8 + 38.1875 * e6 + 16.5 * e4 + 9.0 * e2,
            2 : 748.280486437374 * e18 + 529.637537261631 * e16 + 363.242660745893 * e14 + 240.380771446228 * e12 + 149.09169921875 * e10 + 97.048095703125 * e8 + 32.296875 * e6 + 43.890625 * e4,
            3 : 834.960173141855 * e18 + 595.336160783179 * e16 + 424.487224633488 * e14 + 246.567255015432 * e12 + 254.317795138889 * e10 - 13.3680555555556 * e8 + 164.694444444444 * e6,
            4 : 886.845959899865 * e18 + 777.771468742585 * e16 + 183.313497563519 * e14 + 831.843743430244 * e12 - 416.388549804688 * e10 + 532.960510253906 * e8,
            5 : 1937.65933276118 * e18 - 1119.99898856027 * e16 + 3203.3710281808 * e14 - 2226.549453125 * e12 + 1567.17015625 * e10,
            6 : -8683.12656763917 * e18 + 12506.5681938571 * e16 - 8722.04787448412 * e14 + 4308.45518784182 * e12,
            7 : 46104.7598727891 * e18 - 29352.2401665102 * e16 + 11267.6961313067 * e14,
            8 : -89869.025264909 * e18 + 28351.9614085849 * e16,
            9 : 69178.6722268833 * e18,
            },
        2: {
            -9: 69178.6722268833 * e18,
            },
        3: {
            -9: 5012301.89930944 * e18,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-8] = eccentricity_results_bymode[1][8]
    eccentricity_results_bymode[2][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[2][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[3][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[3][9] = eccentricity_results_bymode[0][-9]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc20(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^20 for order-l = 3
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20.
    #     and order-l = 3.

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
            -10: 0.00577876217247681 * e20,
            -9 : 0.00441314074457908 * e20 + 0.00294209382971939 * e18,
            -8 : 0.00381197671927261 * e20 + 0.00254610649694277 * e18 + 0.00143218490453031 * e16,
            -7 : 0.00251161315006112 * e20 + 0.00194395453654713 * e18 + 0.00128999748047367 * e16 + 0.000644998740236836 * e14,
            -6 : 0.0012049009028884 * e20 + 0.0010325740551462 * e18 + 0.000807711062382679 * e16 + 0.000536319187709263 * e14 + 0.000250282287597656 * e12,
            -5 : 0.000332659701317298 * e20 + 0.000309964959009581 * e18 + 0.00027371307319224 * e16 + 0.000220803020282187 * e14 + 0.000150462962962963 * e12 + 6.94444444444444e-5 * e10,
            -4 : 2.348637819608e-5 * e20 + 2.35019862219653e-5 * e18 + 2.29065267392147e-5 * e16 + 2.13616319013861e-5 * e14 + 1.84094464337384e-5 * e12 + 1.35633680555556e-5 * e10 + 6.78168402777778e-6 * e8,
            -2 : 0.00190286507274192 * e20 + 0.00214341433927743 * e18 + 0.00244013611233844 * e16 + 0.00281377997558163 * e14 + 0.00329664018419054 * e12 + 0.00393880208333333 * e10 + 0.00490993923611111 * e8 + 0.00520833333333333 * e6 + 0.015625 * e4,
            -1 : 0.015568215951714 * e20 + 0.0175260220170405 * e18 + 0.0198851794768833 * e16 + 0.0241298776455026 * e14 + 0.0147482638888889 * e12 + 0.136197916666667 * e10 - 0.524305555555556 * e8 + 1.85416666666667 * e6 - 2.5 * e4 + e2,
            0  : 0.0580052081902435 * e20 + 0.0260870858601161 * e18 + 0.369872448363779 * e16 - 1.8579429315061 * e14 + 9.42518173217773 * e12 - 31.233955078125 * e10 + 68.0408935546875 * e8 - 83.21875 * e6 + 49.21875 * e4 - 12.0 * e2 + 1.0,
            1  : -3.13635328681146 * e20 + 18.5461801783682 * e18 - 81.3516651706979 * e16 + 275.77939401455 * e14 - 678.596527777778 * e12 + 1147.2109375 * e10 - 1221.72222222222 * e8 + 736.916666666667 * e6 - 220.0 * e4 + 25.0 * e2,
            2  : 506.15613018119 * e20 - 1700.47068536249 * e18 + 4454.5061560165 * e16 - 8750.68870295334 * e14 + 12278.0423667696 * e12 - 11511.4770507813 * e10 + 6597.1491156684 * e8 - 2027.36979166667 * e6 + 252.015625 * e4,
            3  : -22921.0744634175 * e20 + 49207.8836325509 * e18 - 80931.3161132813 * e16 + 97701.9738448661 * e14 - 81727.431640625 * e12 + 43691.82890625 * e10 - 13126.59375 * e8 + 1660.5625 * e6,
            4  : 417774.266643642 * e20 - 593336.542954374 * e18 + 630833.254275488 * e16 - 476303.109310042 * e14 + 236639.213326379 * e12 - 68154.558710395 * e10 + 8504.77816433377 * e8,
            5  : -3661950.53815625 * e20 + 3486150.09457409 * e18 - 2399019.14743937 * e16 + 1108883.48523754 * e14 - 303761.039259259 * e12 + 36828.8084027778 * e10,
            6  : 17077504.0535998 * e20 - 10798182.5623592 * e18 + 4656601.00533799 * e16 - 1209861.33836051 * e14 + 141428.145432472 * e12,
            7  : -44446051.9008859 * e20 + 17941640.86237 * e18 - 4419448.48776258 * e16 + 496174.887629126 * e14,
            8  : 64487056.9576431 * e20 - 15071136.1662944 * e18 + 1622243.68803112 * e16,
            9  : -48596163.8377771 * e20 + 5012301.89930944 * e18,
            10 : 14784814.8073487 * e20,
            },
        1: {
            -10: 384.877878640526 * e20,
            -9 : 101.533325796563 * e20 + 199.482423600843 * e18,
            -8 : 268.926495555016 * e20 + 108.392941165549 * e18 + 103.143449445698 * e16,
            -7 : 301.848924331354 * e20 + 181.40987406828 * e18 + 85.0223375717474 * e16 + 53.1922769850128 * e14,
            -6 : 359.344685646951 * e20 + 219.49750673412 * e18 + 123.839627952332 * e16 + 58.8906035768082 * e14 + 27.3566682721362 * e12,
            -5 : 419.911594749673 * e20 + 265.473238732686 * e18 + 156.5581383653 * e16 + 83.769890666336 * e14 + 38.0982465277778 * e12 + 14.0312673611111 * e10,
            -4 : 485.233585844873 * e20 + 315.430263128451 * e18 + 193.491497041434 * e16 + 109.782940444946 * e14 + 55.6978223419189 * e12 + 23.6063720703125 * e10 + 7.18072509765625 * e8,
            -3 : 555.708327452235 * e20 + 369.911178495056 * e18 + 234.328517037313 * e16 + 139.202924296599 * e14 + 75.7906635802469 * e12 + 36.3644097222222 * e10 + 14.2152777777778 * e8 + 3.67361111111111 * e6,
            -2 : 631.972787820299 * e20 + 429.477882236213 * e18 + 279.575049269658 * e16 + 172.38474787697 * e14 + 99.0341271930271 * e12 + 51.6582790798611 * e10 + 23.4019368489583 * e8 + 8.421875 * e6 + 1.890625 * e4,
            -1 : -e2 / (e2 - 1.0)**5,
            0  : 817.993238180167 * e20 + 577.946203070403 * e18 + 395.424782489565 * e16 + 260.314515388356 * e14 + 163.495271295618 * e12 + 96.8403244357639 * e10 + 53.2151557074653 * e8 + 26.4756944444444 * e6 + 11.46875 * e4 + 4.0 * e2 + 1.0,
            1  : 923.055064126159 * e20 + 662.550597040598 * e18 + 462.1278311504 * e16 + 311.5595890687 * e14 + 201.607039930556 * e12 + 124.026996527778 * e10 + 71.4791666666667 * e8 + 38.1875 * e6 + 16.5 * e4 + 9.0 * e2,
            2  : 1029.59788221084 * e20 + 748.280486437374 * e18 + 529.637537261631 * e16 + 363.242660745893 * e14 + 240.380771446228 * e12 + 149.09169921875 * e10 + 97.048095703125 * e8 + 32.296875 * e6 + 43.890625 * e4,
            3  : 1136.95977696348 * e20 + 834.960173141855 * e18 + 595.336160783179 * e16 + 424.487224633488 * e14 + 246.567255015432 * e12 + 254.317795138889 * e10 - 13.3680555555556 * e8 + 164.694444444444 * e6,
            4  : 1253.08691891711 * e20 + 886.845959899865 * e18 + 777.771468742585 * e16 + 183.313497563519 * e14 + 831.843743430244 * e12 - 416.388549804688 * e10 + 532.960510253906 * e8,
            5  : 993.076754640416 * e20 + 1937.65933276118 * e18 - 1119.99898856027 * e16 + 3203.3710281808 * e14 - 2226.549453125 * e12 + 1567.17015625 * e10,
            6  : 7426.33341910896 * e20 - 8683.12656763917 * e18 + 12506.5681938571 * e16 - 8722.04787448412 * e14 + 4308.45518784182 * e12,
            7  : -42652.0074389436 * e20 + 46104.7598727891 * e18 - 29352.2401665102 * e16 + 11267.6961313067 * e14,
            8  : 158366.763809994 * e20 - 89869.025264909 * e18 + 28351.9614085849 * e16,
            9  : -257564.462565875 * e20 + 69178.6722268833 * e18,
            10 : 164611.09281725 * e20,
            },
        2: {
            -10: 164611.09281725 * e20,
            },
        3: {
            -10: 14784814.8073487 * e20,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[2][-9] = eccentricity_results_bymode[1][9]
    eccentricity_results_bymode[2][-8] = eccentricity_results_bymode[1][8]
    eccentricity_results_bymode[2][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[2][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[2][10] = eccentricity_results_bymode[1][-10]
    eccentricity_results_bymode[3][-9] = eccentricity_results_bymode[0][9]
    eccentricity_results_bymode[3][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[3][9] = eccentricity_results_bymode[0][-9]
    eccentricity_results_bymode[3][10] = eccentricity_results_bymode[0][-10]

    return eccentricity_results_bymode
