""" Eccentricity functions (squared) for various truncations of e at tidal order-l = 2
"""

from typing import Dict, TYPE_CHECKING

from . import EccenOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def eccentricity_funcs_trunc2(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^2
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20 (reduced to ^2).
    #     and order-l = 2.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -1: 0.25 * e2,
            0 : -5.0 * e2 + 1.0,
            1 : 12.25 * e2,
            },
        1: {
            -1: 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -1: 12.25 * e2,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc4(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^4
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20 (reduced to ^4).
    #     and order-l = 2.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -1: -0.0625 * e4 + 0.25 * e2,
            0 : 7.875 * e4 - 5.0 * e2 + 1.0,
            1 : -53.8125 * e4 + 12.25 * e2,
            2 : 72.25 * e4,
            },
        1: {
            -2: 5.0625 * e4,
            -1: 5.0625 * e4 + 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -2: 72.25 * e4,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc6(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^6
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20 (reduced to ^6).
    #     and order-l = 2.

    # Performance and readability improvements
    e = eccentricity
    e2 = e * e
    e4 = e * e * e * e
    e6 = e * e * e * e * e * e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -3: 0.000434027777777778 * e6,
            -1: 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0 : -4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1 : 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2 : -325.833333333333 * e6 + 72.25 * e4,
            3 : 309.906684027778 * e6,
            },
        1: {
            -3: 10.97265625 * e6,
            -2: 7.875 * e6 + 5.0625 * e4,
            -1: 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -3: 309.906684027778 * e6,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc8(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^8
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20 (reduced to ^8).
    #     and order-l = 2.

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
            -3: 0.000596788194444444 * e8 + 0.000434027777777778 * e6,
            -1: 0.00613064236111111 * e8 + 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0 : 1.25043402777778 * e8 - 4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1 : -64.76318359375 * e8 + 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2 : 580.215277777778 * e8 - 325.833333333333 * e6 + 72.25 * e4,
            3 : -1491.08208550347 * e8 + 309.906684027778 * e6,
            4 : 1109.72265625 * e8,
            },
        1: {
            -4: 23.16015625 * e8,
            -3: 10.17041015625 * e8 + 10.97265625 * e6,
            -2: 12.9765625 * e8 + 7.875 * e6 + 5.0625 * e4,
            -1: 13.86865234375 * e8 + 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -4: 1109.72265625 * e8,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[1][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[0][-4]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc10(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^10
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20 (reduced to ^10).
    #     and order-l = 2.

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
            -5: 0.0040045166015625 * e10,
            -4: 0.00243055555555556 * e10 + 0.00173611111111111 * e8,
            -3: 0.000629679361979167 * e10 + 0.000596788194444444 * e8 + 0.000434027777777778 * e6,
            -1: 0.00536905924479167 * e10 + 0.00613064236111111 * e8 + 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0 : -0.181302083333333 * e10 + 1.25043402777778 * e8 - 4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1 : 28.4081359863281 * e10 - 64.76318359375 * e8 + 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2 : -547.1625 * e10 + 580.215277777778 * e8 - 325.833333333333 * e6 + 72.25 * e4,
            3 : 2986.78270975749 * e10 - 1491.08208550347 * e8 + 309.906684027778 * e6,
            4 : -5757.64921875 * e10 + 1109.72265625 * e8,
            5 : 3536.12958502875 * e10,
            },
        1: {
            -5: 47.9664459228516 * e10,
            -4: 7.76015625 * e10 + 23.16015625 * e8,
            -3: 18.3712188720703 * e10 + 10.17041015625 * e8 + 10.97265625 * e6,
            -2: 18.7921875 * e10 + 12.9765625 * e8 + 7.875 * e6 + 5.0625 * e4,
            -1: 19.7799133300781 * e10 + 13.86865234375 * e8 + 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -5: 3536.12958502875 * e10,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[1][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[1][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[0][-5]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc12(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^12
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20 (reduced to ^12).
    #     and order-l = 2.

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
            -6: 0.00790123456790123 * e12,
            -5: 0.00500564575195313 * e12 + 0.0040045166015625 * e10,
            -4: 0.00274594907407407 * e12 + 0.00243055555555556 * e10 + 0.00173611111111111 * e8,
            -3: 0.000607874364028742 * e12 + 0.000629679361979167 * e10 + 0.000596788194444444 * e8 + 0.000434027777777778 * e6,
            -1: 0.00439822726779514 * e12 + 0.00536905924479167 * e10 + 0.00613064236111111 * e8 + 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0 : 0.046263744212963 * e12 - 0.181302083333333 * e10 + 1.25043402777778 * e8 - 4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1 : -8.03940399169922 * e12 + 28.4081359863281 * e10 - 64.76318359375 * e8 + 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2 : 320.252213541667 * e12 - 547.1625 * e10 + 580.215277777778 * e8 - 325.833333333333 * e6 + 72.25 * e4,
            3 : -3357.87871443195 * e12 + 2986.78270975749 * e10 - 1491.08208550347 * e8 + 309.906684027778 * e6,
            4 : 12888.0920507812 * e12 - 5757.64921875 * e10 + 1109.72265625 * e8,
            5 : -19815.8197198091 * e12 + 3536.12958502875 * e10,
            6 : 10383.8930574846 * e12,
            },
        1: {
            -6: 97.948134765625 * e12,
            -5: -11.2431221008301 * e12 + 47.9664459228516 * e10,
            -4: 27.6902734375 * e12 + 7.76015625 * e10 + 23.16015625 * e8,
            -3: 24.530770111084 * e12 + 18.3712188720703 * e10 + 10.17041015625 * e8 + 10.97265625 * e6,
            -2: 25.660205078125 * e12 + 18.7921875 * e10 + 12.9765625 * e8 + 7.875 * e6 + 5.0625 * e4,
            -1: 26.6969418334961 * e12 + 19.7799133300781 * e10 + 13.86865234375 * e8 + 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -6: 10383.8930574846 * e12,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[1][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[1][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[1][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[0][-6]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc14(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^14
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 14.
    #     and order-l = 2.

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
            -7: 0.0146655734223904 * e14,
            -6: 0.00790123456790123 * e14 + 0.00790123456790123 * e12,
            -5: 0.00583694049290248 * e14 + 0.00500564575195313 * e12 + 0.0040045166015625 * e10,
            -4: 0.00280891754850088 * e14 + 0.00274594907407407 * e12 + 0.00243055555555556 * e10 + 0.00173611111111111 * e8,
            -3: 0.000566547921935927 * e14 + 0.000607874364028742 * e12 + 0.000629679361979167 * e10 + 0.000596788194444444 * e8 + 0.000434027777777778 * e6,
            -1: 0.00369286426160701 * e14 + 0.00439822726779514 * e12 + 0.00536905924479167 * e10 + 0.00613064236111111 * e8 + 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0 : 0.0170164280990174 * e14 + 0.046263744212963 * e12 - 0.181302083333333 * e10 + 1.25043402777778 * e8 - 4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1 : 1.68084948539734 * e14 - 8.03940399169922 * e12 + 28.4081359863281 * e10 - 64.76318359375 * e8 + 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2 : -127.884387676367 * e14 + 320.252213541667 * e12 - 547.1625 * e10 + 580.215277777778 * e8 - 325.833333333333 * e6 + 72.25 * e4,
            3 : 2441.29427500082 * e14 - 3357.87871443195 * e12 + 2986.78270975749 * e10 - 1491.08208550347 * e8 + 309.906684027778 * e6,
            4 : -16774.5880691964 * e14 + 12888.0920507812 * e12 - 5757.64921875 * e10 + 1109.72265625 * e8,
            5 : 49132.362332862 * e14 - 19815.8197198091 * e12 + 3536.12958502875 * e10,
            6 : -62736.2370623898 * e14 + 10383.8930574846 * e12,
            7 : 28704.305901533 * e14,
            },
        1: {
            -7: 197.837228013145 * e14,
            -6: -77.3016629464286 * e14 + 97.948134765625 * e12,
            -5: 51.1812637457772 * e14 - 11.2431221008301 * e12 + 47.9664459228516 * e10,
            -4: 30.1392534722222 * e14 + 27.6902734375 * e12 + 7.76015625 * e10 + 23.16015625 * e8,
            -3: 32.4878563785553 * e14 + 24.530770111084 * e12 + 18.3712188720703 * e10 + 10.17041015625 * e8 + 10.97265625 * e6,
            -2: 33.5304290674603 * e14 + 25.660205078125 * e12 + 18.7921875 * e10 + 12.9765625 * e8 + 7.875 * e6 + 5.0625 * e4,
            -1: 34.6187963317689 * e14 + 26.6969418334961 * e12 + 19.7799133300781 * e10 + 13.86865234375 * e8 + 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -7: 28704.305901533 * e14,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[1][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[1][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[1][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[1][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[0][-7]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc16(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^16
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20 (reduced to ^16).
    #     and order-l = 2.

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
            -8: 0.0264788444674745 * e16,
            -7: 0.0100825817278934 * e16 + 0.0146655734223904 * e14,
            -6: 0.00973544973544974 * e16 + 0.00790123456790123 * e14 + 0.00790123456790123 * e12,
            -5: 0.00612018406391144 * e16 + 0.00583694049290248 * e14 + 0.00500564575195313 * e12 + 0.0040045166015625 * e10,
            -4: 0.00274395547426146 * e16 + 0.00280891754850088 * e14 + 0.00274594907407407 * e12 + 0.00243055555555556 * e10 + 0.00173611111111111 * e8,
            -3: 0.000520591959751472 * e16 + 0.000566547921935927 * e14 + 0.000607874364028742 * e12 + 0.000629679361979167 * e10 + 0.000596788194444444 * e8 + 0.000434027777777778 * e6,
            -1: 0.00315436635052836 * e16 + 0.00369286426160701 * e14 + 0.00439822726779514 * e12 + 0.00536905924479167 * e10 + 0.00613064236111111 * e8 + 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0 : 0.0164815543730513 * e16 + 0.0170164280990174 * e14 + 0.046263744212963 * e12 - 0.181302083333333 * e10 + 1.25043402777778 * e8 - 4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1 : -0.204321659037045 * e16 + 1.68084948539734 * e14 - 8.03940399169922 * e12 + 28.4081359863281 * e10 - 64.76318359375 * e8 + 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2 : 37.3756312692901 * e16 - 127.884387676367 * e14 + 320.252213541667 * e12 - 547.1625 * e10 + 580.215277777778 * e8 - 325.833333333333 * e6 + 72.25 * e4,
            3 : -1249.34950189655 * e16 + 2441.29427500082 * e14 - 3357.87871443195 * e12 + 2986.78270975749 * e10 - 1491.08208550347 * e8 + 309.906684027778 * e6,
            4 : 14510.8509811837 * e16 - 16774.5880691964 * e14 + 12888.0920507812 * e12 - 5757.64921875 * e10 + 1109.72265625 * e8,
            5 : -72539.8904461506 * e16 + 49132.362332862 * e14 - 19815.8197198091 * e12 + 3536.12958502875 * e10,
            6 : 170853.057341565 * e16 - 62736.2370623898 * e14 + 10383.8930574846 * e12,
            7 : -186401.272142464 * e16 + 28704.305901533 * e14,
            8 : 75740.821221434 * e16,
            },
        1: {
            -8: 396.129942116251 * e16,
            -7: -263.885636272629 * e16 + 197.837228013145 * e14,
            -6: 123.027789057518 * e16 - 77.3016629464286 * e14 + 97.948134765625 * e12,
            -5: 30.2102372724385 * e16 + 51.1812637457772 * e14 - 11.2431221008301 * e12 + 47.9664459228516 * e10,
            -4: 40.4940519787016 * e16 + 30.1392534722222 * e14 + 27.6902734375 * e12 + 7.76015625 * e10 + 23.16015625 * e8,
            -3: 41.3038247908013 * e16 + 32.4878563785553 * e14 + 24.530770111084 * e12 + 18.3712188720703 * e10 + 10.17041015625 * e8 + 10.97265625 * e6,
            -2: 42.4070859491257 * e16 + 33.5304290674603 * e14 + 25.660205078125 * e12 + 18.7921875 * e10 + 12.9765625 * e8 + 7.875 * e6 + 5.0625 * e4,
            -1: 43.544745513067 * e16 + 34.6187963317689 * e14 + 26.6969418334961 * e12 + 19.7799133300781 * e10 + 13.86865234375 * e8 + 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -8: 75740.821221434 * e16,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[1][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[1][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[1][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[1][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[1][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[2][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[0][-8]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc18(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^18
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20 (reduced to ^18).
    #     and order-l = 2.

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
            -9: 0.0471735687549127 * e18,
            -8: 0.00882628148915816 * e18 + 0.0264788444674745 * e16,
            -7: 0.0147817395313134 * e18 + 0.0100825817278934 * e16 + 0.0146655734223904 * e14,
            -6: 0.0102005356326344 * e18 + 0.00973544973544974 * e16 + 0.00790123456790123 * e14 + 0.00790123456790123 * e12,
            -5: 0.00613748774008483 * e18 + 0.00612018406391144 * e16 + 0.00583694049290248 * e14 + 0.00500564575195313 * e12 + 0.0040045166015625 * e10,
            -4: 0.00261827591842421 * e18 + 0.00274395547426146 * e16 + 0.00280891754850088 * e14 + 0.00274594907407407 * e12 + 0.00243055555555556 * e10 + 0.00173611111111111 * e8,
            -3: 0.000476017669862863 * e18 + 0.000520591959751472 * e16 + 0.000566547921935927 * e14 + 0.000607874364028742 * e12 + 0.000629679361979167 * e10 + 0.000596788194444444 * e8 + 0.000434027777777778 * e6,
            -1: 0.00273372545506666 * e18 + 0.00315436635052836 * e16 + 0.00369286426160701 * e14 + 0.00439822726779514 * e12 + 0.00536905924479167 * e10 + 0.00613064236111111 * e8 + 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0 : 0.0142331724345241 * e18 + 0.0164815543730513 * e16 + 0.0170164280990174 * e14 + 0.046263744212963 * e12 - 0.181302083333333 * e10 + 1.25043402777778 * e8 - 4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1 : 0.0651670998787241 * e18 - 0.204321659037045 * e16 + 1.68084948539734 * e14 - 8.03940399169922 * e12 + 28.4081359863281 * e10 - 64.76318359375 * e8 + 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2 : -8.22944731937271 * e18 + 37.3756312692901 * e16 - 127.884387676367 * e14 + 320.252213541667 * e12 - 547.1625 * e10 + 580.215277777778 * e8 - 325.833333333333 * e6 + 72.25 * e4,
            3 : 476.906338866677 * e18 - 1249.34950189655 * e16 + 2441.29427500082 * e14 - 3357.87871443195 * e12 + 2986.78270975749 * e10 - 1491.08208550347 * e8 + 309.906684027778 * e6,
            4 : -9025.59107639858 * e18 + 14510.8509811837 * e16 - 16774.5880691964 * e14 + 12888.0920507812 * e12 - 5757.64921875 * e10 + 1109.72265625 * e8,
            5 : 72569.1358528209 * e18 - 72539.8904461506 * e16 + 49132.362332862 * e14 - 19815.8197198091 * e12 + 3536.12958502875 * e10,
            6 : -281837.593544038 * e18 + 170853.057341565 * e16 - 62736.2370623898 * e14 + 10383.8930574846 * e12,
            7 : 553402.330039802 * e18 - 186401.272142464 * e16 + 28704.305901533 * e14,
            8 : -526830.350208301 * e18 + 75740.821221434 * e16,
            9 : 192609.936997079 * e18,
            },
        1: {
            -9: 787.554483571802 * e18,
            -8: -741.199753151776 * e18 + 396.129942116251 * e16,
            -7: 344.931125093217 * e18 - 263.885636272629 * e16 + 197.837228013145 * e14,
            -6: 1.03552826550542 * e18 + 123.027789057518 * e16 - 77.3016629464286 * e14 + 97.948134765625 * e12,
            -5: 51.6717954338243 * e18 + 30.2102372724385 * e16 + 51.1812637457772 * e14 - 11.2431221008301 * e12 + 47.9664459228516 * e10,
            -4: 49.954815828838 * e18 + 40.4940519787016 * e16 + 30.1392534722222 * e14 + 27.6902734375 * e12 + 7.76015625 * e10 + 23.16015625 * e8,
            -3: 51.1473705865602 * e18 + 41.3038247908013 * e16 + 32.4878563785553 * e14 + 24.530770111084 * e12 + 18.3712188720703 * e10 + 10.17041015625 * e8 + 10.97265625 * e6,
            -2: 52.2890643555057 * e18 + 42.4070859491257 * e16 + 33.5304290674603 * e14 + 25.660205078125 * e12 + 18.7921875 * e10 + 12.9765625 * e8 + 7.875 * e6 + 5.0625 * e4,
            -1: 53.474220224043 * e18 + 43.544745513067 * e16 + 34.6187963317689 * e14 + 26.6969418334961 * e12 + 19.7799133300781 * e10 + 13.86865234375 * e8 + 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0 : -1 / (e2 - 1.0)**3,
            },
        2: {
            -9: 192609.936997079 * e18
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[1][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[1][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[1][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[1][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[1][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[1][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[2][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[2][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[2][9] = eccentricity_results_bymode[0][-9]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc20(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^20
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20.
    #     and order-l = 2.

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
            -10: 0.0834975554373556 * e20,
            -9 : -0.00235867843774564 * e20 + 0.0471735687549127 * e18,
            -8 : 0.0225805701430963 * e20 + 0.00882628148915816 * e18 + 0.0264788444674745 * e16,
            -7 : 0.0147818721410268 * e20 + 0.0147817395313134 * e18 + 0.0100825817278934 * e16 + 0.0146655734223904 * e14,
            -6 : 0.0103928174836464 * e20 + 0.0102005356326344 * e18 + 0.00973544973544974 * e16 + 0.00790123456790123 * e14 + 0.00790123456790123 * e12,
            -5 : 0.00599541358490836 * e20 + 0.00613748774008483 * e18 + 0.00612018406391144 * e16 + 0.00583694049290248 * e14 + 0.00500564575195313 * e12 + 0.0040045166015625 * e10,
            -4 : 0.00246818683772672 * e20 + 0.00261827591842421 * e18 + 0.00274395547426146 * e16 + 0.00280891754850088 * e14 + 0.00274594907407407 * e12 + 0.00243055555555556 * e10 + 0.00173611111111111 * e8,
            -3 : 0.000435026152113318 * e20 + 0.000476017669862863 * e18 + 0.000520591959751472 * e16 + 0.000566547921935927 * e14 + 0.000607874364028742 * e12 + 0.000629679361979167 * e10 + 0.000596788194444444 * e8 + 0.000434027777777778 * e6,
            -1 : 0.00239823363820855 * e20 + 0.00273372545506666 * e18 + 0.00315436635052836 * e16 + 0.00369286426160701 * e14 + 0.00439822726779514 * e12 + 0.00536905924479167 * e10 + 0.00613064236111111 * e8 + 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0  : 0.0125270277012167 * e20 + 0.0142331724345241 * e18 + 0.0164815543730513 * e16 + 0.0170164280990174 * e14 + 0.046263744212963 * e12 - 0.181302083333333 * e10 + 1.25043402777778 * e8 - 4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1  : 0.0295605874348395 * e20 + 0.0651670998787241 * e18 - 0.204321659037045 * e16 + 1.68084948539734 * e14 - 8.03940399169922 * e12 + 28.4081359863281 * e10 - 64.76318359375 * e8 + 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2  : 1.52120554635196 * e20 - 8.22944731937271 * e18 + 37.3756312692901 * e16 - 127.884387676367 * e14 + 320.252213541667 * e12 - 547.1625 * e10 + 580.215277777778 * e8 - 325.833333333333 * e6 + 72.25 * e4,
            3  : -141.220232909408 * e20 + 476.906338866677 * e18 - 1249.34950189655 * e16 + 2441.29427500082 * e14 - 3357.87871443195 * e12 + 2986.78270975749 * e10 - 1491.08208550347 * e8 + 309.906684027778 * e6,
            4  : 4254.98096976673 * e20 - 9025.59107639858 * e18 + 14510.8509811837 * e16 - 16774.5880691964 * e14 + 12888.0920507812 * e12 - 5757.64921875 * e10 + 1109.72265625 * e8,
            5  : -53016.0214441785 * e20 + 72569.1358528209 * e18 - 72539.8904461506 * e16 + 49132.362332862 * e14 - 19815.8197198091 * e12 + 3536.12958502875 * e10,
            6  : 319552.923519184 * e20 - 281837.593544038 * e18 + 170853.057341565 * e16 - 62736.2370623898 * e14 + 10383.8930574846 * e12,
            7  : -1008041.84574145 * e20 + 553402.330039802 * e18 - 186401.272142464 * e16 + 28704.305901533 * e14,
            8  : 1693960.27679947 * e20 - 526830.350208301 * e18 + 75740.821221434 * e16,
            9  : -1430025.62556884 * e20 + 192609.936997079 * e18,
            10 : 475327.098089377 * e20,
            },
        1: {
            -10: 1556.52012510487 * e20,
            -9 : -1892.32611401541 * e20 + 787.554483571802 * e18,
            -8 : 997.159033551745 * e20 - 741.199753151776 * e18 + 396.129942116251 * e16,
            -7 : -145.098734249892 * e20 + 344.931125093217 * e18 - 263.885636272629 * e16 + 197.837228013145 * e14,
            -6 : 77.5418329808761 * e20 + 1.03552826550542 * e18 + 123.027789057518 * e16 - 77.3016629464286 * e14 + 97.948134765625 * e12,
            -5 : 58.9354729324144 * e20 + 51.6717954338243 * e18 + 30.2102372724385 * e16 + 51.1812637457772 * e14 - 11.2431221008301 * e12 + 47.9664459228516 * e10,
            -4 : 60.826377047433 * e20 + 49.954815828838 * e18 + 40.4940519787016 * e16 + 30.1392534722222 * e14 + 27.6902734375 * e12 + 7.76015625 * e10 + 23.16015625 * e8,
            -3 : 61.9943894349361 * e20 + 51.1473705865602 * e18 + 41.3038247908013 * e16 + 32.4878563785553 * e14 + 24.530770111084 * e12 + 18.3712188720703 * e10 + 10.17041015625 * e8 + 10.97265625 * e6,
            -2 : 63.1757855906994 * e20 + 52.2890643555057 * e18 + 42.4070859491257 * e16 + 33.5304290674603 * e14 + 25.660205078125 * e12 + 18.7921875 * e10 + 12.9765625 * e8 + 7.875 * e6 + 5.0625 * e4,
            -1 : 64.4067689087038 * e20 + 53.474220224043 * e18 + 43.544745513067 * e16 + 34.6187963317689 * e14 + 26.6969418334961 * e12 + 19.7799133300781 * e10 + 13.86865234375 * e8 + 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0  : -1 / (e2 - 1.0)**3,
            },
        2: {
            -10: 475327.098089377 * e20
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[1][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[1][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[1][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[1][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[1][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[1][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[1][10] = eccentricity_results_bymode[1][-10]
    eccentricity_results_bymode[2][-9] = eccentricity_results_bymode[0][9]
    eccentricity_results_bymode[2][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[2][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[2][9] = eccentricity_results_bymode[0][-9]
    eccentricity_results_bymode[2][10] = eccentricity_results_bymode[0][-10]

    return eccentricity_results_bymode


@njit(cacheable=True)
def eccentricity_funcs_trunc22(eccentricity: 'FloatArray') -> 'EccenOutput':
    """ Calculates the eccentricity functions (by mode) truncated to e^22
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 22.
    #     and order-l = 2.

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
    e22 = e**22

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -11: 0.147354657384016 * e22,
            -10: -0.0379534342897071 * e22 + 0.0834975554373556 * e20,
            -9 : 0.0383312049297731 * e22 - 0.00235867843774564 * e20 + 0.0471735687549127 * e18,
            -8 : 0.0191369830469475 * e22 + 0.0225805701430963 * e20 + 0.00882628148915816 * e18 + 0.0264788444674745 * e16,
            -7 : 0.0153470288705224 * e22 + 0.0147818721410268 * e20 + 0.0147817395313134 * e18 + 0.0100825817278934 * e16 + 0.0146655734223904 * e14,
            -6 : 0.0102980729441047 * e22 + 0.0103928174836464 * e20 + 0.0102005356326344 * e18 + 0.00973544973544974 * e16 + 0.00790123456790123 * e14 + 0.00790123456790123 * e12,
            -5 : 0.00576885846778482 * e22 + 0.00599541358490836 * e20 + 0.00613748774008483 * e18 + 0.00612018406391144 * e16 + 0.00583694049290248 * e14 + 0.00500564575195313 * e12 + 0.0040045166015625 * e10,
            -4 : 0.0023126098625253 * e22 + 0.00246818683772672 * e20 + 0.00261827591842421 * e18 + 0.00274395547426146 * e16 + 0.00280891754850088 * e14 + 0.00274594907407407 * e12 + 0.00243055555555556 * e10 + 0.00173611111111111 * e8,
            -3 : 0.000398196429816899 * e22 + 0.000435026152113318 * e20 + 0.000476017669862863 * e18 + 0.000520591959751472 * e16 + 0.000566547921935927 * e14 + 0.000607874364028742 * e12 + 0.000629679361979167 * e10 + 0.000596788194444444 * e8 + 0.000434027777777778 * e6,
            -1 : 0.00212579981153931 * e22 + 0.00239823363820855 * e20 + 0.00273372545506666 * e18 + 0.00315436635052836 * e16 + 0.00369286426160701 * e14 + 0.00439822726779514 * e12 + 0.00536905924479167 * e10 + 0.00613064236111111 * e8 + 0.0169270833333333 * e6 - 0.0625 * e4 + 0.25 * e2,
            0  : 0.0111233338634036 * e22 + 0.0125270277012167 * e20 + 0.0142331724345241 * e18 + 0.0164815543730513 * e16 + 0.0170164280990174 * e14 + 0.046263744212963 * e12 - 0.181302083333333 * e10 + 1.25043402777778 * e8 - 4.30555555555556 * e6 + 7.875 * e4 - 5.0 * e2 + 1.0,
            1  : 0.0289420044285592 * e22 + 0.0295605874348395 * e20 + 0.0651670998787241 * e18 - 0.204321659037045 * e16 + 1.68084948539734 * e14 - 8.03940399169922 * e12 + 28.4081359863281 * e10 - 64.76318359375 * e8 + 85.83984375 * e6 - 53.8125 * e4 + 12.25 * e2,
            2  : -0.15239674828622 * e22 + 1.52120554635196 * e20 - 8.22944731937271 * e18 + 37.3756312692901 * e16 - 127.884387676367 * e14 + 320.252213541667 * e12 - 547.1625 * e10 + 580.215277777778 * e8 - 325.833333333333 * e6 + 72.25 * e4,
            3  : 33.6461456561839 * e22 - 141.220232909408 * e20 + 476.906338866677 * e18 - 1249.34950189655 * e16 + 2441.29427500082 * e14 - 3357.87871443195 * e12 + 2986.78270975749 * e10 - 1491.08208550347 * e8 + 309.906684027778 * e6,
            4  : -1578.45613287013 * e22 + 4254.98096976673 * e20 - 9025.59107639858 * e18 + 14510.8509811837 * e16 - 16774.5880691964 * e14 + 12888.0920507812 * e12 - 5757.64921875 * e10 + 1109.72265625 * e8,
            5  : 29723.8389568398 * e22 - 53016.0214441785 * e20 + 72569.1358528209 * e18 - 72539.8904461506 * e16 + 49132.362332862 * e14 - 19815.8197198091 * e12 + 3536.12958502875 * e10,
            6  : -267726.354561033 * e22 + 319552.923519184 * e20 - 281837.593544038 * e18 + 170853.057341565 * e16 - 62736.2370623898 * e14 + 10383.8930574846 * e12,
            7  : 1275965.34396805 * e22 - 1008041.84574145 * e20 + 553402.330039802 * e18 - 186401.272142464 * e16 + 28704.305901533 * e14,
            8  : -3375305.75044233 * e22 + 1693960.27679947 * e20 - 526830.350208301 * e18 + 75740.821221434 * e16,
            9  : 4951498.21534923 * e22 - 1430025.62556884 * e20 + 192609.936997079 * e18,
            10 : -3754295.53573552 * e22 + 475327.098089377 * e20,
            11 : 1144199.37368505 * e22,
            },
        1: {
            -11: 3060.95007304803 * e22,
            -10: -4560.85485925432 * e22 + 1556.52012510487 * e20,
            -9 : 2810.89451879814 * e22 - 1892.32611401541 * e20 + 787.554483571802 * e18,
            -8 : -697.674811232914 * e22 + 997.159033551745 * e20 - 741.199753151776 * e18 + 396.129942116251 * e16,
            -7 : 171.391111449986 * e22 - 145.098734249892 * e20 + 344.931125093217 * e18 - 263.885636272629 * e16 + 197.837228013145 * e14,
            -6 : 64.0482643178002 * e22 + 77.5418329808761 * e20 + 1.03552826550542 * e18 + 123.027789057518 * e16 - 77.3016629464286 * e14 + 97.948134765625 * e12,
            -5 : 71.560544965011 * e22 + 58.9354729324144 * e20 + 51.6717954338243 * e18 + 30.2102372724385 * e16 + 51.1812637457772 * e14 - 11.2431221008301 * e12 + 47.9664459228516 * e10,
            -4 : 72.6354802480823 * e22 + 60.826377047433 * e20 + 49.954815828838 * e18 + 40.4940519787016 * e16 + 30.1392534722222 * e14 + 27.6902734375 * e12 + 7.76015625 * e10 + 23.16015625 * e8,
            -3 : 73.8469179207771 * e22 + 61.9943894349361 * e20 + 51.1473705865602 * e18 + 41.3038247908013 * e16 + 32.4878563785553 * e14 + 24.530770111084 * e12 + 18.3712188720703 * e10 + 10.17041015625 * e8 + 10.97265625 * e6,
            -2 : 75.066743458041 * e22 + 63.1757855906994 * e20 + 52.2890643555057 * e18 + 42.4070859491257 * e16 + 33.5304290674603 * e14 + 25.660205078125 * e12 + 18.7921875 * e10 + 12.9765625 * e8 + 7.875 * e6 + 5.0625 * e4,
            -1 : 76.3420270998437 * e22 + 64.4067689087038 * e20 + 53.474220224043 * e18 + 43.544745513067 * e16 + 34.6187963317689 * e14 + 26.6969418334961 * e12 + 19.7799133300781 * e10 + 13.86865234375 * e8 + 8.96484375 * e6 + 5.0625 * e4 + 2.25 * e2,
            0  : -1 / (e2 - 1.0)**3,
            },
        2: {
            -11: 1144199.37368505 * e22,
            }
        }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[1][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[1][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[1][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[1][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[1][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[1][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[1][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[1][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[1][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[1][10] = eccentricity_results_bymode[1][-10]
    eccentricity_results_bymode[1][11] = eccentricity_results_bymode[1][-11]
    eccentricity_results_bymode[2][-10] = eccentricity_results_bymode[0][10]
    eccentricity_results_bymode[2][-9] = eccentricity_results_bymode[0][9]
    eccentricity_results_bymode[2][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[2][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[2][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[2][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[2][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[2][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[2][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[2][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[2][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[2][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[2][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[2][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[2][5] = eccentricity_results_bymode[0][-5]
    eccentricity_results_bymode[2][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[2][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[2][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[2][9] = eccentricity_results_bymode[0][-9]
    eccentricity_results_bymode[2][10] = eccentricity_results_bymode[0][-10]
    eccentricity_results_bymode[2][11] = eccentricity_results_bymode[0][-11]

    return eccentricity_results_bymode
