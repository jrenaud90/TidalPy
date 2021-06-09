""" Eccentricity functions (squared) for various truncations of e at tidal order-l = 5
"""

from . import EccenOutput
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray



@njit(cacheable=True)
def eccentricity_funcs_trunc2(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^2 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 2.
    #     and order-l = 5.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -1: 4.0*e2,
            0: 1.0 - 35.0*e2,
            1: 64.0*e2,
        },
        1: {
            0: 1.0 - 3.0*e2,
            1: 36.0*e2,
        },
        2: {
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 13.0*e2 + 1.0,
            1: 16.0*e2,
        },
        3: {
            -1: 16.0*e2,
        },
        4: {
            -1: 36.0*e2,
        },
        5: {
            -1: 64.0*e2,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc4(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^4 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 4.
    #     and order-l = 5.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -2: 1.265625*e4,
            -1: -46.0*e4 + 4.0*e2,
            0: 439.21875*e4 - 35.0*e2 + 1.0,
            1: -1416.0*e4 + 64.0*e2,
            2: 1396.890625*e4,
        },
        1: {
            -3: -0.25*e**6/(e2 - 1.0)**9,
            -2: 0.140625*e4,
            0: 11.71875*e4 - 3.0*e2 + 1.0,
            1: -126.0*e4 + 36.0*e2,
            2: 489.515625*e4,
        },
        2: {
            -2: 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 85.96875*e4 + 13.0*e2 + 1.0,
            1: 116.0*e4 + 16.0*e2,
            2: 118.265625*e4,
        },
        3: {
            -2: 118.265625*e4,
        },
        4: {
            -2: 489.515625*e4,
        },
        5: {
            -2: 1396.890625*e4,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc6(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^6 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 6.
    #     and order-l = 5.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e
    e6 = e*e*e*e*e*e

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -3: 0.0277777777777778*e6,
            -2: -5.0625*e6 + 1.265625*e4,
            -1: 199.583333333333*e6 - 46.0*e4 + 4.0*e2,
            0: -2507.63888888889*e6 + 439.21875*e4 - 35.0*e2 + 1.0,
            1: 12416.25*e6 - 1416.0*e4 + 64.0*e2,
            2: -25334.0208333333*e6 + 1396.890625*e4,
            3: 17733.3611111111*e6,
        },
        1: {
            -3: -0.25*e6/(e2 - 1.0)**9,
            -2: 1.3125*e6 + 0.140625*e4,
            -1: 2.25*e6,
            0: 0.25*e6 + 11.71875*e4 - 3.0*e2 + 1.0,
            1: 296.25*e6 - 126.0*e4 + 36.0*e2,
            2: -1958.0625*e6 + 489.515625*e4,
            3: 4160.25*e6,
        },
        2: {
            -3: 38.0277777777778*e6,
            -2: 116.604166666667*e6 + 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 391.527777777778*e6 + 85.96875*e4 + 13.0*e2 + 1.0,
            1: 521.583333333333*e6 + 116.0*e4 + 16.0*e2,
            2: 581.8125*e6 + 118.265625*e4,
            3: 616.694444444444*e6,
        },
        3: {
            -3: 616.694444444444*e6,
        },
        4: {
            -3: 4160.25*e6,
        },
        5: {
            -3: 17733.3611111111*e6,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[0][-3]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc8(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^8 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 8.
    #     and order-l = 5.

    # Performance and readability improvements
    e = eccentricity
    e2 = e*e
    e4 = e*e*e*e
    e6 = e*e*e*e*e*e
    e8 = e**8

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -4: 6.78168402777778e-6*e8,
            -3: -0.0138888888888889*e8 + 0.0277777777777778*e6,
            -2: 6.980712890625*e8 - 5.0625*e6 + 1.265625*e4,
            -1: -420.513888888889*e8 + 199.583333333333*e6 - 46.0*e4 + 4.0*e2,
            0: 7688.9645046658*e8 - 2507.63888888889*e6 + 439.21875*e4 - 35.0*e2 + 1.0,
            1: -56928.5*e8 + 12416.25*e6 - 1416.0*e4 + 64.0*e2,
            2: 191421.670925564*e8 - 25334.0208333333*e6 + 1396.890625*e4,
            3: -290347.722222222*e8 + 17733.3611111111*e6,
            4: 160544.211975098*e8,
        },
        1: {
            -4: 0.46197509765625*e8,
            -3: -0.25*e6/(e2 - 1.0)**9,
            -2: 6.805908203125*e8 + 1.3125*e6 + 0.140625*e4,
            -1: 12.0*e8 + 2.25*e6,
            0: 31.0682373046875*e8 + 0.25*e6 + 11.71875*e4 - 3.0*e2 + 1.0,
            1: -209.375*e8 + 296.25*e6 - 126.0*e4 + 36.0*e2,
            2: 4245.96899414063*e8 - 1958.0625*e6 + 489.515625*e4,
            3: -18721.125*e8 + 4160.25*e6,
            4: 26597.0230102539*e8,
        },
        2: {
            -4: 101.726135253906*e8,
            -3: 292.402777777778*e8 + 38.0277777777778*e6,
            -2: 577.690131293403*e8 + 116.604166666667*e6 + 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 1394.9820827908*e8 + 391.527777777778*e6 + 85.96875*e4 + 13.0*e2 + 1.0,
            1: 1800.27777777778*e8 + 521.583333333333*e6 + 116.0*e4 + 16.0*e2,
            2: 2086.38598632813*e8 + 581.8125*e6 + 118.265625*e4,
            3: 2137.73611111111*e8 + 616.694444444444*e6,
            4: 2623.6271226671*e8,
        },
        3: {
            -4: 2623.6271226671*e8,
        },
        4: {
            -4: 26597.0230102539*e8,
        },
        5: {
            -4: 160544.211975098*e8,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[0][-4]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc10(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^10 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 10.
    #     and order-l = 5.

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
            -4: 9.49435763888889e-6*e10 + 6.78168402777778e-6*e8,
            -3: 0.003125*e10 - 0.0138888888888889*e8 + 0.0277777777777778*e6,
            -2: -4.227978515625*e10 + 6.980712890625*e8 - 5.0625*e6 + 1.265625*e4,
            -1: 483.995833333333*e10 - 420.513888888889*e8 + 199.583333333333*e6 - 46.0*e4 + 4.0*e2,
            0: -13917.3656005859*e10 + 7688.9645046658*e8 - 2507.63888888889*e6 + 439.21875*e4 - 35.0*e2 + 1.0,
            1: 155069.1*e10 - 56928.5*e8 + 12416.25*e6 - 1416.0*e4 + 64.0*e2,
            2: -800837.513460286*e10 + 191421.670925564*e8 - 25334.0208333333*e6 + 1396.890625*e4,
            3: 2043397.8*e10 - 290347.722222222*e8 + 17733.3611111111*e6,
            4: -2488891.37178955*e10 + 160544.211975098*e8,
            5: 1150166.87673611*e10,
        },
        1: {
            -5: 0.855625*e10,
            -4: 3.92625732421875*e10 + 0.46197509765625*e8,
            -3: -0.25*e6/(e2 - 1.0)**9,
            -2: 25.805029296875*e10 + 6.805908203125*e8 + 1.3125*e6 + 0.140625*e4,
            -1: 43.5625*e10 + 12.0*e8 + 2.25*e6,
            0: 77.0075244140625*e10 + 31.0682373046875*e8 + 0.25*e6 + 11.71875*e4 - 3.0*e2 + 1.0,
            1: 321.05*e10 - 209.375*e8 + 296.25*e6 - 126.0*e4 + 36.0*e2,
            2: -4464.33813476563*e10 + 4245.96899414063*e8 - 1958.0625*e6 + 489.515625*e4,
            3: 42220.490625*e10 - 18721.125*e8 + 4160.25*e6,
            4: -132985.11505127*e10 + 26597.0230102539*e8,
            5: 140212.8025*e10,
        },
        2: {
            -5: 257.870069444444*e10,
            -4: 683.779284667969*e10 + 101.726135253906*e8,
            -3: 1301.596875*e10 + 292.402777777778*e8 + 38.0277777777778*e6,
            -2: 2104.82892252604*e10 + 577.690131293403*e8 + 116.604166666667*e6 + 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 4169.41127441406*e10 + 1394.9820827908*e8 + 391.527777777778*e6 + 85.96875*e4 + 13.0*e2 + 1.0,
            1: 5212.66458333333*e10 + 1800.27777777778*e8 + 521.583333333333*e6 + 116.0*e4 + 16.0*e2,
            2: 6063.48588867187*e10 + 2086.38598632813*e8 + 581.8125*e6 + 118.265625*e4,
            3: 6619.446875*e10 + 2137.73611111111*e8 + 616.694444444444*e6,
            4: 6226.46247694227*e10 + 2623.6271226671*e8,
            5: 9756.500625*e10,
        },
        3: {
            -5: 9756.500625*e10,
        },
        4: {
            -5: 140212.8025*e10,
        },
        5: {
            -5: 1150166.87673611*e10,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[0][-4]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc12(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^12 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 12.
    #     and order-l = 5.

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
            -4: 1.11078332971644e-5*e12 + 9.49435763888889e-6*e10 + 6.78168402777778e-6*e8,
            -3: -0.000713734567901235*e12 + 0.003125*e10 - 0.0138888888888889*e8 + 0.0277777777777778*e6,
            -2: 1.52299861907959*e12 - 4.227978515625*e10 + 6.980712890625*e8 - 5.0625*e6 + 1.265625*e4,
            -1: -332.984861111111*e12 + 483.995833333333*e10 - 420.513888888889*e8 + 199.583333333333*e6 - 46.0*e4 + 4.0*e2,
            0: 15894.4846605372*e12 - 13917.3656005859*e10 + 7688.9645046658*e8 - 2507.63888888889*e6 + 439.21875*e4 - 35.0*e2 + 1.0,
            1: -270838.2628125*e12 + 155069.1*e10 - 56928.5*e8 + 12416.25*e6 - 1416.0*e4 + 64.0*e2,
            2: 2096098.77213904*e12 - 800837.513460286*e10 + 191421.670925564*e8 - 25334.0208333333*e6 + 1396.890625*e4,
            3: -8235941.68966049*e12 + 2043397.8*e10 - 290347.722222222*e8 + 17733.3611111111*e6,
            4: 16952893.2292628*e12 - 2488891.37178955*e10 + 160544.211975098*e8,
            5: -17345375.0636574*e12 + 1150166.87673611*e10,
            6: 6933474.2367738*e12,
        },
        1: {
            -6: 1.58947086334229*e12,
            -5: 6.83729166666667*e12 + 0.855625*e10,
            -4: 18.7349378204346*e12 + 3.92625732421875*e10 + 0.46197509765625*e8,
            -3: -0.25*e6/(e2 - 1.0)**9,
            -2: 79.7755010604858*e12 + 25.805029296875*e10 + 6.805908203125*e8 + 1.3125*e6 + 0.140625*e4,
            -1: 127.3984375*e12 + 43.5625*e10 + 12.0*e8 + 2.25*e6,
            0: 203.985735931396*e12 + 77.0075244140625*e10 + 31.0682373046875*e8 + 0.25*e6 + 11.71875*e4 - 3.0*e2 + 1.0,
            1: 265.9675*e12 + 321.05*e10 - 209.375*e8 + 296.25*e6 - 126.0*e4 + 36.0*e2,
            2: 4121.18733882904*e12 - 4464.33813476563*e10 + 4245.96899414063*e8 - 1958.0625*e6 + 489.515625*e4,
            3: -53539.6359375*e12 + 42220.490625*e10 - 18721.125*e8 + 4160.25*e6,
            4: 323084.552558263*e12 - 132985.11505127*e10 + 26597.0230102539*e8,
            5: -771170.41375*e12 + 140212.8025*e10,
            6: 642633.015056191*e12,
        },
        2: {
            -6: 628.482485182491*e12,
            -5: 1515.75054398148*e12 + 257.870069444444*e10,
            -4: 2763.49007514954*e12 + 683.779284667969*e10 + 101.726135253906*e8,
            -3: 4361.44135802469*e12 + 1301.596875*e10 + 292.402777777778*e8 + 38.0277777777778*e6,
            -2: 6285.66432501475*e12 + 2104.82892252604*e10 + 577.690131293403*e8 + 116.604166666667*e6 + 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 10915.6940242796*e12 + 4169.41127441406*e10 + 1394.9820827908*e8 + 391.527777777778*e6 + 85.96875*e4 + 13.0*e2 + 1.0,
            1: 13273.1937847222*e12 + 5212.66458333333*e10 + 1800.27777777778*e8 + 521.583333333333*e6 + 116.0*e4 + 16.0*e2,
            2: 15343.7658376694*e12 + 6063.48588867187*e10 + 2086.38598632813*e8 + 581.8125*e6 + 118.265625*e4,
            3: 16925.4747492284*e12 + 6619.446875*e10 + 2137.73611111111*e8 + 616.694444444444*e6,
            4: 18126.4285817323*e12 + 6226.46247694227*e10 + 2623.6271226671*e8,
            5: 14286.5690625*e12 + 9756.500625*e10,
            6: 32959.0393131402*e12,
        },
        3: {
            -6: 32959.0393131402*e12,
        },
        4: {
            -6: 642633.015056191*e12,
        },
        5: {
            -6: 6933474.2367738*e12,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[0][-6]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc14(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^14 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 14.
    #     and order-l = 5.

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
            -7: 3.93675988914084e-8*e14,
            -6: 1.48012945050705e-9*e14 + 4.7095027970679e-10*e12,
            -4: 1.19306507144235e-5*e14 + 1.11078332971644e-5*e12 + 9.49435763888889e-6*e10 + 6.78168402777778e-6*e8,
            -3: -0.000270612874779541*e14 - 0.000713734567901235*e12 + 0.003125*e10 - 0.0138888888888889*e8 + 0.0277777777777778*e6,
            -2: -0.342169876098633*e14 + 1.52299861907959*e12 - 4.227978515625*e10 + 6.980712890625*e8 - 5.0625*e6 + 1.265625*e4,
            -1: 151.515158523479*e14 - 332.984861111111*e12 + 483.995833333333*e10 - 420.513888888889*e8 + 199.583333333333*e6 - 46.0*e4 + 4.0*e2,
            0: -12153.2662728873*e14 + 15894.4846605372*e12 - 13917.3656005859*e10 + 7688.9645046658*e8 - 2507.63888888889*e6 + 439.21875*e4 - 35.0*e2 + 1.0,
            1: 320604.250189732*e14 - 270838.2628125*e12 + 155069.1*e10 - 56928.5*e8 + 12416.25*e6 - 1416.0*e4 + 64.0*e2,
            2: -3687965.39941759*e14 + 2096098.77213904*e12 - 800837.513460286*e10 + 191421.670925564*e8 - 25334.0208333333*e6 + 1396.890625*e4,
            3: 21455761.7771062*e14 - 8235941.68966049*e12 + 2043397.8*e10 - 290347.722222222*e8 + 17733.3611111111*e6,
            4: -67677073.0728521*e14 + 16952893.2292628*e12 - 2488891.37178955*e10 + 160544.211975098*e8,
            5: 116838612.276596*e14 - 17345375.0636574*e12 + 1150166.87673611*e10,
            6: -103417649.657522*e14 + 6933474.2367738*e12,
            7: 36590541.423349*e14,
        },
        1: {
            -7: 2.96024314413265*e14,
            -6: 11.8825654711042*e14 + 1.58947086334229*e12,
            -5: 31.0688678075397*e14 + 6.83729166666667*e12 + 0.855625*e10,
            -4: 66.069752713158*e14 + 18.7349378204346*e12 + 3.92625732421875*e10 + 0.46197509765625*e8,
            -3: -0.25*e6/(e2 - 1.0)**9,
            -2: 213.03739320422*e14 + 79.7755010604858*e12 + 25.805029296875*e10 + 6.805908203125*e8 + 1.3125*e6 + 0.140625*e4,
            -1: 322.906640625*e14 + 127.3984375*e12 + 43.5625*e10 + 12.0*e8 + 2.25*e6,
            0: 479.959438682089*e14 + 203.985735931396*e12 + 77.0075244140625*e10 + 31.0682373046875*e8 + 0.25*e6 + 11.71875*e4 - 3.0*e2 + 1.0,
            1: 746.14517671131*e14 + 265.9675*e12 + 321.05*e10 - 209.375*e8 + 296.25*e6 - 126.0*e4 + 36.0*e2,
            2: -763.934376035418*e14 + 4121.18733882904*e12 - 4464.33813476563*e10 + 4245.96899414063*e8 - 1958.0625*e6 + 489.515625*e4,
            3: 49300.2275223214*e14 - 53539.6359375*e12 + 42220.490625*e10 - 18721.125*e8 + 4160.25*e6,
            4: -467668.906286028*e14 + 323084.552558263*e12 - 132985.11505127*e10 + 26597.0230102539*e8,
            5: 2035839.23679315*e14 - 771170.41375*e12 + 140212.8025*e10,
            6: -3855798.09033714*e14 + 642633.015056191*e12,
            7: 2647923.19051162*e14,
        },
        2: {
            -7: 1486.30902503189*e14,
            -6: 3212.81600939702*e14 + 628.482485182491*e12,
            -5: 5598.00777743331*e14 + 1515.75054398148*e12 + 257.870069444444*e10,
            -4: 8592.93173158918*e14 + 2763.49007514954*e12 + 683.779284667969*e10 + 101.726135253906*e8,
            -3: 12173.6873135058*e14 + 4361.44135802469*e12 + 1301.596875*e10 + 292.402777777778*e8 + 38.0277777777778*e6,
            -2: 16286.6587067233*e14 + 6285.66432501475*e12 + 2104.82892252604*e10 + 577.690131293403*e8 + 116.604166666667*e6 + 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 25764.5106606622*e14 + 10915.6940242796*e12 + 4169.41127441406*e10 + 1394.9820827908*e8 + 391.527777777778*e6 + 85.96875*e4 + 13.0*e2 + 1.0,
            1: 30603.6518586723*e14 + 13273.1937847222*e12 + 5212.66458333333*e10 + 1800.27777777778*e8 + 521.583333333333*e6 + 116.0*e4 + 16.0*e2,
            2: 35037.963668093*e14 + 15343.7658376694*e12 + 6063.48588867187*e10 + 2086.38598632813*e8 + 581.8125*e6 + 118.265625*e4,
            3: 38769.9656162505*e14 + 16925.4747492284*e12 + 6619.446875*e10 + 2137.73611111111*e8 + 616.694444444444*e6,
            4: 41413.6438768511*e14 + 18126.4285817323*e12 + 6226.46247694227*e10 + 2623.6271226671*e8,
            5: 45600.1819252232*e14 + 14286.5690625*e12 + 9756.500625*e10,
            6: 21856.1744357573*e14 + 32959.0393131402*e12,
            7: 103585.25130106*e14,
        },
        3: {
            -7: 103585.25130106*e14,
        },
        4: {
            -7: 2647923.19051162*e14,
        },
        5: {
            -7: 36590541.423349*e14,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[0][-7]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc16(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^16 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 16.
    #     and order-l = 5.

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
            -8: 4.04035102347938e-7*e16,
            -7: 1.37786596119929e-7*e16 + 3.93675988914084e-8*e14,
            -6: 2.85070589430298e-9*e16 + 1.48012945050705e-9*e14 + 4.7095027970679e-10*e12,
            -4: 1.22658803244327e-5*e16 + 1.19306507144235e-5*e14 + 1.11078332971644e-5*e12 + 9.49435763888889e-6*e10 + 6.78168402777778e-6*e8,
            -3: -0.000333374669312169*e16 - 0.000270612874779541*e14 - 0.000713734567901235*e12 + 0.003125*e10 - 0.0138888888888889*e8 + 0.0277777777777778*e6,
            -2: 0.0672655477268355*e16 - 0.342169876098633*e14 + 1.52299861907959*e12 - 4.227978515625*e10 + 6.980712890625*e8 - 5.0625*e6 + 1.265625*e4,
            -1: -48.6809549152416*e16 + 151.515158523479*e14 - 332.984861111111*e12 + 483.995833333333*e10 - 420.513888888889*e8 + 199.583333333333*e6 - 46.0*e4 + 4.0*e2,
            0: 6589.41925143572*e16 - 12153.2662728873*e14 + 15894.4846605372*e12 - 13917.3656005859*e10 + 7688.9645046658*e8 - 2507.63888888889*e6 + 439.21875*e4 - 35.0*e2 + 1.0,
            1: -269624.439096381*e16 + 320604.250189732*e14 - 270838.2628125*e12 + 155069.1*e10 - 56928.5*e8 + 12416.25*e6 - 1416.0*e4 + 64.0*e2,
            2: 4590017.59108127*e16 - 3687965.39941759*e14 + 2096098.77213904*e12 - 800837.513460286*e10 + 191421.670925564*e8 - 25334.0208333333*e6 + 1396.890625*e4,
            3: -38752628.33052*e16 + 21455761.7771062*e14 - 8235941.68966049*e12 + 2043397.8*e10 - 290347.722222222*e8 + 17733.3611111111*e6,
            4: 178650977.781291*e16 - 67677073.0728521*e14 + 16952893.2292628*e12 - 2488891.37178955*e10 + 160544.211975098*e8,
            5: -469216881.940332*e16 + 116838612.276596*e14 - 17345375.0636574*e12 + 1150166.87673611*e10,
            6: 697964847.197242*e16 - 103417649.657522*e14 + 6933474.2367738*e12,
            7: -545639315.352849*e16 + 36590541.423349*e14,
            8: 173709188.475204*e16,
        },
        1: {
            -8: 5.52377184242424*e16,
            -7: 20.5906135735544*e16 + 2.96024314413265*e14,
            -6: 51.3008514210703*e16 + 11.8825654711042*e14 + 1.58947086334229*e12,
            -5: 105.267217261905*e16 + 31.0688678075397*e14 + 6.83729166666667*e12 + 0.855625*e10,
            -4: 191.765672045863*e16 + 66.069752713158*e14 + 18.7349378204346*e12 + 3.92625732421875*e10 + 0.46197509765625*e8,
            -3: -0.25*e6/(e2 - 1.0)**9,
            -2: 509.062509645342*e16 + 213.03739320422*e14 + 79.7755010604858*e12 + 25.805029296875*e10 + 6.805908203125*e8 + 1.3125*e6 + 0.140625*e4,
            -1: 737.067096664186*e16 + 322.906640625*e14 + 127.3984375*e12 + 43.5625*e10 + 12.0*e8 + 2.25*e6,
            0: 1040.75103759889*e16 + 479.959438682089*e14 + 203.985735931396*e12 + 77.0075244140625*e10 + 31.0682373046875*e8 + 0.25*e6 + 11.71875*e4 - 3.0*e2 + 1.0,
            1: 1468.32180103812*e16 + 746.14517671131*e14 + 265.9675*e12 + 321.05*e10 - 209.375*e8 + 296.25*e6 - 126.0*e4 + 36.0*e2,
            2: 2831.67690130749*e16 - 763.934376035418*e14 + 4121.18733882904*e12 - 4464.33813476563*e10 + 4245.96899414063*e8 - 1958.0625*e6 + 489.515625*e4,
            3: -26996.8343554687*e16 + 49300.2275223214*e14 - 53539.6359375*e12 + 42220.490625*e10 - 18721.125*e8 + 4160.25*e6,
            4: 474366.056637128*e16 - 467668.906286028*e14 + 323084.552558263*e12 - 132985.11505127*e10 + 26597.0230102539*e8,
            5: -3290742.72052145*e16 + 2035839.23679315*e14 - 771170.41375*e12 + 140212.8025*e10,
            6: 11057452.3430296*e16 - 3855798.09033714*e14 + 642633.015056191*e12,
            7: -17211500.7383255*e16 + 2647923.19051162*e14,
            8: 10033260.6033253*e16,
        },
        2: {
            -8: 3432.10793723475*e16,
            -7: 6537.82863400829*e16 + 1486.30902503189*e14,
            -6: 10900.8992366276*e16 + 3212.81600939702*e14 + 628.482485182491*e12,
            -5: 16252.3247433725*e16 + 5598.00777743331*e14 + 1515.75054398148*e12 + 257.870069444444*e10,
            -4: 22583.2049408161*e16 + 8592.93173158918*e14 + 2763.49007514954*e12 + 683.779284667969*e10 + 101.726135253906*e8,
            -3: 29831.8829444651*e16 + 12173.6873135058*e14 + 4361.44135802469*e12 + 1301.596875*e10 + 292.402777777778*e8 + 38.0277777777778*e6,
            -2: 37900.9833138942*e16 + 16286.6587067233*e14 + 6285.66432501475*e12 + 2104.82892252604*e10 + 577.690131293403*e8 + 116.604166666667*e6 + 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 55948.5613034558*e16 + 25764.5106606622*e14 + 10915.6940242796*e12 + 4169.41127441406*e10 + 1394.9820827908*e8 + 391.527777777778*e6 + 85.96875*e4 + 13.0*e2 + 1.0,
            1: 65167.7195738753*e16 + 30603.6518586723*e14 + 13273.1937847222*e12 + 5212.66458333333*e10 + 1800.27777777778*e8 + 521.583333333333*e6 + 116.0*e4 + 16.0*e2,
            2: 73840.9226120055*e16 + 35037.963668093*e14 + 15343.7658376694*e12 + 6063.48588867187*e10 + 2086.38598632813*e8 + 581.8125*e6 + 118.265625*e4,
            3: 81534.5894062087*e16 + 38769.9656162505*e14 + 16925.4747492284*e12 + 6619.446875*e10 + 2137.73611111111*e8 + 616.694444444444*e6,
            4: 87877.1691519353*e16 + 41413.6438768511*e14 + 18126.4285817323*e12 + 6226.46247694227*e10 + 2623.6271226671*e8,
            5: 91103.6252511161*e16 + 45600.1819252232*e14 + 14286.5690625*e12 + 9756.500625*e10,
            6: 113019.456586192*e16 + 21856.1744357573*e14 + 32959.0393131402*e12,
            7: -7173.74738964278*e16 + 103585.25130106*e14,
            8: 307710.509016398*e16,
        },
        3: {
            -8: 307710.509016398*e16,
        },
        4: {
            -8: 10033260.6033253*e16,
        },
        5: {
            -8: 173709188.475204*e16,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[2][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[5][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[5][8] = eccentricity_results_bymode[0][-8]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc18(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^18 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 18.
    #     and order-l = 5.

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
    e18 = e**18

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -9: 1.99073685258283e-6*e18,
            -8: 1.48146204194244e-6*e18 + 4.04035102347938e-7*e16,
            -7: 2.9060942792755e-7*e18 + 1.37786596119929e-7*e16 + 3.93675988914084e-8*e14,
            -6: 4.37260637420107e-9*e18 + 2.85070589430298e-9*e16 + 1.48012945050705e-9*e14 + 4.7095027970679e-10*e12,
            -4: 1.22969144527851e-5*e18 + 1.22658803244327e-5*e16 + 1.19306507144235e-5*e14 + 1.11078332971644e-5*e12 + 9.49435763888889e-6*e10 + 6.78168402777778e-6*e8,
            -3: -0.000328334135812676*e18 - 0.000333374669312169*e16 - 0.000270612874779541*e14 - 0.000713734567901235*e12 + 0.003125*e10 - 0.0138888888888889*e8 + 0.0277777777777778*e6,
            -2: -0.00186072814829495*e18 + 0.0672655477268355*e16 - 0.342169876098633*e14 + 1.52299861907959*e12 - 4.227978515625*e10 + 6.980712890625*e8 - 5.0625*e6 + 1.265625*e4,
            -1: 11.8489052884327*e18 - 48.6809549152416*e16 + 151.515158523479*e14 - 332.984861111111*e12 + 483.995833333333*e10 - 420.513888888889*e8 + 199.583333333333*e6 - 46.0*e4 + 4.0*e2,
            0: -2655.40794819775*e18 + 6589.41925143572*e16 - 12153.2662728873*e14 + 15894.4846605372*e12 - 13917.3656005859*e10 + 7688.9645046658*e8 - 2507.63888888889*e6 + 439.21875*e4 - 35.0*e2 + 1.0,
            1: 168016.913806252*e18 - 269624.439096381*e16 + 320604.250189732*e14 - 270838.2628125*e12 + 155069.1*e10 - 56928.5*e8 + 12416.25*e6 - 1416.0*e4 + 64.0*e2,
            2: -4210423.34450057*e18 + 4590017.59108127*e16 - 3687965.39941759*e14 + 2096098.77213904*e12 - 800837.513460286*e10 + 191421.670925564*e8 - 25334.0208333333*e6 + 1396.890625*e4,
            3: 50941743.7211054*e18 - 38752628.33052*e16 + 21455761.7771062*e14 - 8235941.68966049*e12 + 2043397.8*e10 - 290347.722222222*e8 + 17733.3611111111*e6,
            4: -334189030.336258*e18 + 178650977.781291*e16 - 67677073.0728521*e14 + 16952893.2292628*e12 - 2488891.37178955*e10 + 160544.211975098*e8,
            5: 1267314609.88367*e18 - 469216881.940332*e16 + 116838612.276596*e14 - 17345375.0636574*e12 + 1150166.87673611*e10,
            6: -2845886807.44977*e18 + 697964847.197242*e16 - 103417649.657522*e14 + 6933474.2367738*e12,
            7: 3720318545.96604*e18 - 545639315.352849*e16 + 36590541.423349*e14,
            8: -2608495171.35896*e18 + 173709188.475204*e16,
            9: 756475110.403021*e18,
        },
        1: {
            -9: 10.3215904516103*e18,
            -8: 35.5339163313267*e18 + 5.52377184242424*e16,
            -7: 84.2818562290367*e18 + 20.5906135735544*e16 + 2.96024314413265*e14,
            -6: 166.820565000049*e18 + 51.3008514210703*e16 + 11.8825654711042*e14 + 1.58947086334229*e12,
            -5: 295.461157732805*e18 + 105.267217261905*e16 + 31.0688678075397*e14 + 6.83729166666667*e12 + 0.855625*e10,
            -4: 484.61432480672*e18 + 191.765672045863*e16 + 66.069752713158*e14 + 18.7349378204346*e12 + 3.92625732421875*e10 + 0.46197509765625*e8,
            -3: -0.25*e6/(e2 - 1.0)**9,
            -2: 1114.28187382151*e18 + 509.062509645342*e16 + 213.03739320422*e14 + 79.7755010604858*e12 + 25.805029296875*e10 + 6.805908203125*e8 + 1.3125*e6 + 0.140625*e4,
            -1: 1551.17784935361*e18 + 737.067096664186*e16 + 322.906640625*e14 + 127.3984375*e12 + 43.5625*e10 + 12.0*e8 + 2.25*e6,
            0: 2106.07947314963*e18 + 1040.75103759889*e16 + 479.959438682089*e14 + 203.985735931396*e12 + 77.0075244140625*e10 + 31.0682373046875*e8 + 0.25*e6 + 11.71875*e4 - 3.0*e2 + 1.0,
            1: 2848.63945591087*e18 + 1468.32180103812*e16 + 746.14517671131*e14 + 265.9675*e12 + 321.05*e10 - 209.375*e8 + 296.25*e6 - 126.0*e4 + 36.0*e2,
            2: 3641.36552870558*e18 + 2831.67690130749*e16 - 763.934376035418*e14 + 4121.18733882904*e12 - 4464.33813476563*e10 + 4245.96899414063*e8 - 1958.0625*e6 + 489.515625*e4,
            3: 19573.9413752093*e18 - 26996.8343554687*e16 + 49300.2275223214*e14 - 53539.6359375*e12 + 42220.490625*e10 - 18721.125*e8 + 4160.25*e6,
            4: -337574.379257527*e18 + 474366.056637128*e16 - 467668.906286028*e14 + 323084.552558263*e12 - 132985.11505127*e10 + 26597.0230102539*e8,
            5: 3720517.10718485*e18 - 3290742.72052145*e16 + 2035839.23679315*e14 - 771170.41375*e12 + 140212.8025*e10,
            6: -19711492.5183945*e18 + 11057452.3430296*e16 - 3855798.09033714*e14 + 642633.015056191*e12,
            7: 53437190.1369566*e18 - 17211500.7383255*e16 + 2647923.19051162*e14,
            8: -70232824.2232771*e18 + 10033260.6033253*e16,
            9: 35525471.4447657*e18,
        },
        2: {
            -9: 7772.81056512025*e18,
            -8: 12777.1761532472*e18 + 3432.10793723475*e16,
            -7: 20507.3945790318*e18 + 6537.82863400829*e16 + 1486.30902503189*e14,
            -6: 29690.0716470885*e18 + 10900.8992366276*e16 + 3212.81600939702*e14 + 628.482485182491*e12,
            -5: 40430.6329508264*e18 + 16252.3247433725*e16 + 5598.00777743331*e14 + 1515.75054398148*e12 + 257.870069444444*e10,
            -4: 52652.3725199578*e18 + 22583.2049408161*e16 + 8592.93173158918*e14 + 2763.49007514954*e12 + 683.779284667969*e10 + 101.726135253906*e8,
            -3: 66242.9428388668*e18 + 29831.8829444651*e16 + 12173.6873135058*e14 + 4361.44135802469*e12 + 1301.596875*e10 + 292.402777777778*e8 + 38.0277777777778*e6,
            -2: 81043.7529962832*e18 + 37900.9833138942*e16 + 16286.6587067233*e14 + 6285.66432501475*e12 + 2104.82892252604*e10 + 577.690131293403*e8 + 116.604166666667*e6 + 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 113445.440020706*e18 + 55948.5613034558*e16 + 25764.5106606622*e14 + 10915.6940242796*e12 + 4169.41127441406*e10 + 1394.9820827908*e8 + 391.527777777778*e6 + 85.96875*e4 + 13.0*e2 + 1.0,
            1: 129987.850262358*e18 + 65167.7195738753*e16 + 30603.6518586723*e14 + 13273.1937847222*e12 + 5212.66458333333*e10 + 1800.27777777778*e8 + 521.583333333333*e6 + 116.0*e4 + 16.0*e2,
            2: 145821.104084816*e18 + 73840.9226120055*e16 + 35037.963668093*e14 + 15343.7658376694*e12 + 6063.48588867187*e10 + 2086.38598632813*e8 + 581.8125*e6 + 118.265625*e4,
            3: 160339.529118675*e18 + 81534.5894062087*e16 + 38769.9656162505*e14 + 16925.4747492284*e12 + 6619.446875*e10 + 2137.73611111111*e8 + 616.694444444444*e6,
            4: 172964.719689798*e18 + 87877.1691519353*e16 + 41413.6438768511*e14 + 18126.4285817323*e12 + 6226.46247694227*e10 + 2623.6271226671*e8,
            5: 183681.21565624*e18 + 91103.6252511161*e16 + 45600.1819252232*e14 + 14286.5690625*e12 + 9756.500625*e10,
            6: 179230.245756046*e18 + 113019.456586192*e16 + 21856.1744357573*e14 + 32959.0393131402*e12,
            7: 299171.214392795*e18 - 7173.74738964278*e16 + 103585.25130106*e14,
            8: -232270.021336349*e18 + 307710.509016398*e16,
            9: 873615.635049072*e18,
        },
        3: {
            -9: 873615.635049072*e18,
        },
        4: {
            -9: 35525471.4447657*e18,
        },
        5: {
            -9: 756475110.403021*e18,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-8] = eccentricity_results_bymode[2][8]
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[2][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[3][9] = eccentricity_results_bymode[2][-9]
    eccentricity_results_bymode[4][-8] = eccentricity_results_bymode[1][8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[4][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[5][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[5][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[5][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[5][9] = eccentricity_results_bymode[0][-9]

    return eccentricity_results_bymode

@njit(cacheable=True)
def eccentricity_funcs_trunc20(eccentricity: FloatArray) -> EccenOutput:
    """ Calculates the eccentricity functions (by mode) truncated to e^20 for order-l = 5
    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    eccentricity_results_bymode : EccenOutput
    """
    # Eccentricity functions calculated at truncation level 20.
    #     and order-l = 5.

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
    e18 = e**18
    e20 = e**20

    # Unique results.
    eccentricity_results_bymode = {
        0: {
            -10: 6.90675590533522e-6*e20,
            -9: 7.36572635455645e-6*e20 + 1.99073685258283e-6*e18,
            -8: 3.2746483885436e-6*e20 + 1.48146204194244e-6*e18 + 4.04035102347938e-7*e16,
            -7: 4.81141316265875e-7*e20 + 2.9060942792755e-7*e18 + 1.37786596119929e-7*e16 + 3.93675988914084e-8*e14,
            -6: 5.88621987412319e-9*e20 + 4.37260637420107e-9*e18 + 2.85070589430298e-9*e16 + 1.48012945050705e-9*e14 + 4.7095027970679e-10*e12,
            -4: 1.21393115087935e-5*e20 + 1.22969144527851e-5*e18 + 1.22658803244327e-5*e16 + 1.19306507144235e-5*e14 + 1.11078332971644e-5*e12 + 9.49435763888889e-6*e10 + 6.78168402777778e-6*e8,
            -3: -0.000317243018794238*e20 - 0.000328334135812676*e18 - 0.000333374669312169*e16 - 0.000270612874779541*e14 - 0.000713734567901235*e12 + 0.003125*e10 - 0.0138888888888889*e8 + 0.0277777777777778*e6,
            -2: 0.00490398866237065*e20 - 0.00186072814829495*e18 + 0.0672655477268355*e16 - 0.342169876098633*e14 + 1.52299861907959*e12 - 4.227978515625*e10 + 6.980712890625*e8 - 5.0625*e6 + 1.265625*e4,
            -1: -2.14203790926097*e20 + 11.8489052884327*e18 - 48.6809549152416*e16 + 151.515158523479*e14 - 332.984861111111*e12 + 483.995833333333*e10 - 420.513888888889*e8 + 199.583333333333*e6 - 46.0*e4 + 4.0*e2,
            0: 827.553105287544*e20 - 2655.40794819775*e18 + 6589.41925143572*e16 - 12153.2662728873*e14 + 15894.4846605372*e12 - 13917.3656005859*e10 + 7688.9645046658*e8 - 2507.63888888889*e6 + 439.21875*e4 - 35.0*e2 + 1.0,
            1: -80410.9710603725*e20 + 168016.913806252*e18 - 269624.439096381*e16 + 320604.250189732*e14 - 270838.2628125*e12 + 155069.1*e10 - 56928.5*e8 + 12416.25*e6 - 1416.0*e4 + 64.0*e2,
            2: 2947564.74626804*e20 - 4210423.34450057*e18 + 4590017.59108127*e16 - 3687965.39941759*e14 + 2096098.77213904*e12 - 800837.513460286*e10 + 191421.670925564*e8 - 25334.0208333333*e6 + 1396.890625*e4,
            3: -50605644.6607078*e20 + 50941743.7211054*e18 - 38752628.33052*e16 + 21455761.7771062*e14 - 8235941.68966049*e12 + 2043397.8*e10 - 290347.722222222*e8 + 17733.3611111111*e6,
            4: 464319803.322212*e20 - 334189030.336258*e18 + 178650977.781291*e16 - 67677073.0728521*e14 + 16952893.2292628*e12 - 2488891.37178955*e10 + 160544.211975098*e8,
            5: -2465484232.85409*e20 + 1267314609.88367*e18 - 469216881.940332*e16 + 116838612.276596*e14 - 17345375.0636574*e12 + 1150166.87673611*e10,
            6: 7907260944.0293*e20 - 2845886807.44977*e18 + 697964847.197242*e16 - 103417649.657522*e14 + 6933474.2367738*e12,
            7: -15488694242.1656*e20 + 3720318545.96604*e18 - 545639315.352849*e16 + 36590541.423349*e14,
            8: 18065549483.8969*e20 - 2608495171.35896*e18 + 173709188.475204*e16,
            9: -11496079833.088*e20 + 756475110.403021*e18,
            10: 3066159585.39932*e20,
        },
        1: {
            -10: 19.3056071469455*e20,
            -9: 60.988530343192*e20 + 10.3215904516103*e18,
            -8: 137.644779547591*e20 + 35.5339163313267*e18 + 5.52377184242424*e16,
            -7: 262.811285443633*e20 + 84.2818562290367*e18 + 20.5906135735544*e16 + 2.96024314413265*e14,
            -6: 452.588529654399*e20 + 166.820565000049*e18 + 51.3008514210703*e16 + 11.8825654711042*e14 + 1.58947086334229*e12,
            -5: 725.622863280938*e20 + 295.461157732805*e18 + 105.267217261905*e16 + 31.0688678075397*e14 + 6.83729166666667*e12 + 0.855625*e10,
            -4: 1103.12373977024*e20 + 484.61432480672*e18 + 191.765672045863*e16 + 66.069752713158*e14 + 18.7349378204346*e12 + 3.92625732421875*e10 + 0.46197509765625*e8,
            -3: -0.25*e6/(e2 - 1.0)**9,
            -2: 2271.43691735276*e20 + 1114.28187382151*e18 + 509.062509645342*e16 + 213.03739320422*e14 + 79.7755010604858*e12 + 25.805029296875*e10 + 6.805908203125*e8 + 1.3125*e6 + 0.140625*e4,
            -1: 3057.66822882589*e20 + 1551.17784935361*e18 + 737.067096664186*e16 + 322.906640625*e14 + 127.3984375*e12 + 43.5625*e10 + 12.0*e8 + 2.25*e6,
            0: 4023.09978256634*e20 + 2106.07947314963*e18 + 1040.75103759889*e16 + 479.959438682089*e14 + 203.985735931396*e12 + 77.0075244140625*e10 + 31.0682373046875*e8 + 0.25*e6 + 11.71875*e4 - 3.0*e2 + 1.0,
            1: 5254.32881435914*e20 + 2848.63945591087*e18 + 1468.32180103812*e16 + 746.14517671131*e14 + 265.9675*e12 + 321.05*e10 - 209.375*e8 + 296.25*e6 - 126.0*e4 + 36.0*e2,
            2: 6914.39061146225*e20 + 3641.36552870558*e18 + 2831.67690130749*e16 - 763.934376035418*e14 + 4121.18733882904*e12 - 4464.33813476563*e10 + 4245.96899414063*e8 - 1958.0625*e6 + 489.515625*e4,
            3: 3590.76369825215*e20 + 19573.9413752093*e18 - 26996.8343554687*e16 + 49300.2275223214*e14 - 53539.6359375*e12 + 42220.490625*e10 - 18721.125*e8 + 4160.25*e6,
            4: 206050.425892116*e20 - 337574.379257527*e18 + 474366.056637128*e16 - 467668.906286028*e14 + 323084.552558263*e12 - 132985.11505127*e10 + 26597.0230102539*e8,
            5: -3086886.66758528*e20 + 3720517.10718485*e18 - 3290742.72052145*e16 + 2035839.23679315*e14 - 771170.41375*e12 + 140212.8025*e10,
            6: 24707723.4716718*e20 - 19711492.5183945*e18 + 11057452.3430296*e16 - 3855798.09033714*e14 + 642633.015056191*e12,
            7: -104125036.448503*e20 + 53437190.1369566*e18 - 17211500.7383255*e16 + 2647923.19051162*e14,
            8: 235064478.987245*e20 - 70232824.2232771*e18 + 10033260.6033253*e16,
            9: -266441035.835743*e20 + 35525471.4447657*e18,
            10: 118939066.613989*e20,
        },
        2: {
            -10: 17321.4943023806*e20,
            -9: 23899.2020714321*e20 + 7772.81056512025*e18,
            -8: 37418.6514182396*e20 + 12777.1761532472*e18 + 3432.10793723475*e16,
            -7: 52606.4011111488*e20 + 20507.3945790318*e18 + 6537.82863400829*e16 + 1486.30902503189*e14,
            -6: 70207.1738967631*e20 + 29690.0716470885*e18 + 10900.8992366276*e16 + 3212.81600939702*e14 + 628.482485182491*e12,
            -5: 90083.8466750485*e20 + 40430.6329508264*e18 + 16252.3247433725*e16 + 5598.00777743331*e14 + 1515.75054398148*e12 + 257.870069444444*e10,
            -4: 112107.5428384*e20 + 52652.3725199578*e18 + 22583.2049408161*e16 + 8592.93173158918*e14 + 2763.49007514954*e12 + 683.779284667969*e10 + 101.726135253906*e8,
            -3: 136095.104403583*e20 + 66242.9428388668*e18 + 29831.8829444651*e16 + 12173.6873135058*e14 + 4361.44135802469*e12 + 1301.596875*e10 + 292.402777777778*e8 + 38.0277777777778*e6,
            -2: 161805.953637066*e20 + 81043.7529962832*e18 + 37900.9833138942*e16 + 16286.6587067233*e14 + 6285.66432501475*e12 + 2104.82892252604*e10 + 577.690131293403*e8 + 116.604166666667*e6 + 13.140625*e4,
            -1: -0.25*e2*(3.0*e2 + 4.0)**2/(e2 - 1.0)**9,
            0: 217203.816618063*e20 + 113445.440020706*e18 + 55948.5613034558*e16 + 25764.5106606622*e14 + 10915.6940242796*e12 + 4169.41127441406*e10 + 1394.9820827908*e8 + 391.527777777778*e6 + 85.96875*e4 + 13.0*e2 + 1.0,
            1: 245457.759619914*e20 + 129987.850262358*e18 + 65167.7195738753*e16 + 30603.6518586723*e14 + 13273.1937847222*e12 + 5212.66458333333*e10 + 1800.27777777778*e8 + 521.583333333333*e6 + 116.0*e4 + 16.0*e2,
            2: 272819.779203866*e20 + 145821.104084816*e18 + 73840.9226120055*e16 + 35037.963668093*e14 + 15343.7658376694*e12 + 6063.48588867187*e10 + 2086.38598632813*e8 + 581.8125*e6 + 118.265625*e4,
            3: 298468.49874363*e20 + 160339.529118675*e18 + 81534.5894062087*e16 + 38769.9656162505*e14 + 16925.4747492284*e12 + 6619.446875*e10 + 2137.73611111111*e8 + 616.694444444444*e6,
            4: 321626.796457893*e20 + 172964.719689798*e18 + 87877.1691519353*e16 + 41413.6438768511*e14 + 18126.4285817323*e12 + 6226.46247694227*e10 + 2623.6271226671*e8,
            5: 341395.52988588*e20 + 183681.21565624*e18 + 91103.6252511161*e16 + 45600.1819252232*e14 + 14286.5690625*e12 + 9756.500625*e10,
            6: 362872.734287827*e20 + 179230.245756046*e18 + 113019.456586192*e16 + 21856.1744357573*e14 + 32959.0393131402*e12,
            7: 291645.290021157*e20 + 299171.214392795*e18 - 7173.74738964278*e16 + 103585.25130106*e14,
            8: 891387.133750462*e20 - 232270.021336349*e18 + 307710.509016398*e16,
            9: -1229182.34681983*e20 + 873615.635049072*e18,
            10: 2389663.95823154*e20,
        },
        3: {
            -10: 2389663.95823154*e20,
        },
        4: {
            -10: 118939066.613989*e20,
        },
        5: {
            -10: 3066159585.39932*e20,
        }
    }

    # Duplicate results are stored as dictionary lookups to previous calculations
    #    Generally leads to a 30--50% speed-up when working with large arrays.
    eccentricity_results_bymode[3][-9] = eccentricity_results_bymode[2][9]
    eccentricity_results_bymode[3][-8] = eccentricity_results_bymode[2][8]
    eccentricity_results_bymode[3][-7] = eccentricity_results_bymode[2][7]
    eccentricity_results_bymode[3][-6] = eccentricity_results_bymode[2][6]
    eccentricity_results_bymode[3][-5] = eccentricity_results_bymode[2][5]
    eccentricity_results_bymode[3][-4] = eccentricity_results_bymode[2][4]
    eccentricity_results_bymode[3][-3] = eccentricity_results_bymode[2][3]
    eccentricity_results_bymode[3][-2] = eccentricity_results_bymode[2][2]
    eccentricity_results_bymode[3][-1] = eccentricity_results_bymode[2][1]
    eccentricity_results_bymode[3][0] = eccentricity_results_bymode[2][0]
    eccentricity_results_bymode[3][1] = eccentricity_results_bymode[2][-1]
    eccentricity_results_bymode[3][2] = eccentricity_results_bymode[2][-2]
    eccentricity_results_bymode[3][3] = eccentricity_results_bymode[2][-3]
    eccentricity_results_bymode[3][4] = eccentricity_results_bymode[2][-4]
    eccentricity_results_bymode[3][5] = eccentricity_results_bymode[2][-5]
    eccentricity_results_bymode[3][6] = eccentricity_results_bymode[2][-6]
    eccentricity_results_bymode[3][7] = eccentricity_results_bymode[2][-7]
    eccentricity_results_bymode[3][8] = eccentricity_results_bymode[2][-8]
    eccentricity_results_bymode[3][9] = eccentricity_results_bymode[2][-9]
    eccentricity_results_bymode[3][10] = eccentricity_results_bymode[2][-10]
    eccentricity_results_bymode[4][-9] = eccentricity_results_bymode[1][9]
    eccentricity_results_bymode[4][-8] = eccentricity_results_bymode[1][8]
    eccentricity_results_bymode[4][-7] = eccentricity_results_bymode[1][7]
    eccentricity_results_bymode[4][-6] = eccentricity_results_bymode[1][6]
    eccentricity_results_bymode[4][-5] = eccentricity_results_bymode[1][5]
    eccentricity_results_bymode[4][-4] = eccentricity_results_bymode[1][4]
    eccentricity_results_bymode[4][-3] = eccentricity_results_bymode[1][3]
    eccentricity_results_bymode[4][-2] = eccentricity_results_bymode[1][2]
    eccentricity_results_bymode[4][-1] = eccentricity_results_bymode[1][1]
    eccentricity_results_bymode[4][0] = eccentricity_results_bymode[1][0]
    eccentricity_results_bymode[4][1] = eccentricity_results_bymode[1][-1]
    eccentricity_results_bymode[4][2] = eccentricity_results_bymode[1][-2]
    eccentricity_results_bymode[4][3] = eccentricity_results_bymode[1][-3]
    eccentricity_results_bymode[4][4] = eccentricity_results_bymode[1][-4]
    eccentricity_results_bymode[4][5] = eccentricity_results_bymode[1][-5]
    eccentricity_results_bymode[4][6] = eccentricity_results_bymode[1][-6]
    eccentricity_results_bymode[4][7] = eccentricity_results_bymode[1][-7]
    eccentricity_results_bymode[4][8] = eccentricity_results_bymode[1][-8]
    eccentricity_results_bymode[4][9] = eccentricity_results_bymode[1][-9]
    eccentricity_results_bymode[4][10] = eccentricity_results_bymode[1][-10]
    eccentricity_results_bymode[5][-9] = eccentricity_results_bymode[0][9]
    eccentricity_results_bymode[5][-8] = eccentricity_results_bymode[0][8]
    eccentricity_results_bymode[5][-7] = eccentricity_results_bymode[0][7]
    eccentricity_results_bymode[5][-6] = eccentricity_results_bymode[0][6]
    eccentricity_results_bymode[5][-5] = eccentricity_results_bymode[0][5]
    eccentricity_results_bymode[5][-4] = eccentricity_results_bymode[0][4]
    eccentricity_results_bymode[5][-3] = eccentricity_results_bymode[0][3]
    eccentricity_results_bymode[5][-2] = eccentricity_results_bymode[0][2]
    eccentricity_results_bymode[5][-1] = eccentricity_results_bymode[0][1]
    eccentricity_results_bymode[5][0] = eccentricity_results_bymode[0][0]
    eccentricity_results_bymode[5][1] = eccentricity_results_bymode[0][-1]
    eccentricity_results_bymode[5][2] = eccentricity_results_bymode[0][-2]
    eccentricity_results_bymode[5][3] = eccentricity_results_bymode[0][-3]
    eccentricity_results_bymode[5][4] = eccentricity_results_bymode[0][-4]
    eccentricity_results_bymode[5][6] = eccentricity_results_bymode[0][-6]
    eccentricity_results_bymode[5][7] = eccentricity_results_bymode[0][-7]
    eccentricity_results_bymode[5][8] = eccentricity_results_bymode[0][-8]
    eccentricity_results_bymode[5][9] = eccentricity_results_bymode[0][-9]
    eccentricity_results_bymode[5][10] = eccentricity_results_bymode[0][-10]

    return eccentricity_results_bymode