""" Inclination functions (squared) for tidal order-l = 4. These are exact (no truncation on I)
"""

import numpy as np

from TidalPy.performance import njit


@njit
def calc_inclination_off(inclination):

    # Inclination Functions Calculated for l = 4, Inclination == off.
    ones_ = np.ones_like(inclination)

    # Output is structured (p0_(m0_result, m1_result, ...), p1_(...), ...)
    inclination_results = {
        0 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_,
            3 : 0 * ones_,
            4 : 11025. * ones_
        },
        1 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 56.25 * ones_,
            3 : 0 * ones_,
            4 : 0 * ones_
        },
        2 : {
            0 : 0.140625 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_,
            3 : 0 * ones_,
            4 : 0 * ones_
        },
        3 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_,
            3 : 0 * ones_,
            4 : 0 * ones_
        },
        4 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_,
            3 : 0 * ones_,
            4 : 0 * ones_
        }
    }

    return inclination_results


@njit
def calc_inclination(inclination):

    # Inclination Functions Calculated for l = 4.
    # Optimizations
    i = inclination
    i_half = i / 2.
    i_double = 2. * i
    sin_i = np.sin(i)
    cos_i = np.cos(i)
    sin_i_half = np.sin(i_half)
    cos_i_half = np.cos(i_half)
    sin_i_double = np.sin(i_double)

    # Output is structured (p0_(m0_result, m1_result, ...), p1_(...), ...)
    inclination_results = {
        0 : {
            0 : 19.140625*sin_i_half**8*cos_i_half**8,
            1 : 306.25*sin_i_half**6*cos_i_half**10,
            2 : 2756.25*sin_i_half**4*cos_i_half**12,
            3 : 43.06640625*(cos_i + 1.0)**6*sin_i**2,
            4 : 11025.0*cos_i_half**16
        },
        1 : {
            0 : 4.78515625*(sin_i**2 - 0.8571428571428571428571429)**2*(cos_i + 1.0)**2*sin_i_half**4,
            1 : 1406.25*(-0.0125*(cos_i + 1.0)**3*sin_i - 0.6666666666666666666666667*sin_i_half**5*cos_i_half**3 + sin_i_half**3*cos_i_half**5)**2,
            2 : 225.0*(cos_i + 1.0)**4*(-0.875*sin_i**2 - 0.875*cos_i + 1)**2,
            3 : 99225.0*(0.02083333333333333333333333*(cos_i + 1.0)**3*sin_i - sin_i_half**3*cos_i_half**5)**2,
            4 : 176400.0*sin_i_half**4*cos_i_half**12
        },
        2 : {
            0 : 1089.0*(0.6022727272727272727272727*sin_i_half**8 - sin_i_half**6 + 0.4090909090909090909090909*sin_i_half**4 + 0.1931818181818181818181818*cos_i_half**8 - 0.1818181818181818181818182*cos_i_half**6)**2,
            1 : 2025.0*(0.01041666666666666666666667*(cos_i + 1.0)**3*sin_i - 0.1666666666666666666666667*sin_i_half**7*cos_i_half + sin_i_half**5*cos_i_half**3 - sin_i_half**3*cos_i_half**5)**2,
            2 : 1550.390625*(sin_i**2 - 0.8571428571428571428571429)**2*(cos_i + 1.0)**2*sin_i_half**4,
            3 : 99225.0*sin_i_half**6*cos_i_half**6*cos_i**2,
            4 : 396900.0*sin_i_half**8*cos_i_half**8
        },
        3 : {
            0 : 4.78515625*(sin_i**2 - 0.8571428571428571428571429)**2*(cos_i + 1.0)**2*sin_i_half**4,
            1 : 306.25*(-sin_i**2 + 0.5*cos_i + 0.9285714285714285714285714)**2*sin_i_half**6*cos_i_half**2,
            2 : 3600.0*(-0.875*sin_i**2 + 0.875*cos_i + 1)**2*sin_i_half**8,
            3 : 172.265625*(sin_i + sin_i_double)**2*(cos_i - 1.0)**4,
            4 : 176400.0*sin_i_half**12*cos_i_half**4
        },
        4 : {
            0 : 19.140625*sin_i_half**8*cos_i_half**8,
            1 : 306.25*sin_i_half**10*cos_i_half**6,
            2 : 2756.25*sin_i_half**12*cos_i_half**4,
            3 : 11025.0*sin_i_half**14*cos_i_half**2,
            4 : 11025.0*sin_i_half**16
        }
    }

    return inclination_results