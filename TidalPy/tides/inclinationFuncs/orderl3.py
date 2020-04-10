""" Inclination functions (squared) for tidal order-l = 3. These are exact (no truncation on I)
"""

import numpy as np

from TidalPy.performance import njit

@njit
def calc_inclination_off(inclination):

    # Inclination Functions Calculated for l = 3, Inclination == off.
    ones_ = np.ones_like(inclination)

    # Output is structured (p0_(m0_result, m1_result, ...), p1_(...), ...)
    inclination_results = {
        0 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_,
            3 : 225. * ones_
        },
        1 : {
            0 : 0 * ones_,
            1 : 2.25 * ones_,
            2 : 0 * ones_,
            3 : 0 * ones_
        },
        2 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_,
            3 : 0 * ones_
        },
        3 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_,
            3 : 0 * ones_
        }
    }

    return inclination_results


@njit
def calc_inclination(inclination):

    # Inclination Functions Calculated for l = 3.
    # Optimizations
    i = inclination
    i_half = i / 2.
    sin_i = np.sin(i)
    cos_i = np.cos(i)
    sin_i_half = np.sin(i_half)
    cos_i_half = np.cos(i_half)

    # Output is structured (p0_(m0_result, m1_result, ...), p1_(...), ...)
    inclination_results = {
        0 : {
            0 : 0.09765625*sin_i**6,
            1 : 56.25*sin_i_half**4*cos_i_half**8,
            2 : 3.515625*(cos_i + 1.0)**4*sin_i**2,
            3 : 225.0*cos_i_half**12
        },
        1 : {
            0 : 0.87890625*(sin_i**2 - 0.8)**2*sin_i**2,
            1 : 31.640625*(cos_i**2 - 0.6666666666666666666666667*cos_i - 0.06666666666666666666666667)**2*cos_i_half**4,
            2 : 31.640625*(-sin_i**2 + 0.6666666666666666666666667*cos_i + 0.6666666666666666666666667)**2*sin_i**2,
            3 : 2025.0*sin_i_half**4*cos_i_half**8
        },
        2 : {
            0 : 0.87890625*(0.8 - sin_i**2)**2*sin_i**2,
            1 : 31.640625*(-sin_i**2 + 0.6666666666666666666666667*cos_i + 0.9333333333333333333333333)**2*sin_i_half**4,
            2 : 225.0*(-sin_i_half**5*cos_i_half + 0.25*sin_i**3)**2,
            3 : 2025.0*sin_i_half**8*cos_i_half**4
        },
        3 : {
            0 : 0.09765625*sin_i**6,
            1 : 56.25*sin_i_half**8*cos_i_half**4,
            2 : 225.0*sin_i_half**10*cos_i_half**2,
            3 : 225.0*sin_i_half**12
        }
    }

    return inclination_results