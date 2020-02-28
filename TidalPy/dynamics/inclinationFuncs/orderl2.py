import numpy as np

from ...performance import njit

@njit
def calc_inclination_off(inclination):

    # Inclination Functions Calculated for l = 2, Inclination == off.
    ones_ = np.ones_like(inclination)

    # Output is structured (p0_(m0_result, m1_result, ...), p1_(...), ...)
    inclination_results = {
        0 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 9.000000000000000000000000 * ones_
        },
        1 : {
            0 : 0.2500000000000000000000000 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_
        },
        2 : {
            0 : 0 * ones_,
            1 : 0 * ones_,
            2 : 0 * ones_
        }
    }

    return inclination_results


@njit
def calc_inclination(inclination):

    # Inclination Functions Calculated for l = 2.
    # Optimizations
    i = inclination
    i_half = i / 2.
    i_double = 2. * i
    sin_i = np.sin(i)
    sin_i_half = np.sin(i_half)
    cos_i_half = np.cos(i_half)
    sin_i_double = np.sin(i_double)

    # Output is structured (p0_(m0_result, m1_result, ...), p1_(...), ...)
    inclination_results = {
        0 : {
            0 : 0.140625*sin_i**4,
            1 : 9.0*sin_i_half**2*cos_i_half**6,
            2 : 9.0*cos_i_half**8
        },
        1 : {
            0 : (sin_i_half**4 - sin_i_half**2 - 0.5*sin_i**2 + 0.5)**2,
            1 : 0.5625*sin_i_double**2,
            2 : 2.25*sin_i**4
        },
        2 : {
            0 : 0.140625*sin_i**4,
            1 : 9.0*sin_i_half**6*cos_i_half**2,
            2 : 9.0*sin_i_half**8
        }
    }

    return inclination_results

