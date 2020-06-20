""" Inclination functions (squared) for tidal order-l = 4. These are exact (no truncation on I)
"""

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray



@njit(cacheable=True)
def calc_inclination_off(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp (assuming I=0) for l = 4"""

    # Inclination Functions Calculated for l = 4, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (0, 2) : 0.140625 * ones_,
        (2, 1) : 56.25 * ones_,
        (4, 0) : 11025. * ones_,
    }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp for l = 4"""

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

    inclination_results = {
        (0, 0) : 19.140625*sin_i_half**8*cos_i_half**8,
        (0, 1) : 0.09765625*(7.0*sin_i**2 - 6.0)**2*(cos_i + 1.0)**2*sin_i_half**4,
        (0, 2) : 0.140625*(53.0*sin_i_half**8 - 88.0*sin_i_half**6 + 36.0*sin_i_half**4 + 17.0*cos_i_half**8 - 16.0*cos_i_half**6)**2,
        (0, 3) : 0.09765625*(7.0*sin_i**2 - 6.0)**2*(cos_i + 1.0)**2*sin_i_half**4,
        (0, 4) : 19.140625*sin_i_half**8*cos_i_half**8,
        (1, 0) : 306.25*sin_i_half**6*cos_i_half**10,
        (1, 1) : (-0.46875*(cos_i + 1.0)**3*sin_i - 25.0*sin_i_half**5*cos_i_half**3 + 37.5*sin_i_half**3*cos_i_half**5)**2,
        (1, 2) : (0.46875*(cos_i + 1.0)**3*sin_i - 7.5*sin_i_half**7*cos_i_half + 45.0*sin_i_half**5*cos_i_half**3 - 45.0*sin_i_half**3*cos_i_half**5)**2,
        (1, 3) : 1.5625*(-14.0*sin_i**2 + 7.0*cos_i + 13.0)**2*sin_i_half**6*cos_i_half**2,
        (1, 4) : 306.25*sin_i_half**10*cos_i_half**6,
        (2, 0) : 2756.25*sin_i_half**4*cos_i_half**12,
        (2, 1) : 3.515625*(cos_i + 1.0)**4*(-7.0*sin_i**2 - 7.0*cos_i + 8.0)**2,
        (2, 2) : 31.640625*(7.0*sin_i**2 - 6.0)**2*(cos_i + 1.0)**2*sin_i_half**4,
        (2, 3) : (-52.5*sin_i**2 + 52.5*cos_i + 60.0)**2*sin_i_half**8,
        (2, 4) : 2756.25*sin_i_half**12*cos_i_half**4,
        (3, 0) : 43.06640625*(cos_i + 1.0)**6*sin_i**2,
        (3, 1) : (6.5625*(cos_i + 1.0)**3*sin_i - 315.0*sin_i_half**3*cos_i_half**5)**2,
        (3, 2) : 99225.0*sin_i_half**6*cos_i_half**6*cos_i**2,
        (3, 3) : 172.265625*(sin_i + sin_i_double)**2*(cos_i - 1.0)**4,
        (3, 4) : 11025.0*sin_i_half**14*cos_i_half**2,
        (4, 0) : 11025.0*cos_i_half**16,
        (4, 1) : 176400.0*sin_i_half**4*cos_i_half**12,
        (4, 2) : 396900.0*sin_i_half**8*cos_i_half**8,
        (4, 3) : 176400.0*sin_i_half**12*cos_i_half**4,
        (4, 4) : 11025.0*sin_i_half**16
    }

    return inclination_results
