""" Inclination functions (squared) for tidal order-l = 4. These are exact (no truncation on I)
"""

from typing import TYPE_CHECKING

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def calc_inclination_off(inclination: 'FloatArray') -> 'InclinOutput':
    """Calculate F^2_lmp (assuming I=0) for l = 4"""

    # Inclination Functions Calculated for l = 4, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (0, 2): 0.140625 * ones_,
        (2, 1): 56.25 * ones_,
        (4, 0): 11025. * ones_,
        }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: 'FloatArray') -> 'InclinOutput':
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
        (0, 1) : 4.78515625*(sin_i**2 - 0.8571428571428571428571429)**2*(cos_i + 1.0)**2*sin_i_half**4,
        (0, 2) : 1089.0*(0.6022727272727272727272727*sin_i_half**8 - sin_i_half**6 + 0.4090909090909090909090909*sin_i_half**4 + 0.1931818181818181818181818*cos_i_half**8 - 0.1818181818181818181818182*cos_i_half**6)**2,
        (0, 3) : 4.78515625*(sin_i**2 - 0.8571428571428571428571429)**2*(cos_i + 1.0)**2*sin_i_half**4,
        (0, 4) : 19.140625*sin_i_half**8*cos_i_half**8,
        (1, 0) : 306.25*sin_i_half**6*cos_i_half**10,
        (1, 1) : 1406.25*(-0.0125*(cos_i + 1.0)**3*sin_i - 0.6666666666666666666666667*sin_i_half**5*cos_i_half**3 + sin_i_half**3*cos_i_half**5)**2,
        (1, 2) : 2025.0*(0.01041666666666666666666667*(cos_i + 1.0)**3*sin_i - 0.1666666666666666666666667*sin_i_half**7*cos_i_half + sin_i_half**5*cos_i_half**3 - sin_i_half**3*cos_i_half**5)**2,
        (1, 3) : 306.25*(-sin_i**2 + 0.5*cos_i + 0.9285714285714285714285714)**2*sin_i_half**6*cos_i_half**2,
        (1, 4) : 306.25*sin_i_half**10*cos_i_half**6,
        (2, 0) : 2756.25*sin_i_half**4*cos_i_half**12,
        (2, 1) : 225.0*(cos_i + 1.0)**4*(0.875*sin_i**2 + 0.875*cos_i - 1)**2,
        (2, 2) : 1550.390625*(sin_i**2 - 0.8571428571428571428571429)**2*(cos_i + 1.0)**2*sin_i_half**4,
        (2, 3) : 3600.0*(0.875*sin_i**2 - 0.875*cos_i - 1)**2*sin_i_half**8,
        (2, 4) : 2756.25*sin_i_half**12*cos_i_half**4,
        (3, 0) : 43.06640625*(cos_i + 1.0)**6*sin_i**2,
        (3, 1) : 99225.0*(-0.02083333333333333333333333*(cos_i + 1.0)**3*sin_i + sin_i_half**3*cos_i_half**5)**2,
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
