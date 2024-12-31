""" Inclination functions (squared) for tidal order-l = 3. These are exact (no truncation on I)
"""

from typing import TYPE_CHECKING

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def calc_inclination_off(inclination: 'FloatArray') -> 'InclinOutput':
    """Calculate F^2_lmp (assuming I=0) for l = 3"""

    # Inclination Functions Calculated for l = 3, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (1, 1): 2.25 * ones_,
        (3, 0): 225. * ones_,
        }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: 'FloatArray') -> 'InclinOutput':
    """Calculate F^2_lmp for l = 3"""

    # Inclination Functions Calculated for l = 3.
    # Optimizations
    i = inclination
    i_half = i / 2.
    sin_i = np.sin(i)
    cos_i = np.cos(i)
    sin_i_half = np.sin(i_half)
    cos_i_half = np.cos(i_half)

    inclination_results = {
        (0, 0) : 0.09765625*sin_i**6,
        (0, 1) : 0.87890625*(sin_i**2 - 0.8)**2*sin_i**2,
        (0, 2) : 0.87890625*(0.8 - sin_i**2)**2*sin_i**2,
        (0, 3) : 0.09765625*sin_i**6,
        (1, 0) : 56.25*sin_i_half**4*cos_i_half**8,
        (1, 1) : 31.640625*(-cos_i**2 + 0.6666666666666666666666667*cos_i + 0.06666666666666666666666667)**2*cos_i_half**4,
        (1, 2) : 31.640625*(sin_i**2 - 0.6666666666666666666666667*cos_i - 0.9333333333333333333333333)**2*sin_i_half**4,
        (1, 3) : 56.25*sin_i_half**8*cos_i_half**4,
        (2, 0) : 3.515625*(cos_i + 1.0)**4*sin_i**2,
        (2, 1) : 31.640625*(sin_i**2 - 0.6666666666666666666666667*cos_i - 0.6666666666666666666666667)**2*sin_i**2,
        (2, 2) : 225.0*(sin_i_half**5*cos_i_half - 0.25*sin_i**3)**2,
        (2, 3) : 225.0*sin_i_half**10*cos_i_half**2,
        (3, 0) : 225.0*cos_i_half**12,
        (3, 1) : 2025.0*sin_i_half**4*cos_i_half**8,
        (3, 2) : 2025.0*sin_i_half**8*cos_i_half**4,
        (3, 3) : 225.0*sin_i_half**12
    }

    return inclination_results