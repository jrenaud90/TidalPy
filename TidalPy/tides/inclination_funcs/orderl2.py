""" Inclination functions (squared) for tidal order-l = 2. These are exact (no truncation on I)
"""

from typing import TYPE_CHECKING

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def calc_inclination_off(inclination: 'FloatArray') -> 'InclinOutput':
    """Calculate F^2_lmp (assuming I=0) for l = 2"""

    # Inclination Functions Calculated for l = 2, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (0, 1): 0.25 * ones_,
        (2, 0): 9. * ones_,
        }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: 'FloatArray') -> 'InclinOutput':
    """Calculate F^2_lmp for l = 2"""

    # Inclination Functions Calculated for l = 2.
    # Optimizations
    i = inclination
    i_half = i / 2.
    i_double = 2. * i
    sin_i = np.sin(i)
    sin_i_half = np.sin(i_half)
    cos_i_half = np.cos(i_half)
    sin_i_double = np.sin(i_double)

    inclination_results = {
        (0, 0) : 0.140625*sin_i**4,
        (0, 1) : (-sin_i_half**4 + sin_i_half**2 + 0.5*sin_i**2 - 0.5)**2,
        (0, 2) : 0.140625*sin_i**4,
        (1, 0) : 9.0*sin_i_half**2*cos_i_half**6,
        (1, 1) : 0.5625*sin_i_double**2,
        (1, 2) : 9.0*sin_i_half**6*cos_i_half**2,
        (2, 0) : 9.0*cos_i_half**8,
        (2, 1) : 2.25*sin_i**4,
        (2, 2) : 9.0*sin_i_half**8
    }

    return inclination_results
