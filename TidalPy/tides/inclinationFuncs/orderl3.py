""" Inclination functions (squared) for tidal order-l = 3. These are exact (no truncation on I)
"""

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray



@njit(cacheable=True)
def calc_inclination_off(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp (assuming I=0) for l = 3"""

    # Inclination Functions Calculated for l = 3, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (1, 1) : 2.25 * ones_,
        (3, 0) : 225. * ones_,
    }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: FloatArray) -> InclinOutput:
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
        (0, 1) : 0.03515625*(5.0*sin_i**2 - 4.0)**2*sin_i**2,
        (0, 2) : 0.03515625*(4.0 - 5.0*sin_i**2)**2*sin_i**2,
        (0, 3) : 0.09765625*sin_i**6,
        (1, 0) : 56.25*sin_i_half**4*cos_i_half**8,
        (1, 1) : 0.140625*(15.0*cos_i**2 - 10.0*cos_i - 1.0)**2*cos_i_half**4,
        (1, 2) : 0.140625*(-15.0*sin_i**2 + 10.0*cos_i + 14.0)**2*sin_i_half**4,
        (1, 3) : 56.25*sin_i_half**8*cos_i_half**4,
        (2, 0) : 3.515625*(cos_i + 1.0)**4*sin_i**2,
        (2, 1) : 3.515625*(-3.0*sin_i**2 + 2.0*cos_i + 2.0)**2*sin_i**2,
        (2, 2) : 225.0*(-sin_i_half**5*cos_i_half + 0.25*sin_i**3)**2,
        (2, 3) : 225.0*sin_i_half**10*cos_i_half**2,
        (3, 0) : 225.0*cos_i_half**12,
        (3, 1) : 2025.0*sin_i_half**4*cos_i_half**8,
        (3, 2) : 2025.0*sin_i_half**8*cos_i_half**4,
        (3, 3) : 225.0*sin_i_half**12
    }

    return inclination_results
