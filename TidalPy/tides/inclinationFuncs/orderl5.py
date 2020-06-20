""" Inclination functions (squared) for tidal order-l = 5. These are exact (no truncation on I)
"""

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray



@njit(cacheable=True)
def calc_inclination_off(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp (assuming I=0) for l = 5"""

    # Inclination Functions Calculated for l = 5, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (1, 2) : 3.515625 * ones_,
        (3, 1) : 2756.25 * ones_,
        (5, 0) : 893025. * ones_,
    }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp for l = 5"""

    # Inclination Functions Calculated for l = 5.
    # Optimizations
    i = inclination
    i_half = i / 2.
    i_double = 2. * i
    i_triple = 3. * i
    sin_i = np.sin(i)
    cos_i = np.cos(i)
    sin_i_half = np.sin(i_half)
    cos_i_half = np.cos(i_half)
    cos_i_double = np.cos(i_double)
    cos_i_triple = np.cos(i_triple)

    inclination_results = {
        (0, 0) : 62.015625*sin_i_half**10*cos_i_half**10,
        (0, 1) : 0.0186920166015625*(9.0*sin_i**2 - 8.0)**2*sin_i**6,
        (0, 2) : (-0.05859375*(cos_i + 1.0)**4*sin_i - 1.875*sin_i_half**9*cos_i_half + 18.75*sin_i_half**7*cos_i_half**3 - 37.5*sin_i_half**5*cos_i_half**5 + 18.75*sin_i_half**3*cos_i_half**7)**2,
        (0, 3) : (0.05859375*(cos_i + 1.0)**4*sin_i + 1.875*sin_i_half**9*cos_i_half - 18.75*sin_i_half**7*cos_i_half**3 + 37.5*sin_i_half**5*cos_i_half**5 - 18.75*sin_i_half**3*cos_i_half**7)**2,
        (0, 4) : 0.0186920166015625*(8.0 - 9.0*sin_i**2)**2*sin_i**6,
        (0, 5) : 62.015625*sin_i_half**10*cos_i_half**10,
        (1, 0) : 1550.390625*sin_i_half**8*cos_i_half**12,
        (1, 1) : 0.01051425933837890625*(cos_i + 1.0)**4*(-65.0*cos_i + 42.0*cos_i_double - 15.0*cos_i_triple + 38.0)**2,
        (1, 2) : 3.515625*(-185.0*sin_i_half**10 + 445.0*sin_i_half**8 - 350.0*sin_i_half**6 + 90.0*sin_i_half**4 + 25.0*cos_i_half**10 - 24.0*cos_i_half**8)**2,
        (1, 3) : 3.515625*(115.0*sin_i_half**8 - 204.0*sin_i_half**6 + 90.0*sin_i_half**4 + 95.0*cos_i_half**8 - 80.0*cos_i_half**6)**2*sin_i_half**4,
        (1, 4) : 2.691650390625*(cos_i + 1.0)**2*(-15.0*sin_i**2 + 6.0*cos_i + 14.0)**2*sin_i_half**8,
        (1, 5) : 1550.390625*sin_i_half**12*cos_i_half**8,
        (2, 0) : 24806.25*sin_i_half**6*cos_i_half**14,
        (2, 1) : (-1.640625*(cos_i + 1.0)**4*sin_i - 367.5*sin_i_half**5*cos_i_half**5 + 367.5*sin_i_half**3*cos_i_half**7)**2,
        (2, 2) : (1.640625*(cos_i + 1.0)**4*sin_i - 262.5*sin_i_half**7*cos_i_half**3 + 787.5*sin_i_half**5*cos_i_half**5 - 472.5*sin_i_half**3*cos_i_half**7)**2,
        (2, 3) : 2756.25*(-25.0*sin_i_half**6 + 39.0*sin_i_half**4 - 15.0*sin_i_half**2 + 5.0*cos_i_half**6)**2*sin_i_half**6*cos_i_half**2,
        (2, 4) : 2756.25*(15.0*sin_i_half**4 - 21.0*sin_i_half**2 + 7.0)**2*sin_i_half**10*cos_i_half**2,
        (2, 5) : 24806.25*sin_i_half**14*cos_i_half**6,
        (3, 0) : 223256.25*sin_i_half**4*cos_i_half**16,
        (3, 1) : 2.691650390625*(cos_i + 1.0)**6*(45.0*cos_i**2 - 54.0*cos_i + 13.0)**2,
        (3, 2) : 6.05621337890625*(cos_i + 1.0)**4*(-65.0*cos_i + 42.0*cos_i_double - 15.0*cos_i_triple + 38.0)**2,
        (3, 3) : 1550.390625*(cos_i + 1.0)**2*(-15.0*sin_i**2 + 6.0*cos_i + 14.0)**2*sin_i_half**8,
        (3, 4) : 172.265625*(-45.0*sin_i**2 + 54.0*cos_i + 58.0)**2*sin_i_half**12,
        (3, 5) : 223256.25*sin_i_half**16*cos_i_half**4,
        (4, 0) : 872.0947265625*(cos_i + 1.0)**8*sin_i**2,
        (4, 1) : (29.53125*(cos_i + 1.0)**4*sin_i - 3780.0*sin_i_half**3*cos_i_half**7)**2,
        (4, 2) : (4725.0*cos_i - 945.0)**2*sin_i_half**6*cos_i_half**10,
        (4, 3) : (4725.0*cos_i + 945.0)**2*sin_i_half**10*cos_i_half**6,
        (4, 4) : 223256.25*(5.0*cos_i + 3.0)**2*sin_i_half**14*cos_i_half**2,
        (4, 5) : 893025.0*sin_i_half**18*cos_i_half**2,
        (5, 0) : 893025.0*cos_i_half**20,
        (5, 1) : 22325625.0*sin_i_half**4*cos_i_half**16,
        (5, 2) : 89302500.0*sin_i_half**8*cos_i_half**12,
        (5, 3) : 89302500.0*sin_i_half**12*cos_i_half**8,
        (5, 4) : 22325625.0*sin_i_half**16*cos_i_half**4,
        (5, 5) : 893025.0*sin_i_half**20
    }

    return inclination_results
