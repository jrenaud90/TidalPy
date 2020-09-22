""" Inclination functions (squared) for tidal order-l = 6. These are exact (no truncation on I)
"""

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray


@njit(cacheable=True)
def calc_inclination_off(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp (assuming I=0) for l = 6"""

    # Inclination Functions Calculated for l = 6, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (0, 3) : 0.09765625 * ones_,
        (2, 2) : 172.265625 * ones_,
        (4, 1) : 223256.25 * ones_,
        (6, 0) : 108056025. * ones_,
    }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp for l = 6"""

    # Inclination Functions Calculated for l = 6.
    # Optimizations
    i = inclination
    i_half = i / 2.
    i_double = 2. * i
    i_triple = 3. * i
    sin_i = np.sin(i)
    cos_i = np.cos(i)
    sin_i_half = np.sin(i_half)
    cos_i_half = np.cos(i_half)
    sin_i_double = np.sin(i_double)
    cos_i_double = np.cos(i_double)
    cos_i_triple = np.cos(i_triple)

    inclination_results = {
        (0, 0) : 208.44140625*sin_i_half**12*cos_i_half**12,
        (0, 1) : 0.015140533447265625*(10.0 - 11.0*sin_i**2)**2*sin_i**8,
        (0, 2) : 43.06640625*(-24.0*sin_i_half**10 + 62.0*sin_i_half**8 - 53.0*sin_i_half**6 + 15.0*sin_i_half**4 + 9.0*cos_i_half**10 - 8.0*cos_i_half**8)**2*sin_i_half**4,
        (0, 3) : 0.09765625*(262.0*sin_i_half**12 - 486.0*sin_i_half**10 + 225.0*sin_i_half**8 - 400.0*sin_i_half**6*cos_i_half**6 + 225.0*sin_i_half**4*cos_i_half**8 + 37.0*cos_i_half**12 - 36.0*cos_i_half**10)**2,
        (0, 4) : 43.06640625*(-24.0*sin_i_half**10 + 62.0*sin_i_half**8 - 53.0*sin_i_half**6 + 15.0*sin_i_half**4 + 9.0*cos_i_half**10 - 8.0*cos_i_half**8)**2*sin_i_half**4,
        (0, 5) : 0.015140533447265625*(10.0 - 11.0*sin_i**2)**2*sin_i**8,
        (0, 6) : 208.44140625*sin_i_half**12*cos_i_half**12,
        (1, 0) : 7503.890625*sin_i_half**10*cos_i_half**14,
        (1, 1) : 62.015625*(-66.0*sin_i_half**4 + 55.0*sin_i_half**2 - 10.0)**2*sin_i_half**6*cos_i_half**10,
        (1, 2) : (-0.205078125*(cos_i + 1.0)**5*sin_i - 91.875*sin_i_half**9*cos_i_half**3 + 459.375*sin_i_half**7*cos_i_half**5 - 551.25*sin_i_half**5*cos_i_half**7 + 183.75*sin_i_half**3*cos_i_half**9)**2,
        (1, 3) : (0.205078125*(cos_i + 1.0)**5*sin_i - 13.125*sin_i_half**11*cos_i_half + 196.875*sin_i_half**9*cos_i_half**3 - 656.25*sin_i_half**7*cos_i_half**5 + 656.25*sin_i_half**5*cos_i_half**7 - 196.875*sin_i_half**3*cos_i_half**9)**2,
        (1, 4) : 172.265625*(57.0*sin_i_half**8 - 98.0*sin_i_half**6 + 42.0*sin_i_half**4 + 42.0*cos_i_half**8 - 35.0*cos_i_half**6)**2*sin_i_half**6*cos_i_half**2,
        (1, 5) : 62.015625*(66.0*sin_i_half**4 - 77.0*sin_i_half**2 + 21.0)**2*sin_i_half**10*cos_i_half**6,
        (1, 6) : 7503.890625*sin_i_half**14*cos_i_half**10,
        (2, 0) : 187597.265625*sin_i_half**8*cos_i_half**16,
        (2, 1) : 0.09462833404541015625*(cos_i + 1.0)**6*(-191.0*cos_i + 110.0*cos_i_double - 33.0*cos_i_triple + 114.0)**2,
        (2, 2) : 172.265625*(42.0*(cos_i - 1.0)**2 + 462.0*sin_i_half**8 - 560.0*sin_i_half**6 + 33.0*cos_i_half**8 - 32.0*cos_i_half**6)**2*cos_i_half**8,
        (2, 3) : 68906.25*(-24.0*sin_i_half**10 + 62.0*sin_i_half**8 - 53.0*sin_i_half**6 + 15.0*sin_i_half**4 + 9.0*cos_i_half**10 - 8.0*cos_i_half**8)**2*sin_i_half**4,
        (2, 4) : 172.265625*(201.0*sin_i_half**8 - 368.0*sin_i_half**6 + 168.0*sin_i_half**4 + 294.0*cos_i_half**8 - 224.0*cos_i_half**6)**2*sin_i_half**8,
        (2, 5) : 96.8994140625*(cos_i + 1.0)**2*(-33.0*sin_i**2 + 22.0*cos_i + 34.0)**2*sin_i_half**12,
        (2, 6) : 187597.265625*sin_i_half**16*cos_i_half**8,
        (3, 0) : 3001556.25*sin_i_half**6*cos_i_half**18,
        (3, 1) : (-7.3828125*(cos_i + 1.0)**5*sin_i - 5670.0*sin_i_half**5*cos_i_half**7 + 4252.5*sin_i_half**3*cos_i_half**9)**2,
        (3, 2) : (7.3828125*(cos_i + 1.0)**5*sin_i - 6615.0*sin_i_half**7*cos_i_half**5 + 13230.0*sin_i_half**5*cos_i_half**7 - 5670.0*sin_i_half**3*cos_i_half**9)**2,
        (3, 3) : 155039.0625*(11.0*sin_i_half**4 - 11.0*sin_i_half**2 + 2.0)**2*sin_i_half**4*sin_i_double**2*cos_i_half,
        (3, 4) : 223256.25*(-41.0*sin_i_half**6 + 68.0*sin_i_half**4 - 28.0*sin_i_half**2 + 14.0*cos_i_half**6)**2*sin_i_half**10*cos_i_half**2,
        (3, 5) : 223256.25*(22.0*sin_i_half**4 - 33.0*sin_i_half**2 + 12.0)**2*sin_i_half**14*cos_i_half**2,
        (3, 6) : 3001556.25*sin_i_half**18*cos_i_half**6,
        (4, 0) : 27014006.25*sin_i_half**4*cos_i_half**20,
        (4, 1) : 218.023681640625*(cos_i + 1.0)**8*(-33.0*sin_i**2 - 44.0*cos_i + 46.0)**2,
        (4, 2) : 85.165500640869140625*(cos_i + 1.0)**6*(-191.0*cos_i + 110.0*cos_i_double - 33.0*cos_i_triple + 114.0)**2,
        (4, 3) : (23625.0 - 25987.5*sin_i**2)**2*sin_i_half**8*cos_i_half**8,
        (4, 4) : 87209.47265625*(cos_i + 1.0)**2*(-33.0*sin_i**2 + 22.0*cos_i + 34.0)**2*sin_i_half**12,
        (4, 5) : (-7796.25*sin_i**2 + 10395.0*cos_i + 10867.5)**2*sin_i_half**16,
        (4, 6) : 27014006.25*sin_i_half**20*cos_i_half**4,
        (5, 0) : 26380.865478515625*(cos_i + 1.0)**10*sin_i**2,
        (5, 1) : (162.421875*(cos_i + 1.0)**5*sin_i - 51975.0*sin_i_half**3*cos_i_half**9)**2,
        (5, 2) : 675350156.25*(3.0*cos_i - 1.0)**2*sin_i_half**6*cos_i_half**14,
        (5, 3) : 10805602500.0*sin_i_half**10*cos_i_half**10*cos_i**2,
        (5, 4) : 675350156.25*(3.0*cos_i + 1.0)**2*sin_i_half**14*cos_i_half**6,
        (5, 5) : 105523.4619140625*(cos_i - 1.0)**8*(3.0*cos_i + 2.0)**2*sin_i**2,
        (5, 6) : 108056025.0*sin_i_half**22*cos_i_half**2,
        (6, 0) : 108056025.0*cos_i_half**24,
        (6, 1) : 3890016900.0*sin_i_half**4*cos_i_half**20,
        (6, 2) : 24312605625.0*sin_i_half**8*cos_i_half**16,
        (6, 3) : 43222410000.0*sin_i_half**12*cos_i_half**12,
        (6, 4) : 24312605625.0*sin_i_half**16*cos_i_half**8,
        (6, 5) : 3890016900.0*sin_i_half**20*cos_i_half**4,
        (6, 6) : 108056025.0*sin_i_half**24
    }

    return inclination_results
