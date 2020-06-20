""" Inclination functions (squared) for tidal order-l = 7. These are exact (no truncation on I)
"""

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray



@njit(cacheable=True)
def calc_inclination_off(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp (assuming I=0) for l = 7"""

    # Inclination Functions Calculated for l = 7, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (1, 3) : 4.78515625 * ones_,
        (3, 2) : 13953.515625 * ones_,
        (5, 1) : 27014006.25 * ones_,
        (7, 0) : 18261468225. * ones_,
    }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: FloatArray) -> InclinOutput:
    """Calculate F^2_lmp for l = 7"""

    # Inclination Functions Calculated for l = 7.
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

    # TODO: Optimizations - Lots of powers duplicated in the below. Could be precalculated and turned into references.

    inclination_results = {
        (0, 0) : 718.91015625*sin_i_half**14*cos_i_half**14,
        (0, 1) : 208.44140625*(-13.0*sin_i_half**4 + 13.0*sin_i_half**2 - 3.0)**2*sin_i_half**10*cos_i_half**10,
        (0, 2) : 15.50390625*(-103.0*sin_i_half**8 + 161.0*sin_i_half**6 - 63.0*sin_i_half**4 - 40.0*cos_i_half**8 + 35.0*cos_i_half**6)**2*sin_i_half**6*cos_i_half**6,
        (0, 3) : (-0.01708984375*(cos_i + 1.0)**6*sin_i - 2.1875*sin_i_half**13*cos_i_half + 45.9375*sin_i_half**11*cos_i_half**3 - 229.6875*sin_i_half**9*cos_i_half**5 + 382.8125*sin_i_half**7*cos_i_half**7 - 229.6875*sin_i_half**5*cos_i_half**9 + 45.9375*sin_i_half**3*cos_i_half**11)**2,
        (0, 4) : (0.01708984375*(cos_i + 1.0)**6*sin_i + 2.1875*sin_i_half**13*cos_i_half - 45.9375*sin_i_half**11*cos_i_half**3 + 229.6875*sin_i_half**9*cos_i_half**5 - 382.8125*sin_i_half**7*cos_i_half**7 + 229.6875*sin_i_half**5*cos_i_half**9 - 45.9375*sin_i_half**3*cos_i_half**11)**2,
        (0, 5) : 15.50390625*(103.0*sin_i_half**8 - 161.0*sin_i_half**6 + 63.0*sin_i_half**4 + 40.0*cos_i_half**8 - 35.0*cos_i_half**6)**2*sin_i_half**6*cos_i_half**6,
        (0, 6) : 208.44140625*(13.0*sin_i_half**4 - 13.0*sin_i_half**2 + 3.0)**2*sin_i_half**10*cos_i_half**10,
        (0, 7) : 718.91015625*sin_i_half**14*cos_i_half**14,
        (1, 0) : 35226.59765625*sin_i_half**12*cos_i_half**16,
        (1, 1) : 0.0127222537994384765625*(-65.0*sin_i**2 + 91.0*cos_i**3 - 31.0*cos_i + 60.0)**2*sin_i**8,
        (1, 2) : 15.50390625*(826.0*sin_i_half**8 - 1176.0*sin_i_half**6 + 420.0*sin_i_half**4 + 175.0*cos_i_half**8 - 160.0*cos_i_half**6)**2*sin_i_half**4*cos_i_half**8,
        (1, 3) : 4.78515625*(-1414.0*sin_i_half**14 + 3850.0*sin_i_half**12 - 3486.0*sin_i_half**10 + 1050.0*sin_i_half**8 - 1120.0*sin_i_half**6*cos_i_half**8 + 420.0*sin_i_half**4*cos_i_half**10 + 49.0*cos_i_half**14 - 48.0*cos_i_half**12)**2,
        (1, 4) : 4.78515625*(469.0*sin_i_half**12 - 888.0*sin_i_half**10 + 420.0*sin_i_half**8 - 1120.0*sin_i_half**6*cos_i_half**6 + 1050.0*sin_i_half**4*cos_i_half**8 + 364.0*cos_i_half**12 - 336.0*cos_i_half**10)**2*sin_i_half**4,
        (1, 5) : 15.50390625*(-595.0*sin_i_half**10 + 1595.0*sin_i_half**8 - 1420.0*sin_i_half**6 + 420.0*sin_i_half**4 + 406.0*cos_i_half**10 - 336.0*cos_i_half**8)**2*sin_i_half**8,
        (1, 6) : 0.0127222537994384765625*(-65.0*sin_i**2 - 91.0*cos_i**3 + 31.0*cos_i + 60.0)**2*sin_i**8,
        (1, 7) : 35226.59765625*sin_i_half**16*cos_i_half**12,
        (2, 0) : 1268157.515625*sin_i_half**10*cos_i_half**18,
        (2, 1) : 7503.890625*(-91.0*sin_i_half**4 + 65.0*sin_i_half**2 - 10.0)**2*sin_i_half**6*cos_i_half**14,
        (2, 2) : (-0.9228515625*(cos_i + 1.0)**6*sin_i - 2976.75*sin_i_half**9*cos_i_half**5 + 9922.5*sin_i_half**7*cos_i_half**7 - 8505.0*sin_i_half**5*cos_i_half**9 + 2126.25*sin_i_half**3*cos_i_half**11)**2,
        (2, 3) : (0.9228515625*(cos_i + 1.0)**6*sin_i - 1102.5*sin_i_half**11*cos_i_half**3 + 8268.75*sin_i_half**9*cos_i_half**5 - 16537.5*sin_i_half**7*cos_i_half**7 + 11025.0*sin_i_half**5*cos_i_half**9 - 2362.5*sin_i_half**3*cos_i_half**11)**2,
        (2, 4) : 1550.390625*(-343.0*sin_i_half**10 + 620.0*sin_i_half**8 - 280.0*sin_i_half**6 + 420.0*sin_i_half**4*cos_i_half**6 + 238.0*cos_i_half**10 - 210.0*cos_i_half**8)**2*sin_i_half**6*cos_i_half**2,
        (2, 5) : 558.140625*(455.0*sin_i_half**8 - 810.0*sin_i_half**6 + 360.0*sin_i_half**4 + 546.0*cos_i_half**8 - 420.0*cos_i_half**6)**2*sin_i_half**10*cos_i_half**2,
        (2, 6) : 7503.890625*(91.0*sin_i_half**4 - 117.0*sin_i_half**2 + 36.0)**2*sin_i_half**14*cos_i_half**6,
        (2, 7) : 1268157.515625*sin_i_half**18*cos_i_half**10,
        (3, 0) : 31703937.890625*sin_i_half**8*cos_i_half**20,
        (3, 1) : 0.715626776218414306640625*(cos_i + 1.0)**8*(-629.0*cos_i + 338.0*cos_i_double - 91.0*cos_i_triple + 382.0)**2,
        (3, 2) : 13953.515625*(960.0*sin_i_half**8 - 1020.0*sin_i_half**6 + 270.0*sin_i_half**4 + 41.0*cos_i_half**8 - 40.0*cos_i_half**6)**2*cos_i_half**12,
        (3, 3) : 38759.765625*(826.0*sin_i_half**8 - 1176.0*sin_i_half**6 + 420.0*sin_i_half**4 + 175.0*cos_i_half**8 - 160.0*cos_i_half**6)**2*sin_i_half**4*cos_i_half**8,
        (3, 4) : 38759.765625*(-595.0*sin_i_half**10 + 1595.0*sin_i_half**8 - 1420.0*sin_i_half**6 + 420.0*sin_i_half**4 + 406.0*cos_i_half**10 - 336.0*cos_i_half**8)**2*sin_i_half**8,
        (3, 5) : 13953.515625*(311.0*sin_i_half**8 - 580.0*sin_i_half**6 + 270.0*sin_i_half**4 + 690.0*cos_i_half**8 - 480.0*cos_i_half**6)**2*sin_i_half**12,
        (3, 6) : 2931.207275390625*(cos_i + 1.0)**2*(-91.0*sin_i**2 + 78.0*cos_i + 102.0)**2*sin_i_half**16,
        (3, 7) : 31703937.890625*sin_i_half**20*cos_i_half**8,
        (4, 0) : 507263006.25*sin_i_half**6*cos_i_half**22,
        (4, 1) : (-40.60546875*(cos_i + 1.0)**6*sin_i - 95287.5*sin_i_half**5*cos_i_half**9 + 57172.5*sin_i_half**3*cos_i_half**11)**2,
        (4, 2) : (40.60546875*(cos_i + 1.0)**6*sin_i - 155925.0*sin_i_half**7*cos_i_half**7 + 233887.5*sin_i_half**5*cos_i_half**9 - 77962.5*sin_i_half**3*cos_i_half**11)**2,
        (4, 3) : (-744975.0*sin_i_half**6 + 883575.0*sin_i_half**4 - 259875.0*sin_i_half**2 + 43312.5*cos_i_half**6)**2*sin_i_half**6*cos_i_half**10,
        (4, 4) : (-667012.5*sin_i_half**6 + 987525.0*sin_i_half**4 - 363825.0*sin_i_half**2 + 121275.0*cos_i_half**6)**2*sin_i_half**10*cos_i_half**6,
        (4, 5) : 27014006.25*(-61.0*sin_i_half**6 + 105.0*sin_i_half**4 - 45.0*sin_i_half**2 + 30.0*cos_i_half**6)**2*sin_i_half**14*cos_i_half**2,
        (4, 6) : 3001556.25*(91.0*sin_i_half**4 - 143.0*sin_i_half**2 + 55.0)**2*sin_i_half**18*cos_i_half**2,
        (4, 7) : 507263006.25*sin_i_half**22*cos_i_half**6,
        (5, 0) : 4565367056.25*sin_i_half**4*cos_i_half**24,
        (5, 1) : 1648.8040924072265625*(cos_i + 1.0)**10*(-91.0*sin_i**2 - 130.0*cos_i + 134.0)**2,
        (5, 2) : 927.45230197906494140625*(cos_i + 1.0)**8*(-629.0*cos_i + 338.0*cos_i_double - 91.0*cos_i_triple + 382.0)**2,
        (5, 3) : 41220.1023101806640625*(-65.0*sin_i**2 + 91.0*cos_i**3 - 31.0*cos_i + 60.0)**2*sin_i**8,
        (5, 4) : 41220.1023101806640625*(-65.0*sin_i**2 - 91.0*cos_i**3 + 31.0*cos_i + 60.0)**2*sin_i**8,
        (5, 5) : 3798844.62890625*(cos_i + 1.0)**2*(-91.0*sin_i**2 + 78.0*cos_i + 102.0)**2*sin_i_half**16,
        (5, 6) : 1688375.390625*(-91.0*sin_i**2 + 130.0*cos_i + 134.0)**2*sin_i_half**20,
        (5, 7) : 4565367056.25*sin_i_half**24*cos_i_half**4,
        (6, 0) : 1114591.56646728515625*(cos_i + 1.0)**12*sin_i**2,
        (6, 1) : (1055.7421875*(cos_i + 1.0)**6*sin_i - 810810.0*sin_i_half**3*cos_i_half**11)**2,
        (6, 2) : 41088303506.25*(7.0*cos_i - 3.0)**2*sin_i_half**6*cos_i_half**18,
        (6, 3) : 114134176406.25*(7.0*cos_i - 1.0)**2*sin_i_half**10*cos_i_half**14,
        (6, 4) : 114134176406.25*(7.0*cos_i + 1.0)**2*sin_i_half**14*cos_i_half**10,
        (6, 5) : 41088303506.25*(7.0*cos_i + 3.0)**2*sin_i_half**18*cos_i_half**6,
        (6, 6) : 1114591.56646728515625*(cos_i - 1.0)**10*(7.0*cos_i + 5.0)**2*sin_i**2,
        (6, 7) : 18261468225.0*sin_i_half**26*cos_i_half**2,
        (7, 0) : 18261468225.0*cos_i_half**28,
        (7, 1) : 894811943025.0*sin_i_half**4*cos_i_half**24,
        (7, 2) : 8053307487225.0*sin_i_half**8*cos_i_half**20,
        (7, 3) : 22370298575625.0*sin_i_half**12*cos_i_half**16,
        (7, 4) : 22370298575625.0*sin_i_half**16*cos_i_half**12,
        (7, 5) : 8053307487225.0*sin_i_half**20*cos_i_half**8,
        (7, 6) : 894811943025.0*sin_i_half**24*cos_i_half**4,
        (7, 7) : 18261468225.0*sin_i_half**28
    }

    return inclination_results