""" Inclination functions (squared) for tidal order-l = 5. These are exact (no truncation on I)
"""

from typing import TYPE_CHECKING

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def calc_inclination_off(inclination: 'FloatArray') -> 'InclinOutput':
    """Calculate F^2_lmp (assuming I=0) for l = 5"""

    # Inclination Functions Calculated for l = 5, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (1, 2): 3.515625 * ones_,
        (3, 1): 2756.25 * ones_,
        (5, 0): 893025. * ones_,
        }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: 'FloatArray') -> 'InclinOutput':
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
        (0, 1) : 1.5140533447265625*(0.8888888888888888888888889 - sin_i**2)**2*sin_i**6,
        (0, 2) : 1406.25*(0.0015625*(cos_i + 1.0)**4*sin_i + 0.05*sin_i_half**9*cos_i_half - 0.5*sin_i_half**7*cos_i_half**3 + sin_i_half**5*cos_i_half**5 - 0.5*sin_i_half**3*cos_i_half**7)**2,
        (0, 3) : 1406.25*(-0.0015625*(cos_i + 1.0)**4*sin_i - 0.05*sin_i_half**9*cos_i_half + 0.5*sin_i_half**7*cos_i_half**3 - sin_i_half**5*cos_i_half**5 + 0.5*sin_i_half**3*cos_i_half**7)**2,
        (0, 4) : 1.5140533447265625*(sin_i**2 - 0.8888888888888888888888889)**2*sin_i**6,
        (0, 5) : 62.015625*sin_i_half**10*cos_i_half**10,
        (1, 0) : 1550.390625*sin_i_half**8*cos_i_half**12,
        (1, 1) : 44.42274570465087890625*(cos_i + 1.0)**4*(-cos_i + 0.6461538461538461538461539*cos_i_double - 0.2307692307692307692307692*cos_i_triple + 0.5846153846153846153846154)**2,
        (1, 2) : 696181.640625*(-0.4157303370786516853932584*sin_i_half**10 + sin_i_half**8 - 0.7865168539325842696629214*sin_i_half**6 + 0.2022471910112359550561798*sin_i_half**4 + 0.05617977528089887640449438*cos_i_half**10 - 0.05393258426966292134831461*cos_i_half**8)**2,
        (1, 3) : 146306.25*(0.563725490196078431372549*sin_i_half**8 - sin_i_half**6 + 0.4411764705882352941176471*sin_i_half**4 + 0.4656862745098039215686275*cos_i_half**8 - 0.3921568627450980392156863*cos_i_half**6)**2*sin_i_half**4,
        (1, 4) : 605.621337890625*(cos_i + 1.0)**2*(-sin_i**2 + 0.4*cos_i + 0.9333333333333333333333333)**2*sin_i_half**8,
        (1, 5) : 1550.390625*sin_i_half**12*cos_i_half**8,
        (2, 0) : 24806.25*sin_i_half**6*cos_i_half**14,
        (2, 1) : 135056.25*(-0.004464285714285714285714286*(cos_i + 1.0)**4*sin_i - sin_i_half**5*cos_i_half**5 + sin_i_half**3*cos_i_half**7)**2,
        (2, 2) : 620156.25*(0.002083333333333333333333333*(cos_i + 1.0)**4*sin_i - 0.3333333333333333333333333*sin_i_half**7*cos_i_half**3 + sin_i_half**5*cos_i_half**5 - 0.6*sin_i_half**3*cos_i_half**7)**2,
        (2, 3) : 4192256.25*(-0.6410256410256410256410256*sin_i_half**6 + sin_i_half**4 - 0.3846153846153846153846154*sin_i_half**2 + 0.1282051282051282051282051*cos_i_half**6)**2*sin_i_half**6*cos_i_half**2,
        (2, 4) : 1215506.25*(0.7142857142857142857142857*sin_i_half**4 - sin_i_half**2 + 0.3333333333333333333333333)**2*sin_i_half**10*cos_i_half**2,
        (2, 5) : 24806.25*sin_i_half**14*cos_i_half**6,
        (3, 0) : 223256.25*sin_i_half**4*cos_i_half**16,
        (3, 1) : 7848.8525390625*(cos_i + 1.0)**6*(-0.8333333333333333333333333*cos_i**2 + cos_i - 0.2407407407407407407407407)**2,
        (3, 2) : 25587.50152587890625*(cos_i + 1.0)**4*(cos_i - 0.6461538461538461538461539*cos_i_double + 0.2307692307692307692307692*cos_i_triple - 0.5846153846153846153846154)**2,
        (3, 3) : 348837.890625*(cos_i + 1.0)**2*(-sin_i**2 + 0.4*cos_i + 0.9333333333333333333333333)**2*sin_i_half**8,
        (3, 4) : 579501.5625*(0.7758620689655172413793103*sin_i**2 - 0.9310344827586206896551724*cos_i - 1)**2*sin_i_half**12,
        (3, 5) : 223256.25*sin_i_half**16*cos_i_half**4,
        (4, 0) : 872.0947265625*(cos_i + 1.0)**8*sin_i**2,
        (4, 1) : 14288400.0*(-0.0078125*(cos_i + 1.0)**4*sin_i + sin_i_half**3*cos_i_half**7)**2,
        (4, 2) : 22325625.0*(0.2 - cos_i)**2*sin_i_half**6*cos_i_half**10,
        (4, 3) : 22325625.0*(-cos_i - 0.2)**2*sin_i_half**10*cos_i_half**6,
        (4, 4) : 5581406.25*(cos_i + 0.6)**2*sin_i_half**14*cos_i_half**2,
        (4, 5) : 893025.0*sin_i_half**18*cos_i_half**2,
        (5, 0) : 893025.0*cos_i_half**20,
        (5, 1) : 22325625.0*sin_i_half**4*cos_i_half**16,
        (5, 2) : 89302500.0*sin_i_half**8*cos_i_half**12,
        (5, 3) : 89302500.0*sin_i_half**12*cos_i_half**8,
        (5, 4) : 22325625.0*sin_i_half**16*cos_i_half**4,
        (5, 5) : 893025.0*sin_i_half**20
    }

    return inclination_results

