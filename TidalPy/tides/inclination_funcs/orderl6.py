""" Inclination functions (squared) for tidal order-l = 6. These are exact (no truncation on I)
"""

from typing import TYPE_CHECKING

import numpy as np

from . import InclinOutput
from ...utilities.performance.numba import njit

if TYPE_CHECKING:
    from ...utilities.types import FloatArray


@njit(cacheable=True)
def calc_inclination_off(inclination: 'FloatArray') -> 'InclinOutput':
    """Calculate F^2_lmp (assuming I=0) for l = 6"""

    # Inclination Functions Calculated for l = 6, Inclination == off.
    ones_ = np.ones_like(inclination)

    inclination_results = {
        (0, 3): 0.09765625 * ones_,
        (2, 2): 172.265625 * ones_,
        (4, 1): 223256.25 * ones_,
        (6, 0): 108056025. * ones_,
        }

    return inclination_results


@njit(cacheable=True)
def calc_inclination(inclination: 'FloatArray') -> 'InclinOutput':
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
        (0, 1) : 1.832004547119140625*(sin_i**2 - 0.9090909090909090909090909)**2*sin_i**8,
        (0, 2) : 165547.265625*(0.3870967741935483870967742*sin_i_half**10 - sin_i_half**8 + 0.8548387096774193548387097*sin_i_half**6 - 0.2419354838709677419354839*sin_i_half**4 - 0.1451612903225806451612903*cos_i_half**10 + 0.1290322580645161290322581*cos_i_half**8)**2*sin_i_half**4,
        (0, 3) : 23066.015625*(0.539094650205761316872428*sin_i_half**12 - sin_i_half**10 + 0.462962962962962962962963*sin_i_half**8 - 0.8230452674897119341563786*sin_i_half**6*cos_i_half**6 + 0.462962962962962962962963*sin_i_half**4*cos_i_half**8 + 0.07613168724279835390946502*cos_i_half**12 - 0.07407407407407407407407407*cos_i_half**10)**2,
        (0, 4) : 165547.265625*(0.3870967741935483870967742*sin_i_half**10 - sin_i_half**8 + 0.8548387096774193548387097*sin_i_half**6 - 0.2419354838709677419354839*sin_i_half**4 - 0.1451612903225806451612903*cos_i_half**10 + 0.1290322580645161290322581*cos_i_half**8)**2*sin_i_half**4,
        (0, 5) : 1.832004547119140625*(sin_i**2 - 0.9090909090909090909090909)**2*sin_i**8,
        (0, 6) : 208.44140625*sin_i_half**12*cos_i_half**12,
        (1, 0) : 7503.890625*sin_i_half**10*cos_i_half**14,
        (1, 1) : 270140.0625*(sin_i_half**4 - 0.8333333333333333333333333*sin_i_half**2 + 0.1515151515151515151515152)**2*sin_i_half**6*cos_i_half**10,
        (1, 2) : 303876.5625*(0.0003720238095238095238095238*(cos_i + 1.0)**5*sin_i + 0.1666666666666666666666667*sin_i_half**9*cos_i_half**3 - 0.8333333333333333333333333*sin_i_half**7*cos_i_half**5 + sin_i_half**5*cos_i_half**7 - 0.3333333333333333333333333*sin_i_half**3*cos_i_half**9)**2,
        (1, 3) : 430664.0625*(-0.0003125*(cos_i + 1.0)**5*sin_i + 0.02*sin_i_half**11*cos_i_half - 0.3*sin_i_half**9*cos_i_half**3 + sin_i_half**7*cos_i_half**5 - sin_i_half**5*cos_i_half**7 + 0.3*sin_i_half**3*cos_i_half**9)**2,
        (1, 4) : 1654439.0625*(-0.5816326530612244897959184*sin_i_half**8 + sin_i_half**6 - 0.4285714285714285714285714*sin_i_half**4 - 0.4285714285714285714285714*cos_i_half**8 + 0.3571428571428571428571429*cos_i_half**6)**2*sin_i_half**6*cos_i_half**2,
        (1, 5) : 367690.640625*(-0.8571428571428571428571429*sin_i_half**4 + sin_i_half**2 - 0.2727272727272727272727273)**2*sin_i_half**10*cos_i_half**6,
        (1, 6) : 7503.890625*sin_i_half**14*cos_i_half**10,
        (2, 0) : 187597.265625*sin_i_half**8*cos_i_half**16,
        (2, 1) : 3452.13625431060791015625*(cos_i + 1.0)**6*(-cos_i + 0.5759162303664921465968586*cos_i_double - 0.1727748691099476439790576*cos_i_triple + 0.596858638743455497382199)**2,
        (2, 2) : 54022500.0*(0.075*(cos_i - 1.0)**2 + 0.825*sin_i_half**8 - sin_i_half**6 + 0.05892857142857142857142857*cos_i_half**8 - 0.05714285714285714285714286*cos_i_half**6)**2*cos_i_half**8,
        (2, 3) : 264875625.0*(-0.3870967741935483870967742*sin_i_half**10 + sin_i_half**8 - 0.8548387096774193548387097*sin_i_half**6 + 0.2419354838709677419354839*sin_i_half**4 + 0.1451612903225806451612903*cos_i_half**10 - 0.1290322580645161290322581*cos_i_half**8)**2*sin_i_half**4,
        (2, 4) : 23328900.0*(0.5461956521739130434782609*sin_i_half**8 - sin_i_half**6 + 0.4565217391304347826086956*sin_i_half**4 + 0.7989130434782608695652174*cos_i_half**8 - 0.6086956521739130434782609*cos_i_half**6)**2*sin_i_half**8,
        (2, 5) : 112015.72265625*(cos_i + 1.0)**2*(-0.9705882352941176470588235*sin_i**2 + 0.6470588235294117647058823*cos_i + 1)**2*sin_i_half**12,
        (2, 6) : 187597.265625*sin_i_half**16*cos_i_half**8,
        (3, 0) : 3001556.25*sin_i_half**6*cos_i_half**18,
        (3, 1) : 32148900.0*(-0.001302083333333333333333333*(cos_i + 1.0)**5*sin_i - sin_i_half**5*cos_i_half**7 + 0.75*sin_i_half**3*cos_i_half**9)**2,
        (3, 2) : 175032900.0*(0.0005580357142857142857142857*(cos_i + 1.0)**5*sin_i - 0.5*sin_i_half**7*cos_i_half**5 + sin_i_half**5*cos_i_half**7 - 0.4285714285714285714285714*sin_i_half**3*cos_i_half**9)**2,
        (3, 3) : 18759726.5625*(sin_i_half**4 - sin_i_half**2 + 0.1818181818181818181818182)**2*sin_i_half**4*sin_i_double**2*cos_i_half,
        (3, 4) : 1032336900.0*(-0.6029411764705882352941177*sin_i_half**6 + sin_i_half**4 - 0.4117647058823529411764706*sin_i_half**2 + 0.2058823529411764705882353*cos_i_half**6)**2*sin_i_half**10*cos_i_half**2,
        (3, 5) : 243126056.25*(0.6666666666666666666666667*sin_i_half**4 - sin_i_half**2 + 0.3636363636363636363636364)**2*sin_i_half**14*cos_i_half**2,
        (3, 6) : 3001556.25*sin_i_half**18*cos_i_half**6,
        (4, 0) : 27014006.25*sin_i_half**4*cos_i_half**20,
        (4, 1) : 461338.1103515625*(cos_i + 1.0)**8*(0.7173913043478260869565217*sin_i**2 + 0.9565217391304347826086957*cos_i - 1)**2,
        (4, 2) : 3106922.628879547119140625*(cos_i + 1.0)**6*(cos_i - 0.5759162303664921465968586*cos_i_double + 0.1727748691099476439790576*cos_i_triple - 0.596858638743455497382199)**2,
        (4, 3) : 675350156.25*(sin_i**2 - 0.9090909090909090909090909)**2*sin_i_half**8*cos_i_half**8,
        (4, 4) : 100814150.390625*(cos_i + 1.0)**2*(-0.9705882352941176470588235*sin_i**2 + 0.6470588235294117647058823*cos_i + 1)**2*sin_i_half**12,
        (4, 5) : 118102556.25*(0.7173913043478260869565217*sin_i**2 - 0.9565217391304347826086957*cos_i - 1)**2*sin_i_half**16,
        (4, 6) : 27014006.25*sin_i_half**20*cos_i_half**4,
        (5, 0) : 26380.865478515625*(cos_i + 1.0)**10*sin_i**2,
        (5, 1) : 2701400625.0*(-0.003125*(cos_i + 1.0)**5*sin_i + sin_i_half**3*cos_i_half**9)**2,
        (5, 2) : 6078151406.25*(0.3333333333333333333333333 - cos_i)**2*sin_i_half**6*cos_i_half**14,
        (5, 3) : 10805602500.0*sin_i_half**10*cos_i_half**10*cos_i**2,
        (5, 4) : 24312605625.0*(0.3333333333333333333333333 - cos_i_half**2)**2*sin_i_half**14*cos_i_half**6,
        (5, 5) : 949711.1572265625*(cos_i - 1.0)**8*(cos_i + 0.6666666666666666666666667)**2*sin_i**2,
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
