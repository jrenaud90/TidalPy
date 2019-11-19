""" Functions to calculate the frequency dependence of the Andrade 'zeta' parameter.
"""

import numpy as np

from TidalPy.performance import njit
from TidalPy.types import float_lognat_max


@njit
def off(baseline_zeta, frequency):
    """ Andrade Parameter Frequency Dependence: Off

    --- Parameters ---
    other args: None

    """

    return baseline_zeta


@njit
def exponential(baseline_zeta, frequency, critical_freq):
    """ Andrade Parameter Frequency Dependence: Exponential

    !TPY_args const: critical_freq

    """

    # Calculate the frequency ratio to avoid huge values leading to NANs or INFs.
    _freq_ratio = critical_freq / frequency
    _freq_ratio[_freq_ratio > float_lognat_max] = float_lognat_max

    zeta = baseline_zeta * np.exp(_freq_ratio)
    return zeta


@njit
def jump(baseline_zeta, frequency, critical_freq):
    """ Andrade Parameter Frequency Dependence: Jump

    !TPY_args const: critical_freq

    """

    zeta = baseline_zeta
    zeta[frequency >= critical_freq] = 1.e40

    return zeta


@njit
def elbow(baseline_zeta, frequency, critical_freq):
    """ Andrade Parameter Frequency Dependence: Elbow

    !TPY_args const: critical_freq

    """

    _freq_ratio = critical_freq / frequency

    zeta = baseline_zeta * _freq_ratio
    zeta[frequency >= critical_freq] = baseline_zeta

    return zeta
