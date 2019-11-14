""" Functions to calculate the frequency dependence of the Andrade 'zeta' parameter.
"""

import numpy as np

from ..performance import tpy_vectorize
from ..types import float_lognat_max


@tpy_vectorize(['float64(float64, float64)'])
def off(baseline_zeta, frequency):
    """ Andrade Parameter Frequency Dependence: Off

    --- Parameters ---
    other args: None

    """

    return baseline_zeta


@tpy_vectorize(['float64(float64, float64, float64)'])
def exponential(baseline_zeta, frequency, andrade_critical_freq):
    """ Andrade Parameter Frequency Dependence: Exponential

    !TPY_args const: andrade_critical_freq

    """

    freq_ratio = andrade_critical_freq / frequency
    # Calculate the frequency ratio to avoid huge values leading to NANs or INFs.
    if freq_ratio > float_lognat_max:
        freq_ratio = float_lognat_max

    zeta = baseline_zeta * np.exp(freq_ratio)
    return zeta


@tpy_vectorize(['float64(float64, float64, float64)'])
def jump(baseline_zeta, frequency, andrade_critical_freq):
    """ Andrade Parameter Frequency Dependence: Jump

    !TPY_args const: andrade_critical_freq

    """

    if frequency >= andrade_critical_freq:
        # Just have it jump to something large - Large zeta makes Andrade --> Maxwell, Sundberg --> Burgers
        zeta = 1.e40
    else:
        zeta = baseline_zeta

    return zeta


@tpy_vectorize(['float64(float64, float64, float64)'])
def elbow(baseline_zeta, frequency, andrade_critical_freq):
    """ Andrade Parameter Frequency Dependence: Elbow

    !TPY_args const: andrade_critical_freq

    """

    freq_ratio = andrade_critical_freq / frequency

    if frequency >= andrade_critical_freq:
        zeta = baseline_zeta * freq_ratio
    else:
        zeta = baseline_zeta

    return zeta
