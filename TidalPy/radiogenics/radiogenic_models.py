from typing import Tuple

import numpy as np

from ..utilities.performance.numba import njit
from ..utilities.types import FloatArray

LOG_HALF = np.log(0.5)


@njit(cacheable=True)
def isotope(time: FloatArray, mass: float,
            iso_massfracs_of_isotope: Tuple[float, ...], iso_element_concentrations: Tuple[float, ...],
            iso_halflives: Tuple[float, ...], iso_heat_production: Tuple[float, ...],
            ref_time: float = 4600.) -> FloatArray:
    """ Calculate radiogenic heating based on multiple isotopes

    !TPY_args live: self.time, self.mass
    !TPY_args const: iso_massfracs_of_isotope, iso_element_concentrations, iso_halflives, iso_heat_production, ref_time

    Parameters
    ----------
    time : FloatArray
        Time at which to calculate radiogenic heating at [units must match iso_halflives and ref_time]
    mass : float
        Total mass of radiogenic layer
    iso_massfracs_of_isotope : Tuple[float, ...]
        Mass fraction of isotope in 1 kg of pure element [kg kg-1]
    iso_element_concentrations : Tuple[float, ...]
        Elemental concentration (ppm) at ref_time
    iso_halflives : Tuple[float, ...]
        Isotope half life [units must match time and ref_time]
    iso_heat_production : Tuple[float, ...]
        Isotope heat production rate [Watts kg-1]
    ref_time : float
        Reference time where isotope concentrations were measured [units must match time and iso_halflives]

    Returns
    -------
    radiogenic_heating : FloatArray
        Summed radiogenic heating added for all isotopes [Watts]
    """

    # Start with something fake so that njit knows how to compile. The *0 will make this not matter in the final sum.
    total_specific_heating = time * 0.

    for mass_frac, concentration, halflife, heat_production_rate in \
            zip(iso_massfracs_of_isotope, iso_element_concentrations, iso_halflives, iso_heat_production):
        gamma = LOG_HALF / halflife
        q_iso = mass_frac * concentration * heat_production_rate
        total_specific_heating += q_iso * np.exp(gamma * (time - ref_time))

    # Multiple the specific heating by the total mass (radiogenic mass only)
    radiogenic_heating = total_specific_heating * mass

    return radiogenic_heating


@njit(cacheable=True)
def fixed(time: FloatArray, mass: float,
          fixed_heat_production: float, average_half_life: float,
          ref_time: float = 4600.) -> FloatArray:
    """ Calculate radiogenic heating based on a fixed rate and exponential decay (set at a reference time)

    !TPY_args live: self.time, self.mass
    !TPY_args const: fixed_heat_production, average_half_life, ref_time

    Parameters
    ----------
    time : FloatArray
        Time at which to calculate radiogenic heating at [units must match average_half_life and ref_time]
    mass : float
        Total mass of radiogenic layer
    fixed_heat_production : float
        Fixed heat production rate [Watts kg-1]
    average_half_life : float
        Half life used for the decay of the fixed rate. Set to 0 for no decay [units must match time and ref_time]
    ref_time : float
        Reference time where isotope concentrations were measured [units must match time and iso_halflives]

    Returns
    -------
    radiogenic_heating : FloatArray
        Radiogenic Heating [Watts]

    """

    gamma = LOG_HALF / average_half_life
    radiogenic_heating = mass * fixed_heat_production * np.exp(gamma * (time - ref_time))

    return radiogenic_heating


@njit(cacheable=True)
def off(time: FloatArray, mass: float) -> FloatArray:
    """ Forces radiogenics to be off

    !TPY_args live: self.time, self.mass

    Parameters
    ----------
    time : FloatArray
       Time at which to calculate radiogenic heating at [units must match average_half_life and ref_time]
    mass : float
        Total mass of radiogenic layer

    Returns
    -------
    radiogenic_heating : FloatArray
       Radiogenic heating set to zeros

    """

    shape = 0. * (time + mass)

    # Radiogenic heating is zero. Use the shape as the actual value since it is all zeros
    radiogenic_heating = shape

    return radiogenic_heating

