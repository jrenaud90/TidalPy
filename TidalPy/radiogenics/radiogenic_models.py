from typing import Tuple

import numpy as np

from ..performance import njit


# Default Radiogenic Parameters

LOG_HALF = np.log(0.5)


@njit
def isotope(time: np.ndarray, mass: float,
            iso_massfracs_of_isotope: Tuple[float], iso_element_concentrations: Tuple[float],
            iso_halflives: Tuple[float], iso_heat_production: Tuple[float],
            ref_time: float = 4600.) -> np.ndarray:
    """ Calculates radiogenics based on multiple isotopes

    other args: iso_massfracs_of_isotope, iso_element_concentrations, iso_halflives, iso_heat_production, ref_time

    :param time:                        <ndarray> Time to calculate radiogenics at
    :param mass:                        <float> Total mass of layer
    :param iso_massfracs_of_isotope:     <tuple of floats> Element concentration (ppm) at ref_time
    :param iso_element_concentrations:  <tuple of floats> Isotope concentration in element (kg kg-1)
    :param iso_halflives:               <tuple of floats> Isotope half lives (same units of time as 'time')
    :param iso_heat_production:         <tuple of floats> Isotope specific heat production in Watts per kg
    :param ref_time:                    <float> Reference time of concentration in same units of time as 'time'
    :return:                            <ndarray> Heating (Watts) for each time in 'time'
    """

    total_specific_heating = np.zeros_like(time)
    for mass_frac, concen, halflife, hpr in \
            zip(iso_massfracs_of_isotope, iso_element_concentrations, iso_halflives, iso_heat_production):
        gamma = LOG_HALF / halflife
        q_iso = mass_frac * concen * hpr
        total_specific_heating += q_iso * np.exp(gamma * (time - ref_time))

    return total_specific_heating * mass


@njit
def fixed(time: np.ndarray, mass: float,
          fixed_heat_production: float, average_half_life: float,
          ref_time: float = 4600.) -> np.ndarray:
    """ Calculates radiogenics based on a fixed rate and fixed exponential decay

    other args: fixed_heat_production, average_half_life, ref_time

    :param time:                    <ndarray> Time to calculate radiogenics at
    :param mass:                    <float> Total mass of layer
    :param fixed_heat_production:   <float> Heat production has Watts kg-1
    :param average_half_life:       <float> Fixed half life in same units of time as 'time'
    :param ref_time:                <float> Reference time of concentration used for fixed_hpr, same units of time as 'time'
    :return:                        <ndarray> Heating (Watts) for each time in 'time'
    """

    gamma = LOG_HALF / average_half_life

    return mass * fixed_heat_production * np.exp(gamma * (time - ref_time))


@njit
def off(time: np.ndarray, mass: float) -> np.ndarray:
    """ Forces radiogenics to be off

    :param time:                    <ndarray> Time to calculate radiogenics at
    :param mass:                    <float> Total mass of layer
    :return:                        <ndarray> Heating (Watts) for each time in 'time'
    """

    return np.zeros_like(time)
