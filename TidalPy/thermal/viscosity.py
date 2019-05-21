from numba import njit
import numpy as np
from ..types import float_lognat_max
from scipy.constants import R


@njit
def arrhenius(temperature: np.ndarray, pressure: float,
              arrhenius_coeff: float, additional_temp_dependence: bool, stress: float, stress_expo: float,
              grain_size: float, grain_size_expo: float,
              molar_activation_energy: float, molar_activation_volume: float):
    """ Generic Arrhenius Viscosity Function

    See, for example, Moore (2006)

    --- Parameters ---
        other args: arrhenius_coeff, additional_temp_dependence, stress, stress_expo, grain_size, grain_size_expo, molar_activation_energy, molar_activation_volume

    """

    exponent = (molar_activation_energy + pressure*molar_activation_volume) / (temperature * R)

    exponent[exponent > float_lognat_max] = float_lognat_max
    exponent[exponent < -float_lognat_max] = -float_lognat_max

    viscosity = arrhenius_coeff * stress**(1-stress_expo) * grain_size**grain_size_expo * np.exp(exponent)
    if additional_temp_dependence:
        viscosity *= temperature

    return viscosity


@njit
def reference(temperature: np.ndarray, pressure: float,
              reference_viscosity: float, reference_temperature: float,
              molar_activation_energy: float, molar_activation_volume: float):
    """ Arrhenius-like Viscosity Function uses a reference viscosity

    See, for example, Henning (2009)

    --- Parameters ---
        other args: reference_viscosity, reference_temperature, molar_activation_energy, molar_activation_volume

    """

    exponent = ((molar_activation_energy + pressure * molar_activation_volume) / R) * \
               (1. / temperature - 1. / reference_temperature)

    exponent[exponent > float_lognat_max] = float_lognat_max
    exponent[exponent < -float_lognat_max] = -float_lognat_max

    viscosity = reference_viscosity * np.exp(exponent)

    return viscosity

@njit
def constant(temperature: np.ndarray, pressure: float,
             reference_viscosity: float):
    """ Ignores other input and returns the reference viscosity value

    --- Parameters ---
        other args: reference_viscosity

    """

# Put New Models Here!
