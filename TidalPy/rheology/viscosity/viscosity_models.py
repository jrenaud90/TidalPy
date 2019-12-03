import numpy as np
from scipy.constants import R

from TidalPy.performance import njit
from TidalPy.types import float_lognat_max


@njit
def arrhenius(temperature: np.ndarray, pressure: np.ndarray,
              arrhenius_coeff: float, additional_temp_dependence: bool, stress: float, stress_expo: float,
              grain_size: float, grain_size_expo: float,
              molar_activation_energy: float, molar_activation_volume: float) -> np.ndarray:
    """ Generic Arrhenius Viscosity Function

    See, for example, Moore (2006)

    !TPY_args live: self.temperature, self.pressure
    !TPY_args const: arrhenius_coeff, additional_temp_dependence, stress, stress_expo, grain_size, grain_size_expo, molar_activation_energy, molar_activation_volume

    """

    exponent = (molar_activation_energy + pressure * molar_activation_volume) / (temperature * R)

    exponent[exponent > float_lognat_max] = float_lognat_max
    exponent[exponent < -float_lognat_max] = -float_lognat_max

    viscosity = arrhenius_coeff * stress**(1 - stress_expo) * (grain_size**grain_size_expo) * np.exp(exponent)
    if additional_temp_dependence:
        viscosity *= temperature

    return viscosity


@njit
def reference(temperature: np.ndarray, pressure: np.ndarray,
              reference_viscosity: float, reference_temperature: float,
              molar_activation_energy: float, molar_activation_volume: float) -> np.ndarray:
    """ Arrhenius-like Viscosity Function uses a reference viscosity

    See, for example, Henning (2009)

    !TPY_args live: self.temperature, self.pressure
    !TPY_args const: reference_viscosity, reference_temperature, molar_activation_energy, molar_activation_volume

    """

    temp_diff = (1. / temperature) - (1. / reference_temperature)
    exponent = ((molar_activation_energy + pressure * molar_activation_volume) / R) * temp_diff

    exponent[exponent > float_lognat_max] = float_lognat_max
    exponent[exponent < -float_lognat_max] = -float_lognat_max

    viscosity = reference_viscosity * np.exp(exponent)

    return viscosity


@njit
def constant(temperature: np.ndarray, pressure: np.ndarray,
             reference_viscosity: float) -> np.ndarray:
    """ Ignores other input and returns the reference viscosity value

    !TPY_args live: self.temperature, self.pressure
    !TPY_args const: reference_viscosity

    """

    return reference_viscosity * np.ones_like(temperature)

# Put New Models Here!
