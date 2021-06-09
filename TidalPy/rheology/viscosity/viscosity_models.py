import numpy as np
from scipy.constants import R

from ...utilities.performance.numba import njit
from ...utilities.types import float_lognat_max


@njit()
def arrhenius(temperature: float, pressure: float,
              arrhenius_coeff: float, additional_temp_dependence: bool, stress: float, stress_expo: float,
              grain_size: float, grain_size_expo: float,
              molar_activation_energy: float, molar_activation_volume: float) -> float:
    """ Solid Viscosity Function: Generic Arrhenius Relationship - NonArray Only

    See, for example, Moore (2006)

    !TPY_args live: self.temperature, self.pressure
    !TPY_args const: arrhenius_coeff, additional_temp_dependence, stress, stress_expo, grain_size, grain_size_expo, molar_activation_energy, molar_activation_volume


    Parameters
    ----------
    temperature : float
        Layer/Material temperature [K]
    pressure : float
        Layer/Material pressure [Pa]
        This could be bottom, surface, or (probably best option) middle or average pressure.
            It should be defined at the same location as temperature.
    arrhenius_coeff : float
        Overall coefficient [s m^(-grain_size) Pa^(stress_expo); also see note in additional_temp_dependence]
    additional_temp_dependence : bool
        If True, the result will be multiplied by temperature
        If True, this will also change the units of the arrhenius coeff to have an additional [K-1]
    stress : float
        Tidal stress [Pa]
    stress_expo : float
        Tidal stress will be raised to (1 - stress_expo)
    grain_size : float
        Material average grain size [m]
    grain_size_expo : float
        Material average grain size will be raised to grain_size_expo
    molar_activation_energy : float
        Defines the pure temperature dependence of the Arrhenius exponential component [J mol-1]
    molar_activation_volume : float
        Defines the pressure dependence of the Arrhenius exponential component [m3 mol-1]

    Returns
    -------
    viscosity : float
        Solid body viscosity of the material (before any partial melting) [Pa s]
    """

    exponent = (molar_activation_energy + pressure * molar_activation_volume) / (temperature * R)

    # Check for over/undershoots
    if exponent > float_lognat_max:
        exponent = float_lognat_max
    elif exponent < -float_lognat_max:
        exponent = -float_lognat_max

    viscosity = arrhenius_coeff * stress**(1. - stress_expo) * (grain_size**grain_size_expo) * np.exp(exponent)
    if additional_temp_dependence:
        viscosity *= temperature

    return viscosity

@njit()
def arrhenius_array(temperature: np.ndarray, pressure: np.ndarray,
                    arrhenius_coeff: float, additional_temp_dependence: bool, stress: float, stress_expo: float,
                    grain_size: float, grain_size_expo: float,
                    molar_activation_energy: float, molar_activation_volume: float) -> np.ndarray:
    """ Solid Viscosity Function: Generic Arrhenius Relationship - Arrays Only

    See, for example, Moore (2006)

    !TPY_args live: self.temperature, self.pressure
    !TPY_args const: arrhenius_coeff, additional_temp_dependence, stress, stress_expo, grain_size, grain_size_expo, molar_activation_energy, molar_activation_volume


    Parameters
    ----------
    temperature : np.ndarray
        Layer/Material temperature [K]
    pressure : np.ndarray
        Layer/Material pressure [Pa]
        This could be bottom, surface, or (probably best option) middle or average pressure.
            It should be defined at the same location as temperature.
    arrhenius_coeff : float
        Overall coefficient [s m^(-grain_size) Pa^(stress_expo); also see note in additional_temp_dependence]
    additional_temp_dependence : bool
        If True, the result will be multiplied by temperature
        If True, this will also change the units of the arrhenius coeff to have an additional [K-1]
    stress : float
        Tidal stress [Pa]
    stress_expo : float
        Tidal stress will be raised to (1 - stress_expo)
    grain_size : float
        Material average grain size [m]
    grain_size_expo : float
        Material average grain size will be raised to grain_size_expo
    molar_activation_energy : float
        Defines the pure temperature dependence of the Arrhenius exponential component [J mol-1]
    molar_activation_volume : float
        Defines the pressure dependence of the Arrhenius exponential component [m3 mol-1]

    Returns
    -------
    viscosity : np.ndarray
        Solid body viscosity of the material (before any partial melting) [Pa s]
    """

    exponent = (molar_activation_energy + pressure * molar_activation_volume) / (temperature * R)

    # Check for over/undershoots
    exponent[exponent > float_lognat_max] = float_lognat_max
    exponent[exponent < -float_lognat_max] = -float_lognat_max

    viscosity = arrhenius_coeff * stress**(1. - stress_expo) * (grain_size**grain_size_expo) * np.exp(exponent)
    if additional_temp_dependence:
        viscosity *= temperature

    return viscosity

@njit()
def reference(temperature: float, pressure: float,
              reference_viscosity: float, reference_temperature: float,
              molar_activation_energy: float, molar_activation_volume: float) -> float:
    """ Solid Viscosity Function: Arrhnius-like utilizing a reference viscosity - NonArrays Only

    See, for example, Henning (2009)

    !TPY_args live: self.temperature, self.pressure
    !TPY_args const: reference_viscosity, reference_temperature, molar_activation_energy, molar_activation_volume

    Parameters
    ----------
    temperature : float
        Layer/Material temperature [K]
    pressure : float
        Layer/Material pressure [Pa]
        This could be bottom, surface, or (probably best option) middle or average pressure.
            It should be defined at the same location as temperature.
    reference_viscosity : float
        Material's viscosity at reference temperature [Pa s]
        Warning: the results of this function will become more inaccurate the further the material's temperature moves
            away from the reference temperature.
    reference_temperature : float
        Defining temperature of the reference viscosity [K]
    molar_activation_energy : float
        Defines the pure temperature dependence of the Arrhenius exponential component [J mol-1]
    molar_activation_volume : float
        Defines the pressure dependence of the Arrhenius exponential component [m3 mol-1]

    Returns
    -------
    viscosity : float
        Solid body viscosity of the material (before any partial melting) [Pa s]
    """

    temp_diff = (1. / temperature) - (1. / reference_temperature)
    exponent = ((molar_activation_energy + pressure * molar_activation_volume) / R) * temp_diff

    if exponent > float_lognat_max:
        exponent = float_lognat_max
    elif exponent < -float_lognat_max:
        exponent = -float_lognat_max

    viscosity = reference_viscosity * np.exp(exponent)

    return viscosity

@njit()
def reference_array(temperature: np.ndarray, pressure: np.ndarray,
                    reference_viscosity: float, reference_temperature: float,
                    molar_activation_energy: float, molar_activation_volume: float) -> np.ndarray:
    """ Solid Viscosity Function: Arrhnius-like utilizing a reference viscosity - Arrays Only

    See, for example, Henning (2009)

    !TPY_args live: self.temperature, self.pressure
    !TPY_args const: reference_viscosity, reference_temperature, molar_activation_energy, molar_activation_volume

    Parameters
    ----------
    temperature : np.ndarray
        Layer/Material temperature [K]
    pressure : np.ndarray
        Layer/Material pressure [Pa]
        This could be bottom, surface, or (probably best option) middle or average pressure.
            It should be defined at the same location as temperature.
    reference_viscosity : float
        Material's viscosity at reference temperature [Pa s]
        Warning: the results of this function will become more inaccurate the further the material's temperature moves
            away from the reference temperature.
    reference_temperature : float
        Defining temperature of the reference viscosity [K]
    molar_activation_energy : float
        Defines the pure temperature dependence of the Arrhenius exponential component [J mol-1]
    molar_activation_volume : float
        Defines the pressure dependence of the Arrhenius exponential component [m3 mol-1]

    Returns
    -------
    viscosity : np.ndarray
        Solid body viscosity of the material (before any partial melting) [Pa s]
    """

    temp_diff = (1. / temperature) - (1. / reference_temperature)
    exponent = ((molar_activation_energy + pressure * molar_activation_volume) / R) * temp_diff

    exponent[exponent > float_lognat_max] = float_lognat_max
    exponent[exponent < -float_lognat_max] = -float_lognat_max

    viscosity = reference_viscosity * np.exp(exponent)

    return viscosity

@njit()
def constant(temperature: float, pressure: float,
             reference_viscosity: float) -> float:
    """ Solid Viscosity Function: Constant. Ignores other input and returns the reference viscosity value - NonArrays Only

    Parameters
    ----------
    temperature : float
        Layer/Material temperature [K]
        FOR THIS FUNCTION: this is unused.
    pressure : float
        Layer/Material pressure [Pa]
        This could be bottom, surface, or (probably best option) middle or average pressure.
            It should be defined at the same location as temperature.
        FOR THIS FUNCTION: this is unused.
    reference_viscosity : float
        Constant viscosity that will be returned by this function [Pa s]

    Returns
    -------
    viscosity : float
        Solid body viscosity of the material (before any partial melting) [Pa s]
    """

    return reference_viscosity

@njit()
def constant_array(temperature: np.ndarray, pressure: np.ndarray,
                   reference_viscosity: float) -> np.ndarray:
    """ Solid Viscosity Function: Constant. Ignores other input and returns the reference viscosity value - Arrays Only

    Parameters
    ----------
    temperature : np.ndarray
        Layer/Material temperature [K]
        FOR THIS FUNCTION: Only used to define the viscosity array's shape
    pressure : np.ndarray
        Layer/Material pressure [Pa]
        This could be bottom, surface, or (probably best option) middle or average pressure.
            It should be defined at the same location as temperature.
        FOR THIS FUNCTION: this is unused.
    reference_viscosity : float
        Constant viscosity that will be returned by this function [Pa s]

    Returns
    -------
    viscosity : np.ndarray
        Solid body viscosity of the material (before any partial melting) [Pa s]
    """

    return reference_viscosity * np.ones_like(temperature)

# Put New Models Here!
