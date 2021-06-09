""" Partial Melting Models for Viscosity and Shear Modulus """

from typing import Tuple

import numpy as np

from ...utilities.performance.numba import njit


@njit()
def off(melt_fraction: float, premelt_viscosity: float, premelt_shear: float) -> Tuple[float, float]:
    """ Viscosity and Shear Modulus Partial Melting Model: off - NonArrays Only

    !TPY_args live: self.premelt_viscosity, self.premelt_shear

    Parameters
    ----------
    melt_fraction : float
        Layer/Material volumetric melt fraction [m3 m-3]
    premelt_viscosity : float
        Layer/Material viscosity before partial melting is considered [Pa s]
    premelt_shear : float
        Layer/Material shear modulus before partial melting is considered [Pa]

    Returns
    -------
    postmelt_viscosity : float
        Post melting viscosity
    postmelt_shear_modulus : float
        Post melting shear modulus
    """

    postmelt_viscosity = premelt_viscosity
    postmelt_shear_modulus = premelt_shear

    return postmelt_viscosity, postmelt_shear_modulus

@njit()
def off_array(melt_fraction: np.ndarray, premelt_viscosity: np.ndarray, premelt_shear: float) -> \
        Tuple[np.ndarray, np.ndarray]:
    """ Viscosity and Shear Modulus Partial Melting Model: off - Arrays Only

    !TPY_args live: self.premelt_viscosity, self.premelt_shear

    Parameters
    ----------
    melt_fraction : np.ndarray
        Layer/Material volumetric melt fraction [m3 m-3]
    premelt_viscosity : np.ndarray
        Layer/Material viscosity before partial melting is considered [Pa s]
    premelt_shear : np.ndarray
        Layer/Material shear modulus before partial melting is considered [Pa]

    Returns
    -------
    postmelt_viscosity : np.ndarray
        Post melting viscosity
    postmelt_shear_modulus : np.ndarray
        Post melting shear modulus
    """

    postmelt_viscosity = premelt_viscosity
    postmelt_shear_modulus = premelt_shear * np.ones_like(premelt_viscosity)

    return postmelt_viscosity, postmelt_shear_modulus

@njit()
def spohn(melt_fraction: float, temperature: float, liquid_viscosity: float,
          liquid_shear: float = 1.0e-5,
          fs_visc_power_slope: float = 27000.0, fs_visc_power_phase: float = 1.0,
          fs_shear_power_slope: float = 82000.0, fs_shear_power_phase: float = 40.6) -> Tuple[float, float]:
    """ Viscosity and Shear Modulus Partial Melting Model: spohn - NonArrays Only

    Fischer and Spohn (1990) Partial-Melt Viscosity Function

    !TPY_args live: self.temperature, self.liquid_viscosity
    !TPY_args const: liquid_shear, fs_visc_power_slope, fs_visc_power_phase, fs_shear_power_slope, fs_shear_power_phase

    Parameters
    ----------
    melt_fraction : float
        Layer/Material volumetric melt fraction [m3 m-3]
    temperature : float
        Layer/Material temperature [K]
    liquid_viscosity : float
        Layer/Material viscosity if it were completely molten at this temperature [Pa s]
    liquid_shear : float
        Material's shear modulus assuming pure liquid [Pa]
    fs_visc_power_slope : float
        Fischer & Spohn viscosity exponent multiplier parameter [K]
    fs_visc_power_phase : float
        Fischer & Spohn viscosity exponent additive parameter [K]
    fs_shear_power_slope : float
        Fischer & Spohn viscosity exponent multiplier parameter [K]
    fs_shear_power_phase : float
        Fischer & Spohn viscosity exponent additive parameter [K]

    Returns
    -------
    postmelt_viscosity : float
        Post melting viscosity
    postmelt_shear_modulus : float
        Post melting shear modulus
    """

    postmelt_viscosity = 10.**((fs_visc_power_slope / temperature) - fs_visc_power_phase)
    postmelt_shear_modulus = 10.**((fs_shear_power_slope / temperature) - fs_shear_power_phase)

    # Perform sanity checks
    if postmelt_viscosity < liquid_viscosity:
        postmelt_viscosity = liquid_viscosity
    if postmelt_shear_modulus < liquid_shear:
        postmelt_shear_modulus = liquid_shear

    return postmelt_viscosity, postmelt_shear_modulus

@njit()
def spohn_array(melt_fraction: np.ndarray, temperature: np.ndarray, liquid_viscosity: np.ndarray,
                liquid_shear: float = 1.0e-5,
                fs_visc_power_slope: float = 27000.0, fs_visc_power_phase: float = 1.0,
                fs_shear_power_slope: float = 82000.0, fs_shear_power_phase: float = 40.6) -> \
        Tuple[np.ndarray, np.ndarray]:
    """ Viscosity and Shear Modulus Partial Melting Model: spohn - Arrays Only

    Fischer and Spohn (1990) Partial-Melt Viscosity Function

    !TPY_args live: self.temperature, self.liquid_viscosity
    !TPY_args const: liquid_shear, fs_visc_power_slope, fs_visc_power_phase, fs_shear_power_slope, fs_shear_power_phase

    Parameters
    ----------
    melt_fraction : np.ndarray
        Layer/Material volumetric melt fraction [m3 m-3]
    temperature : np.ndarray
        Layer/Material temperature [K]
    liquid_viscosity : np.ndarray
        Layer/Material viscosity if it were completely molten at this temperature [Pa s]
    liquid_shear : float
        Material's shear modulus assuming pure liquid [Pa]
    fs_visc_power_slope : float
        Fischer & Spohn viscosity exponent multiplier parameter [K]
    fs_visc_power_phase : float
        Fischer & Spohn viscosity exponent additive parameter [K]
    fs_shear_power_slope : float
        Fischer & Spohn viscosity exponent multiplier parameter [K]
    fs_shear_power_phase : float
        Fischer & Spohn viscosity exponent additive parameter [K]

    Returns
    -------
    postmelt_viscosity : np.ndarray
        Post melting viscosity
    postmelt_shear_modulus : np.ndarray
        Post melting shear modulus
    """

    postmelt_viscosity = 10.**((fs_visc_power_slope / temperature) - fs_visc_power_phase)
    postmelt_shear_modulus = 10.**((fs_shear_power_slope / temperature) - fs_shear_power_phase)

    # Perform sanity checks
    postmelt_viscosity[postmelt_viscosity < liquid_viscosity] = liquid_viscosity[postmelt_viscosity < liquid_viscosity]
    postmelt_shear_modulus[postmelt_shear_modulus < liquid_shear] = liquid_shear

    return postmelt_viscosity, postmelt_shear_modulus

@njit()
def henning(melt_fraction: float, temperature: float,
            premelt_viscosity: float, liquid_viscosity: float, premelt_shear: float, solidus: float, liquidus: float,
            liquid_shear: float, crit_melt_frac: float = 0.5, crit_melt_frac_width: float = 0.05,
            hn_visc_slope_1: float = 13.5, hn_visc_falloff_slope: float = 370., hn_shear_param_1: float = 40000.,
            hn_shear_param_2: float = 25., hn_shear_falloff_slope: float = 700.) -> Tuple[float, float]:
    """ Viscosity and Shear Modulus Partial Melting Model: henning - NonArrays Only

    Henning (2009, 2010) Partial-Melt Viscosity Function
    See also Moore's work and Renaud and Henning (2018)

    !TPY_args live: self.temperature, self.premelt_viscosity, self.liquid_viscosity, self.premelt_shear, self.solidus, self.liquidus
    !TPY_args const: liquid_shear, crit_melt_frac, crit_melt_frac_width, hn_visc_slope_1, hn_visc_slope_2, hn_shear_param_1, hn_shear_param_2, hn_shear_falloff_slope

    Parameters
    ----------
    melt_fraction : float
        Layer/Material volumetric melt fraction [m3 m-3]
    temperature : float
        Layer/Material temperature [K]
    premelt_viscosity : float
        Layer/Material viscosity before partial melting is considered [Pa s]
    premelt_shear : float
        Layer/Material shear modulus before partial melting is considered [Pa]
    solidus : float
        Layer/Material solidus temperature
    liquidus : float
        Layer/Material liquidus temperature
    liquid_viscosity : float
        Layer/Material viscosity if it were completely molten at this temperature [Pa s]
    liquid_shear : float
        Material's shear modulus assuming pure liquid [Pa]
    crit_melt_frac : float
        Melt Fraction where the material behaves more like a liquid than a solid [m3 m-3]
    crit_melt_frac_width : float
        Defines the partial melt transition zone [m3 m-3]
        Zone crit_melt_frac and crit_melt_frac + crit_melt_frac_width defines the transition between solid-like and
            liquid-like responses.
    hn_visc_slope_1 : float
        Henning, pre-breakdown, viscosity exponent multiplier parameter
    hn_visc_falloff_slope : float
        Henning, breakdown, viscosity exponent multiplier parameter
    hn_shear_param_1 : float
        Henning, pre-breakdown, shear modulus exponent multiplier parameter 1 [K]
    hn_shear_param_2 : float
        Henning, pre-breakdown, shear modulus exponent multiplier parameter 2
    hn_shear_falloff_slope : float
        Henning, breakdown, shear modulus exponent multiplier parameter

    Returns
    -------
    postmelt_viscosity : float
        Post melting viscosity
    postmelt_shear_modulus : float
        Post melting shear modulus
    """

    crit_melt_frac_plus_width = crit_melt_frac + crit_melt_frac_width
    break_down_temp = solidus + crit_melt_frac * (liquidus - solidus)

    # Initialize post-melt with pre-melt values (if there is no melt fraction then these values will be returned.
    postmelt_viscosity = premelt_viscosity
    postmelt_shear_modulus = premelt_shear

    # Calculate viscosity and shear modulus in the three domains
    if melt_fraction <= 0.:
        # No partial melt.
        pass
    else:
        # Partial melting
        if melt_fraction < crit_melt_frac:
            # Partial melting before critical break-down.
            postmelt_viscosity *= np.exp(-hn_visc_slope_1 * melt_fraction)
            postmelt_shear_modulus *= np.exp((hn_shear_param_1/temperature) - hn_shear_param_2)
        elif crit_melt_frac <= melt_fraction <= crit_melt_frac_plus_width:
            # Get the maximum pre-melt effect
            postmelt_viscosity *= np.exp(-hn_visc_slope_1 * crit_melt_frac)
            postmelt_shear_modulus *= np.exp((hn_shear_param_1/break_down_temp) - hn_shear_param_2)
            # Then apply critical breakdown occurring.
            postmelt_viscosity *= np.exp(-hn_visc_falloff_slope * (melt_fraction - crit_melt_frac))
            postmelt_shear_modulus *= np.exp(-hn_shear_falloff_slope * (melt_fraction - crit_melt_frac))
        elif melt_fraction > crit_melt_frac_plus_width:
            # Past critical breakdown threshold. The material now behaves like a fluid.
            postmelt_viscosity = liquid_viscosity
            postmelt_shear_modulus = liquid_shear

    # Perform sanity corrections
    if postmelt_viscosity < liquid_viscosity:
        postmelt_viscosity = liquid_viscosity
    if postmelt_shear_modulus < liquid_shear:
        postmelt_shear_modulus = liquid_shear

    return postmelt_viscosity, postmelt_shear_modulus

@njit()
def henning_array(melt_fraction: np.ndarray, temperature: np.ndarray,
                  premelt_viscosity: np.ndarray, liquid_viscosity: np.ndarray,
                  premelt_shear: float, solidus: float, liquidus: float,
                  liquid_shear: float, crit_melt_frac: float = 0.5, crit_melt_frac_width: float = 0.05,
                  hn_visc_slope_1: float = 13.5, hn_visc_falloff_slope: float = 370., hn_shear_param_1: float = 40000.,
                  hn_shear_param_2: float = 25., hn_shear_falloff_slope: float = 700.) -> Tuple[np.ndarray, np.ndarray]:
    """ Viscosity and Shear Modulus Partial Melting Model: henning

    Henning (2009, 2010) Partial-Melt Viscosity Function
    See also Moore's work and Renaud and Henning (2018)

    !TPY_args live: self.temperature, self.premelt_viscosity, self.liquid_viscosity, self.premelt_shear
    !TPY_args const: liquid_shear, crit_melt_frac, crit_melt_frac_width, hn_visc_slope_1, hn_visc_slope_2, hn_shear_param_1, hn_shear_param_2, hn_shear_falloff_slope

    Parameters
    ----------
    melt_fraction : np.ndarray
        Layer/Material volumetric melt fraction [m3 m-3]
    temperature : np.ndarray
        Layer/Material temperature [K]
    premelt_viscosity : np.ndarray
        Layer/Material viscosity before partial melting is considered [Pa s]
    premelt_shear : np.ndarray
        Layer/Material shear modulus before partial melting is considered [Pa]
    solidus : float
        Layer/Material solidus temperature
    liquidus : float
        Layer/Material liquidus temperature
    liquid_viscosity : np.ndarray
        Layer/Material viscosity if it were completely molten at this temperature [Pa s]
    liquid_shear : float
        Material's shear modulus assuming pure liquid [Pa]
    crit_melt_frac : float
        Melt Fraction where the material behaves more like a liquid than a solid [m3 m-3]
    crit_melt_frac_width : float
        Defines the partial melt transition zone [m3 m-3]
        Zone crit_melt_frac and crit_melt_frac + crit_melt_frac_width defines the transition between solid-like and
            liquid-like responses.
    hn_visc_slope_1 : float
        Henning, pre-breakdown, viscosity exponent multiplier parameter
    hn_visc_falloff_slope : float
        Henning, breakdown, viscosity exponent multiplier parameter
    hn_shear_param_1 : float
        Henning, pre-breakdown, shear modulus exponent multiplier parameter 1 [K]
    hn_shear_param_2 : float
        Henning, pre-breakdown, shear modulus exponent multiplier parameter 2
    hn_shear_falloff_slope : float
        Henning, breakdown, shear modulus exponent multiplier parameter

    Returns
    -------
    postmelt_viscosity : np.ndarray
        Post melting viscosity
    postmelt_shear_modulus : np.ndarray
        Post melting shear modulus
    """

    crit_melt_frac_plus_width = crit_melt_frac + crit_melt_frac_width
    break_down_temp = solidus + crit_melt_frac * (liquidus - solidus)

    # Calculate ndarray indices that define the three domains
    pre_breakdown_index = np.logical_and(melt_fraction < crit_melt_frac, melt_fraction > 0)
    breakdown_index = np.logical_and(melt_fraction >= crit_melt_frac, melt_fraction <= crit_melt_frac_plus_width)
    molten_index = melt_fraction > crit_melt_frac_plus_width

    # Calculate viscosity in the three domains
    postmelt_viscosity = premelt_viscosity
    postmelt_viscosity[pre_breakdown_index] *= np.exp(-hn_visc_slope_1 * melt_fraction[pre_breakdown_index])
    # For breakdown_index, apply maximum from the pre-critical domain
    postmelt_viscosity[breakdown_index] *= np.exp(-hn_visc_slope_1 * crit_melt_frac)
    # Then apply the breakdown effects
    postmelt_viscosity[breakdown_index] *= np.exp(-hn_visc_falloff_slope * (melt_fraction[breakdown_index] - crit_melt_frac))
    postmelt_viscosity[molten_index] = liquid_viscosity[molten_index]

    # Calculate shear modulus in the three domains
    postmelt_shear_modulus = premelt_shear * np.ones_like(temperature)
    postmelt_shear_modulus[pre_breakdown_index] *= np.exp((hn_shear_param_1 / temperature[pre_breakdown_index]) - hn_shear_param_2)
    # For breakdown_index, apply maximum from the pre-critical domain
    postmelt_shear_modulus[breakdown_index] *= np.exp((hn_shear_param_1 / break_down_temp) - hn_shear_param_2)
    # Then apply the breakdown effects
    postmelt_shear_modulus[breakdown_index] *= np.exp(-hn_shear_falloff_slope * (melt_fraction[breakdown_index] - crit_melt_frac))
    postmelt_shear_modulus[molten_index] = liquid_shear

    # Perform sanity checks
    postmelt_viscosity[postmelt_viscosity < liquid_viscosity] = liquid_viscosity[postmelt_viscosity < liquid_viscosity]
    postmelt_shear_modulus[postmelt_shear_modulus < liquid_shear] = liquid_shear

    return postmelt_viscosity, postmelt_shear_modulus

# Put New Models Below Here!
