""" Partial Melting Models for Viscosity and Shear Modulus """

from typing import Tuple

import numpy as np

from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray


@njit(cacheable=True)
def off(melt_fraction: 'FloatArray', premelt_viscosity: 'FloatArray', premelt_shear: FloatArray) -> \
        Tuple[FloatArray, FloatArray]:
    """ Viscosity and Shear Modulus Partial Melting Model: off

    !TPY_args live: self.premelt_viscosity, self.premelt_shear

    Parameters
    ----------
    melt_fraction : FloatArray
        Layer/Material volumetric melt fraction [m3 m-3]
    premelt_viscosity : FloatArray
        Layer/Material viscosity before partial melting is considered [Pa s]
    premelt_shear : FloatArray
        Layer/Material shear modulus before partial melting is considered [Pa]

    Returns
    -------
    postmelt_viscosity : FloatArray
        Post melting viscosity
    postmelt_shear_modulus : FloatArray
        Post melting shear modulus
    """

    shape = (melt_fraction + premelt_viscosity + premelt_shear) * 0.
    postmelt_viscosity = premelt_viscosity + shape
    postmelt_shear_modulus = premelt_shear + shape

    return postmelt_viscosity, postmelt_shear_modulus

@njit(cacheable=True)
def spohn(
    melt_fraction: 'FloatArray', temperature: 'FloatArray', liquid_viscosity: 'FloatArray',
    liquid_shear: float = 1.0e-5,
    fs_visc_power_slope: float = 27000.0, fs_visc_power_phase: float = 1.0,
    fs_shear_power_slope: float = 82000.0, fs_shear_power_phase: float = 40.6
    ) -> Tuple[FloatArray, FloatArray]:
    """ Viscosity and Shear Modulus Partial Melting Model: spohn

    Fischer and Spohn (1990) Partial-Melt Viscosity Function

    !TPY_args live: self.temperature, self.liquid_viscosity
    !TPY_args const: liquid_shear, fs_visc_power_slope, fs_visc_power_phase, fs_shear_power_slope, fs_shear_power_phase

    Parameters
    ----------
    melt_fraction : FloatArray
        Layer/Material volumetric melt fraction [m3 m-3]
    temperature : FloatArray
        Layer/Material temperature [K]
    liquid_viscosity : FloatArray
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
    postmelt_viscosity : FloatArray
        Post melting viscosity
    postmelt_shear_modulus : FloatArray
        Post melting shear modulus
    """

    shape = (melt_fraction + temperature + liquid_viscosity) * 0.
    postmelt_viscosity = 10.**((fs_visc_power_slope / temperature) - fs_visc_power_phase)
    postmelt_shear_modulus = 10.**((fs_shear_power_slope / temperature) - fs_shear_power_phase)

    # Perform sanity checks
    postmelt_viscosity = (postmelt_viscosity + shape <= liquid_viscosity) * liquid_viscosity + \
                         (postmelt_viscosity + shape > liquid_viscosity) * postmelt_viscosity
    postmelt_shear_modulus = (postmelt_shear_modulus + shape <= liquid_shear) * liquid_shear + \
                             (postmelt_shear_modulus + shape > liquid_shear) * postmelt_shear_modulus

    return postmelt_viscosity, postmelt_shear_modulus

@njit(cacheable=True)
def henning(
    melt_fraction: 'FloatArray', temperature: 'FloatArray',
    premelt_viscosity: 'FloatArray', liquid_viscosity: 'FloatArray', premelt_shear: 'FloatArray',
    solidus: float, liquidus: float,
    liquid_shear: float, crit_melt_frac: float = 0.5, crit_melt_frac_width: float = 0.05,
    hn_visc_slope_1: float = 13.5, hn_visc_falloff_slope: float = 370., hn_shear_param_1: float = 40000.,
    hn_shear_param_2: float = 25., hn_shear_falloff_slope: float = 700.
    ) -> Tuple[FloatArray, FloatArray]:
    """ Viscosity and Shear Modulus Partial Melting Model: henning

    Henning (2009, 2010) Partial-Melt Viscosity Function
    See also Moore's work and Renaud and Henning (2018)

    !TPY_args live: self.temperature, self.premelt_viscosity, self.liquid_viscosity, self.premelt_shear, self.solidus, self.liquidus
    !TPY_args const: liquid_shear, crit_melt_frac, crit_melt_frac_width, hn_visc_slope_1, hn_visc_slope_2, hn_shear_param_1, hn_shear_param_2, hn_shear_falloff_slope

    Parameters
    ----------
    melt_fraction : FloatArray
        Layer/Material volumetric melt fraction [m3 m-3]
    temperature : FloatArray
        Layer/Material temperature [K]
    premelt_viscosity : FloatArray
        Layer/Material viscosity before partial melting is considered [Pa s]
    premelt_shear : FloatArray
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
    postmelt_viscosity : FloatArray
        Post melting viscosity
    postmelt_shear_modulus : FloatArray
        Post melting shear modulus
    """

    shape = (melt_fraction + temperature + premelt_viscosity + liquid_viscosity + premelt_shear) * 0.
    melt_fraction_shape = melt_fraction + shape

    crit_melt_frac_plus_width = crit_melt_frac + crit_melt_frac_width
    break_down_temp = solidus + crit_melt_frac * (liquidus - solidus)

    # Initialize post-melt with pre-melt values (if there is no melt fraction then these values will be returned.
    premelt_viscosity = premelt_viscosity
    premelt_shear_modulus = premelt_shear

    # Calculate viscosity and shear modulus in the three domains
    # # Partial melting
    # if melt_fraction < crit_melt_frac:
    #     # Partial melting before critical break-down.
    # elif crit_melt_frac <= melt_fraction <= crit_melt_frac_plus_width:
    #     # Get the maximum pre-melt effect
    #     # Then apply critical breakdown occurring.
    # elif melt_fraction > crit_melt_frac_plus_width:
    #     # Past critical breakdown threshold. The material now behaves like a fluid.

    postmelt_viscosity = \
        (melt_fraction_shape <= 0.) * \
            premelt_viscosity + \
        (melt_fraction_shape > 0.) * (melt_fraction_shape < crit_melt_frac) * \
            premelt_viscosity * np.exp(-hn_visc_slope_1 * melt_fraction) + \
        (melt_fraction_shape >= crit_melt_frac) * (melt_fraction_shape <= crit_melt_frac_plus_width) * \
            premelt_viscosity * np.exp(-hn_visc_slope_1 * crit_melt_frac) * \
            np.exp(-hn_visc_falloff_slope * (melt_fraction - crit_melt_frac)) + \
        (melt_fraction_shape > crit_melt_frac_plus_width) * \
            liquid_viscosity

    postmelt_shear_modulus = \
        (melt_fraction_shape <= 0.) * \
            premelt_shear_modulus + \
        (melt_fraction_shape > 0.) * (melt_fraction_shape < crit_melt_frac) * \
            premelt_shear_modulus * np.exp((hn_shear_param_1 / temperature) - hn_shear_param_2) + \
        (melt_fraction_shape >= crit_melt_frac) * (melt_fraction_shape <= crit_melt_frac_plus_width) * \
            premelt_shear_modulus * np.exp((hn_shear_param_1 / break_down_temp) - hn_shear_param_2) * \
            np.exp(-hn_shear_falloff_slope * (melt_fraction - crit_melt_frac)) + \
        (melt_fraction_shape > crit_melt_frac_plus_width) * \
            liquid_viscosity

    # Perform sanity checks
    postmelt_viscosity = (postmelt_viscosity + shape <= liquid_viscosity) * liquid_viscosity + \
                         (postmelt_viscosity + shape > liquid_viscosity) * postmelt_viscosity
    postmelt_shear_modulus = (postmelt_shear_modulus + shape <= liquid_shear) * liquid_shear + \
                             (postmelt_shear_modulus + shape > liquid_shear) * postmelt_shear_modulus

    return postmelt_viscosity, postmelt_shear_modulus

# Put New Models Below Here!
