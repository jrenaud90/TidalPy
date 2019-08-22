""" Partial Melting Models for Viscosity and Shear Modulus """

from typing import Tuple

import numpy as np

from ..performance import njit


@njit
def off(temperature: np.ndarray, melt_fraction: np.ndarray,
        premelt_viscosity: np.ndarray, premelt_shear: np.ndarray, liquid_viscosity: np.ndarray):
    """ No Partial Melt

    """

    return premelt_viscosity, premelt_shear


@njit
def spohn(temperature: np.ndarray, melt_fraction: np.ndarray,
          premelt_viscosity: np.ndarray, premelt_shear: np.ndarray, liquid_viscosity: np.ndarray,
          liquid_shear: float, fs_visc_power_slope: float, fs_visc_power_phase: float,
          fs_shear_power_slope: float, fs_shear_power_phase: float) -> Tuple[np.ndarray, np.ndarray]:
    """ Fischer and Spohn (1990) Partial-Melt Viscosity Function

    !TPY_args const: liquid_shear, fs_visc_power_slope, fs_visc_power_phase, fs_shear_power_slope, fs_shear_power_phase

    """

    viscosity = 10.**((fs_visc_power_slope / temperature) - fs_visc_power_phase)
    shear = 10.**((fs_shear_power_slope / temperature) - fs_shear_power_phase)

    # Perform sanity checks
    viscosity[viscosity < liquid_viscosity] = liquid_viscosity[viscosity < liquid_viscosity]
    shear[shear < liquid_shear] = liquid_shear

    return viscosity, shear


@njit
def henning(temperature: np.ndarray, melt_fraction: np.ndarray,
            premelt_viscosity: np.ndarray, premelt_shear: np.ndarray, liquid_viscosity: np.ndarray,
            liquid_shear: float, crit_melt_frac: float, crit_melt_frac_width: float, hn_visc_slope_1: float,
            hn_visc_slope_2: float, hn_shear_param_1: float, hn_shear_param_2: float,
            hn_shear_falloff_slope: float) -> Tuple[np.ndarray, np.ndarray]:
    """ Partial-Melting Viscosity and Shear Function based on Henning+ (2009)

    See also Renaud and Henning (2009)

    !TPY_args const: liquid_shear, crit_melt_frac, crit_melt_frac_width, hn_visc_slope_1, hn_visc_slope_2, hn_shear_param_1, hn_shear_param_2, hn_shear_falloff_slope

    """

    crit_melt_frac_plus_width = crit_melt_frac + crit_melt_frac_width

    # Calculate ndarray indices
    pre_breakdown_index = np.logical_and(crit_melt_frac > melt_fraction, melt_fraction > 0)
    breakdown_index = np.logical_and(melt_fraction >= crit_melt_frac, crit_melt_frac_plus_width >= melt_fraction)
    molten_index = melt_fraction > crit_melt_frac_plus_width

    # Calculate viscosity in the three domains
    # TODO: Should this make a copy of this array before it starts modifying its values?
    viscosity = premelt_viscosity
    viscosity[pre_breakdown_index] *= np.exp(-hn_visc_slope_1 * melt_fraction[pre_breakdown_index])
    viscosity[breakdown_index] *= np.exp(-hn_visc_slope_2 * (melt_fraction[breakdown_index] - crit_melt_frac))
    viscosity[molten_index] = liquid_viscosity[molten_index]
    shear = premelt_shear * np.ones_like(temperature)
    shear[pre_breakdown_index] *= np.exp((hn_shear_param_1 / temperature[pre_breakdown_index]) - hn_shear_param_2)
    shear[breakdown_index] *= np.exp(-hn_shear_falloff_slope * (melt_fraction[breakdown_index] - crit_melt_frac))
    shear[molten_index] = liquid_shear

    # Perform sanity checks
    viscosity[viscosity < liquid_viscosity] = liquid_viscosity[viscosity < liquid_viscosity]
    shear[shear < liquid_shear] = liquid_shear

    return viscosity, shear

# Put New Models Below Here!
