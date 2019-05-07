from numba import njit
import numpy as np
from ..utilities.search import ModelSearcher
from . import viscosity

# Helper Functions
@njit
def melt_fraction(temperature: np.ndarray, solidus: float, liquidus: float) -> np.ndarray:
    """ Calculates the volumetric melt fraction

    :param temperature: <Array> Material/Layer temperature [K]
    :param solidus: <Float> Material/Layer solidus temperature [K]
    :param liquidus: <Float> Material/Layer liquidus temperature [K]
    :return: <Array> Volumetric Melt Fraction [m3 m-3]
    """

    melt_frac = (temperature - solidus) / (liquidus - solidus)

    # Check for over/under-shoots
    melt_frac[melt_frac < 0.] = 0.
    melt_frac[melt_frac > 1.] = 1.

    return melt_frac


@njit
def temp_from_melt(melt_frac: np.ndarray, solidus: float, liquidus: float) -> np.ndarray:
    """ Calculates the temperature from the volumetric melt fraction

    :param melt_frac: <Array> Volumetric Melt Fraction [m3 m-3]
    :param solidus: <Float> Material/Layer solidus temperature [K]
    :param liquidus: <Float> Material/Layer liquidus temperature [K]
    :return: <Array> Temperature [K]
    """

    # Check for over/under-shoots
    melt_frac[melt_frac < 0.] = 0.
    melt_frac[melt_frac > 1.] = 1.

    return melt_frac * (liquidus - solidus) + solidus


# Partial Melt Model Finder
partial_melter_param_defaults = {
    'ice': {

        },

    'rocky': {
            'fs_visc_power_slope': TBD,
            'fs_visc_power_phase': TBD,
            'liquid_shear': TBD,
            'crit_melt_frac': TBD,
            'crit_melt_frac_width': TBD,
            'hn_visc_slope_1': TBD,
            'hn_visc_slope_2': TBD,
            'hn_shear_param_1': TBD,
            'hn_shear_param_2': TBD,
            'hn_shear_falloff_slope': TBD
        },

    'iron': {

        }
    }

find_partial_melter = ModelSearcher(viscosity, partial_melter_param_defaults)
