from numba import njit
import numpy as np
from typing import NamedTuple


class Strength(NamedTuple):

    viscosity: np.ndarray
    shear_modulus: np.ndarray


@njit
def melt_fraction(temperature, solidus, liquidus):
    """ Calculates the volumetric melt fraction

    :param temperature: <Array> Material/Layer temperature [K]
    :param solidus: <Float> Material/Layer solidus temperature [K]
    :param liquidus: <Float> Material/Layer liquidus temperature [K]
    :return: <Array>
    """

    melt_frac = (temperature - solidus) / (liquidus - solidus)

    # Check for over/under-shoots
    melt_frac[melt_frac < 0.] = 0.
    melt_frac[melt_frac > 1.] = 1.

    return melt_frac