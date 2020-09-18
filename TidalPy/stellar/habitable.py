import numpy as np
from scipy.constants import Stefan_Boltzmann as sbc

from ..utilities.performance.numba import njit
from ..utilities.types import FloatArray

INNER_EDGE_TEMP = 273.15
OUTER_EDGE_TEMP = 373.15
INNER_EDGE_TEMP_4 = INNER_EDGE_TEMP**4
OUTER_EDGE_TEMP_4 = OUTER_EDGE_TEMP**4


@njit(cacheable=True)
def calc_equilibrium_temperature(insolation_heating: FloatArray, radius: float, internal_heating: FloatArray = None,
                                 emissivity: float = 1., internal_to_surf_frac: float = 1.0):
    """ Calculates the surface equilibrium temperature of a planet that is heated by stellar and internal heating

    References
    ----------
    .. Henning et al. (2020) Exoplanet Habitable Zone Modifications

    Parameters
    ----------
    insolation_heating : FloatArray
        Heat received from a host star [W]
    radius : float
        Radius of the world [m]
    internal_heating : FloatArray = None
        Heat received from the interior of the world [W]
    emissivity : float = 1.0
        Planet's grey-body emissivity
    internal_to_surf_frac: float = 1.0
        Fraction of the internal heat that reaches the surface (scales the internal heating rate)

    Returns
    -------
    equilibrium_temperature : FloatArray
        World's surface equilibrium temperature [K]

    """

    heating = insolation_heating
    if internal_heating is not None:
        heating += internal_heating * internal_to_surf_frac

    coeff = 4. * np.pi * radius**2 * emissivity * sbc / 2.
    equilibrium_temperature = (heating / coeff)**(1 / 4)

    return equilibrium_temperature
