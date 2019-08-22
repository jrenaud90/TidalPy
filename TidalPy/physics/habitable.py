import numpy as np
from scipy.constants import Stefan_Boltzmann as sbc

from ..performance import njit
from ..types import FloatArray


INNER_EDGE_TEMP = 273.15
OUTER_EDGE_TEMP = 373.15
INNER_EDGE_TEMP_4 = INNER_EDGE_TEMP**4
OUTER_EDGE_TEMP_4 = OUTER_EDGE_TEMP**4


@njit
def equilibrium_temperature(insolation_heating: FloatArray, radius: float, internal_heating: FloatArray = None,
                            emissivity: float = 1., internal_frac: float = .5):
    """ Calculates the surface equilibrium temperature of a planet that is heated by stellar and internal heating

    Based on Henning habitability preprint

    :param insolation_heating: <FloatArray> Heat received from a host star [Watts]
    :param radius:             <float> Planet's radius [m]
    :param internal_heating:   <FloatArray> Heat leaving the interior [Watts]
    :param emissivity:         <float> Planet's greybody emissivity
    :param internal_frac:      <float> Fraction of the internal heat making it to the surface
    :return:                   <float> Planet's surface equilibrium temperature [K]
    """

    heating = insolation_heating
    if insolation_heating is not None:
        heating += insolation_heating

    coeff = 4. * np.pi * radius**2 * emissivity * sbc / 2

    return (heating / coeff)**(1 / 4)
