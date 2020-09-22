import numpy as np
from scipy.special import ellipe
from scipy.constants import Stefan_Boltzmann as sbc

from ..utilities.performance.numba import njit
from ..utilities.types import FloatArray


# Equilibrium Temperature
# TODO: njit does not like the heating addition checking for None.
# @njit(cacheable=True)
def calc_equilibrium_temperature(insolation_heating: FloatArray, radius: float, internal_heating: FloatArray = None,
                                 emissivity: float = 1.):
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

    Returns
    -------
    equilibrium_temperature : FloatArray
        World's surface equilibrium temperature [K]

    """

    heating = insolation_heating
    if internal_heating is not None:
        heating += internal_heating

    # TODO: There was a divisor of 1/2 on this coeff before. I don't recognize it from anywhere. Removed for now.
    # coeff = 4. * np.pi * radius**2 * emissivity * sbc / 2.
    coeff = 4. * np.pi * radius**2 * emissivity * sbc

    equilibrium_temperature = (heating / coeff)**(1 / 4)
    return equilibrium_temperature


# Insolation Heating Functions
# TODO: the ellipe function can't be wrapped by njit. Perhaps one day it will or we could make our own njit-friendly wrapping
# @njit
def equilibrium_insolation_mendez(luminosity: float, semi_major_axis: FloatArray, albedo: float,
                                  radius: float, eccentricity: FloatArray) -> FloatArray:
    """ Calculate the insolation heating using the mendez method for elliptical orbits

    Based on Mendez & Rivera-Valentin (ApJL 837 1, 2017)

    Parameters
    ----------
    luminosity : float
        Stellar luminosity in [W]
    semi_major_axis : FloatArray
        Target planet's semi-major axis to the star (you want to use the host planet's semi-a) in [m]
    albedo : float
        Target planet's geometric albedo
    radius : float
        Target planet's radius in [m]
    eccentricity : FloatArray
        Target planet's semi-major axis relative to the star (you want to use the host planet's e)

    Returns
    -------
    insolation_heating : FloatArray
        Heating at the target planet's surface in [W]

    See Also
    --------
    equilibrium_insolation_williams
    equilibrium_insolation_no_eccentricity
    """

    mendez_correction = (2 * np.sqrt(1 + eccentricity) / np.pi) * \
                        ellipe(np.sqrt(2 * eccentricity / (1. + eccentricity)))
    insolation_heating = mendez_correction * radius**2 * (1. - albedo) * luminosity / (4. * semi_major_axis**2)

    return insolation_heating


@njit(cacheable=True)
def equilibrium_insolation_williams(luminosity: float, semi_major_axis: FloatArray, albedo: float,
                                    radius: float, eccentricity: FloatArray) -> FloatArray:
    """ Calculate the insolation heating using the williams method for elliptical orbits

    Based on *NEED WILLIAMS REF

    Parameters
    ----------
    luminosity : float
        Stellar luminosity in [W]
    semi_major_axis : FloatArray
        Target planet's semi-major axis to the star (you want to use the host planet's semi-a) in [m]
    albedo : float
        Target planet's geometric albedo
    radius : float
        Target planet's radius in [m]
    eccentricity : FloatArray
        Target planet's semi-major axis relative to the star (you want to use the host planet's e)

    Returns
    -------
    insolation_heating : FloatArray
        Heating at the target planet's surface in [W]

    See Also
    --------
    equilibrium_insolation_no_eccentricity
    equilibrium_insolation_mendez
    """

    insolation_heating = radius**2 * (1. - albedo) * luminosity / (4. * semi_major_axis**2 * np.sqrt(1. - eccentricity))

    return insolation_heating


@njit(cacheable=True)
def equilibrium_insolation_no_eccentricity(luminosity: float, semi_major_axis: FloatArray, albedo: float,
                                           radius: float, eccentricity: FloatArray = None) -> FloatArray:
    """ Calculate the insolation heating assuming a circular orbit

    Parameters
    ----------
    luminosity : float
        Stellar luminosity in [W]
    semi_major_axis : FloatArray
        Target planet's semi-major axis to the star (you want to use the host planet's semi-a) in [m]
    albedo : float
        Target planet's geometric albedo
    radius : float
        Target planet's radius in [m]
    eccentricity : FloatArray
        This is not used - It is here to keep the signature of all the insolation heating functions the same

    Returns
    -------
    insolation_heating : FloatArray
        Heating at the target planet's surface in [W]

    See Also
    --------
    equilibrium_insolation_williams
    equilibrium_insolation_mendez
    """

    insolation_heating = radius**2 * (1. - albedo) * luminosity / (4. * semi_major_axis**2)

    return insolation_heating