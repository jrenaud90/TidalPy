import numpy as np
from scipy.constants import Stefan_Boltzmann as sbc
from scipy.special import ellipe

from ..constants import luminosity_solar, mass_solar
from ..utilities.performance.numba import njit
from ..utilities.types import FloatArray


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


@njit
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


@njit
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


@njit
def equilibrium_temperature(surface_heating: FloatArray, radius: float, emissivity: float):
    """ Calculate the surface equilibrium temperature for a provided surface heating

    Parameters
    ----------
    surface_heating : FloatArray
        Summation of all surface heating in [W]
    radius : float
        Target planet's surface radius (or mean radius) in [m]
    emissivity : float
        Target planet's greybody emissivity (=1 for perfect blackbody)

    Returns
    -------
    surf_equilibrium_temperature : FloatArray
        Surface equilibrium temperature in [K]
    """

    surf_equilibrium_temperature = (surface_heating / (4. * np.pi * radius**2 * sbc * emissivity))**(1 / 4)
    return surf_equilibrium_temperature


@njit
def efftemp_from_luminosity(luminosity: float, radius: float):
    """ Calculates a star's effective surface temperature provided a luminosity

    :param luminosity: <float> StarWorld's luminosity [Watts]
    :param radius:     <float> StarWorld's Surface Radius [m]
    :return:           <float> StarWorld's Effective Surface temperature [K]
    """

    return (luminosity / (4. * np.pi * sbc * radius**2))**(1 / 4)


@njit
def luminosity_from_efftemp(effective_temperature: float, radius: float):
    """ Calculates a star's luminosity provided an effective surface temperature

    :param effective_temperature: <Float> StarWorld's Effective Surface temperature [K]
    :param radius:     <float> StarWorld's Surface Radius [m]
    :return:           <float> StarWorld's luminosity [Watt]
    """

    return 4. * np.pi * sbc * radius**2 * effective_temperature**4


# Relationships and Scaling
@njit
def luminosity_from_mass(stellar_mass: float):
    """ Estimates stellar luminosity from a star's mass

    Partially based on Cuntz & Wang 2018 (doi:10.3847/2515-5172/aaaa67) and wikipedia.org/wiki/Mass–luminosity_relation

    :param stellar_mass: <float> StarWorld's mass [kg]
    :return:             <float> StarWorld's luminosity [Watts]
    """

    mass_ratio = stellar_mass / mass_solar

    if mass_ratio < .2:
        return luminosity_solar * 0.23 * mass_ratio**2.3
    if mass_ratio < 0.85:
        a = -141.7 * stellar_mass**4 + 232.4 * stellar_mass**3 - 129.1 * stellar_mass**2 + 33.29 * stellar_mass + 0.215
        return luminosity_solar * mass_ratio**a
    if mass_ratio < 2.:
        return luminosity_solar * mass_ratio**4
    if mass_ratio < 20.:
        return luminosity_solar * 1.4 * mass_ratio**3.5
    return luminosity_solar * 3.2e4 * mass_ratio
