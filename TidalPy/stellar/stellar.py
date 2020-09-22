import numpy as np
from scipy.constants import Stefan_Boltzmann as sbc

from ..constants import luminosity_solar, mass_solar
from ..utilities.performance.numba import njit

@njit(cacheable=True)
def efftemp_from_luminosity(luminosity: float, radius: float):
    """ Calculates a star's effective surface temperature provided a luminosity

    :param luminosity: <float> StarWorld's luminosity [Watts]
    :param radius:     <float> StarWorld's Surface Radius [m]
    :return:           <float> StarWorld's Effective Surface temperature [K]
    """

    return (luminosity / (4. * np.pi * sbc * radius**2))**(1 / 4)


@njit(cacheable=True)
def luminosity_from_efftemp(effective_temperature: float, radius: float):
    """ Calculates a star's luminosity provided an effective surface temperature

    :param effective_temperature: <Float> StarWorld's Effective Surface temperature [K]
    :param radius:     <float> StarWorld's Surface Radius [m]
    :return:           <float> StarWorld's luminosity [Watt]
    """

    return 4. * np.pi * sbc * radius**2 * effective_temperature**4


# Stellar Relationships and Scaling
@njit(cacheable=True)
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
