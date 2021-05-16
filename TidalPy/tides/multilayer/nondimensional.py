""" The functions here are used to non-dimensionalize common variables used for multi-layer tidal calculations.

The scheme is based on that proposed by Martens16 (around page 99)

References
----------
Martens16 : H. Martens, PhD Thesis (CalTech), 2016, DOI: 10.7907/Z9N29TX7

"""
from typing import Tuple

import numpy as np

from ...constants import G
from ...utilities.performance import njit
from ...utilities.types import FloatArray, NumArray


NonDimPhysicalOutput = Tuple[FloatArray, FloatArray, FloatArray, NumArray, NumArray, FloatArray, float]
ReDimPhysicalOutput = Tuple[FloatArray, FloatArray, FloatArray, NumArray, NumArray, FloatArray]


@njit(cacheable=True)
def non_dimensionalize_physicals(radius: FloatArray, gravity: FloatArray, density: FloatArray,
                                 shear_modulus: NumArray, bulk_modulus: NumArray, frequency: FloatArray,
                                 mean_radius: float, bulk_density: float) -> NonDimPhysicalOutput:
    """ A function to non-dimensionalize physical parameters

    Parameters
    ----------
    radius : FloatArray
        Radius of the planet (usually an array) [m]
    gravity : FloatArray
        Acceleration due to gravity of the planet (usually an array) [m s-2]
    density : FloatArray
        Density of the planet (usually an array) [kg m-3]
    shear_modulus : NumArray
        Shear modulus of the planet (usually an array, can be complex) [Pa]
    bulk_modulus : NumArray
        Bulk modulus of the planet (usually an array, can be complex) [Pa]
    frequency : FloatArray
        Forcing frequency of the tidal forces [rad s-1]
    mean_radius : float
        Mean radius of the planet, used in scaling [m]
    bulk_density : float
        Bulk density of the planet, used in scaling [m]

    Returns
    -------
    radius_prime : FloatArray
        Non-dimensional radius of the planet
    gravity_prime : FloatArray
        Non-dimensional Acceleration due to gravity of the planet
    density_prime : FloatArray
        Non-dimensional Density of the planet
    shear_modulus_prime : NumArray
        Non-dimensional Shear modulus of the planet
    bulk_modulus_prime : NumArray
        Non-dimensional Bulk modulus of the planet
    frequency_prime : FloatArray
        Non-dimensional forcing frequency of the tidal forces
    newton_g_prime : float
        Non-dimensional Newton's gravitational constant

    See Also
    --------
    nondimensional.re_dimensionalize_physicals

    """

    # Setup conversions
    second2_conversion = 1. / (np.pi * G * bulk_density)
    second_conversion = np.sqrt(second2_conversion)
    mass_conversion = bulk_density * mean_radius**3
    length_conversion = mean_radius
    density_conversion = bulk_density
    pascal_conversion = mass_conversion / (length_conversion * second2_conversion)

    # Convert variables
    radius_prime = radius / length_conversion
    density_prime = density / density_conversion
    gravity_prime = gravity / (length_conversion / second2_conversion)
    shear_modulus_prime = shear_modulus / pascal_conversion
    bulk_modulus_prime = bulk_modulus / pascal_conversion
    frequency_prime = frequency / (1. / second_conversion)

    # Calculations will also need the universal gravitational constant to be non-dim
    newton_g_prime = 1. / np.pi

    return radius_prime, gravity_prime, density_prime, shear_modulus_prime, bulk_modulus_prime, frequency_prime, \
           newton_g_prime


@njit(cacheable=True)
def re_dimensionalize_physicals(radius_prime: FloatArray, gravity_prime: FloatArray, density_prime: FloatArray,
                                shear_modulus_prime: NumArray, bulk_modulus_prime: NumArray,
                                frequency_prime: FloatArray,
                                mean_radius: float, bulk_density: float) -> ReDimPhysicalOutput:
    """ A function to re-dimensionalize physical parameters that have been previously non-dimensionalized.

    Parameters
    ----------
    radius_prime : FloatArray
        Non-dimensional radius of the planet
    gravity_prime : FloatArray
        Non-dimensional Acceleration due to gravity of the planet
    density_prime : FloatArray
        Non-dimensional Density of the planet
    shear_modulus_prime : NumArray
        Non-dimensional Shear modulus of the planet
    bulk_modulus_prime : NumArray
        Non-dimensional Bulk modulus of the planet
    frequency_prime : FloatArray
        Non-dimensional forcing frequency of the tidal forces
    mean_radius : float
        Mean radius of the planet, used in scaling [m]
    bulk_density : float
        Bulk density of the planet, used in scaling [m]

    Returns
    -------
    radius : FloatArray
        Radius of the planet (usually an array) [m]
    gravity : FloatArray
        Acceleration due to gravity of the planet (usually an array) [m s-2]
    density : FloatArray
        Density of the planet (usually an array) [kg m-3]
    shear_modulus : NumArray
        Shear modulus of the planet (usually an array, can be complex) [Pa]
    bulk_modulus : NumArray
        Bulk modulus of the planet (usually an array, can be complex) [Pa]
    frequency : FloatArray
        Forcing frequency of the tidal forces [rad s-1]

    See Also
    --------
    nondimensional.non_dimensionalize_physicals

    """

    # Setup conversions
    second2_conversion = 1. / (np.pi * G * bulk_density)
    second_conversion = np.sqrt(second2_conversion)
    mass_conversion = bulk_density * mean_radius**3
    length_conversion = mean_radius
    density_conversion = bulk_density
    pascal_conversion = mass_conversion / (length_conversion * second2_conversion)

    # Convert variables
    radius = radius_prime * length_conversion
    density = density_prime * density_conversion
    gravity = gravity_prime / (second2_conversion / length_conversion)
    shear_modulus = shear_modulus_prime * pascal_conversion
    bulk_modulus = bulk_modulus_prime * pascal_conversion
    frequency = frequency_prime / second_conversion

    return radius, gravity, density, shear_modulus, bulk_modulus, frequency


@njit(cacheable=True)
def re_dimensionalize_radial_func(tidal_y_prime: np.ndarray,
                                  mean_radius: float, bulk_density: float) -> np.ndarray:
    """ A function to re-dimensionalize physical parameters that have been previously non-dimensionalized.

    Parameters
    ----------
    tidal_y_prime : np.ndarray
        Non-dimensionalized radial solutions as a function of radius.
    mean_radius : float
        Mean radius of the planet, used in scaling [m]
    bulk_density : float
        Bulk density of the planet, used in scaling [m]

    Returns
    -------
    tidal_y : np.ndarray
        Re-dimensionalized radial solutions as a function of radius.

    """

    # Setup conversions
    second2_conversion = 1. / (np.pi * G * bulk_density)
    mass_conversion = bulk_density * mean_radius**3
    length_conversion = mean_radius

    # Convert displacements
    #    y1, y3 are the radial and tangential displacements with units of [s2 m-1]
    #    y2, y4 are the radial and tangential stresses with units of [kg m-3]
    #    y5 is the tidal potential which is unitl ess and thus needs no conversion.
    #    y6 is a "potential stress" with units of [m-1]
    tidal_y = np.zeros_like(tidal_y_prime)
    tidal_y[0] = tidal_y_prime[0] * (second2_conversion / length_conversion)
    tidal_y[2] = tidal_y_prime[2] * (second2_conversion / length_conversion)

    tidal_y[1] = tidal_y_prime[1] * (mass_conversion / length_conversion**3)
    tidal_y[3] = tidal_y_prime[3] * (mass_conversion / length_conversion**3)

    tidal_y[4] = tidal_y_prime[4]
    tidal_y[5] = tidal_y_prime[5] * (1. / length_conversion)

    return tidal_y
