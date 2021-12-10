""" Tests for extracting useful information out of a multilayer tidal propagation
"""

import numpy as np

import TidalPy
from TidalPy.constants import G
from TidalPy.tides.potential import tidal_potential_simple, tidal_potential_nsr, tidal_potential_obliquity_nsr
from TidalPy.toolbox.conversions import orbital_motion2semi_a

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex)
host_mass = 10. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.01


def test_tidal_potential_simple():
    """ Test the basic tidal potential assuming low eccentricity, no obliquity, and synchronous rotation """
    # Test floats
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_simple(
            radius_array[-1], longitude=0.1, colatitude=0.1, orbital_frequency=orbital_freq,
            eccentricity=eccentricity, time=1000.
            )

    assert type(potential) in [float, np.float, np.float64]
    assert type(potential_partial_theta) in [float, np.float, np.float64]
    assert type(potential_partial_phi) in [float, np.float, np.float64]
    assert type(potential_partial2_theta2) in [float, np.float, np.float64]
    assert type(potential_partial2_phi2) in [float, np.float, np.float64]
    assert type(potential_partial2_theta_phi) in [float, np.float, np.float64]

    # Test arrays
    arr = np.linspace(0.1, 0.8, 10)
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_simple(
            radius=arr, longitude=arr,
            colatitude=arr, orbital_frequency=orbital_freq,
            eccentricity=eccentricity, time=1000.
            )

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == arr.shape
    assert potential_partial_theta.shape == arr.shape
    assert potential_partial_phi.shape == arr.shape
    assert potential_partial2_theta2.shape == arr.shape
    assert potential_partial2_phi2.shape == arr.shape
    assert potential_partial2_theta_phi.shape == arr.shape

    assert type(potential[0]) in [float, np.float, np.float64]
    assert type(potential_partial_theta[0]) in [float, np.float, np.float64]
    assert type(potential_partial_phi[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_theta2[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_phi2[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_theta_phi[0]) in [float, np.float, np.float64]

def test_tidal_potential_obliquity_nsr():
    """ Test the tidal potential equation assuming low eccentricity, no obliquity, and non-synchronous rotation """

    # Test floats
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_obliquity_nsr(
            radius_array[-1] / 2, longitude=0.1, colatitude=0.1, orbital_frequency=orbital_freq,
            eccentricity=eccentricity, time=1000., obliquity=0.5, rotation_rate=1.5 * orbital_freq, periapsis=0.,
            world_radius=radius_array[-1]
            )

    assert type(potential) in [float, np.float, np.float64]
    assert type(potential_partial_theta) in [float, np.float, np.float64]
    assert type(potential_partial_phi) in [float, np.float, np.float64]
    assert type(potential_partial2_theta2) in [float, np.float, np.float64]
    assert type(potential_partial2_phi2) in [float, np.float, np.float64]
    assert type(potential_partial2_theta_phi) in [float, np.float, np.float64]

    # Test arrays
    arr = np.linspace(0.1, 0.8, 10)
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_obliquity_nsr(
            radius=arr, longitude=arr, colatitude=arr, orbital_frequency=orbital_freq,
            eccentricity=eccentricity, time=1000., obliquity=0.5, rotation_rate=1.5 * orbital_freq, periapsis=0.,
            world_radius=radius_array[-1]
            )

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == arr.shape
    assert potential_partial_theta.shape == arr.shape
    assert potential_partial_phi.shape == arr.shape
    assert potential_partial2_theta2.shape == arr.shape
    assert potential_partial2_phi2.shape == arr.shape
    assert potential_partial2_theta_phi.shape == arr.shape

    assert type(potential[0]) in [float, np.float, np.float64]
    assert type(potential_partial_theta[0]) in [float, np.float, np.float64]
    assert type(potential_partial_phi[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_theta2[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_phi2[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_theta_phi[0]) in [float, np.float, np.float64]

def test_tidal_potential_nsr():
    """ Test the tidal potential equation assuming moderate eccentricity, no obliquity, and non-synchronous rotation """

    # Test floats
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_nsr(
            radius_array[-1] / 2, longitude=0.1, colatitude=0.1, orbital_frequency=orbital_freq,
            eccentricity=eccentricity, time=1000., rotation_rate=1.5 * orbital_freq,
            world_radius=radius_array[-1], host_mass=1.0e9, semi_major_axis=1.0e12, use_static=False
            )

    assert type(potential) in [float, np.float, np.float64]
    assert type(potential_partial_theta) in [float, np.float, np.float64]
    assert type(potential_partial_phi) in [float, np.float, np.float64]
    assert type(potential_partial2_theta2) in [float, np.float, np.float64]
    assert type(potential_partial2_phi2) in [float, np.float, np.float64]
    assert type(potential_partial2_theta_phi) in [float, np.float, np.float64]

    # Test arrays
    arr = np.linspace(0.1, 0.8, 10)
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_nsr(
            radius=arr, longitude=arr, colatitude=arr, orbital_frequency=orbital_freq,
            eccentricity=eccentricity, time=1000., rotation_rate=1.5 * orbital_freq,
            world_radius=radius_array[-1], host_mass=1.0e9, semi_major_axis=1.0e12, use_static=False
            )

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == arr.shape
    assert potential_partial_theta.shape == arr.shape
    assert potential_partial_phi.shape == arr.shape
    assert potential_partial2_theta2.shape == arr.shape
    assert potential_partial2_phi2.shape == arr.shape
    assert potential_partial2_theta_phi.shape == arr.shape

    assert type(potential[0]) in [float, np.float, np.float64]
    assert type(potential_partial_theta[0]) in [float, np.float, np.float64]
    assert type(potential_partial_phi[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_theta2[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_phi2[0]) in [float, np.float, np.float64]
    assert type(potential_partial2_theta_phi[0]) in [float, np.float, np.float64]