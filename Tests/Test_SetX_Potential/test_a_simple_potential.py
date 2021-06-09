""" Tests for extracting useful information out of a multilayer tidal propagation
"""

import numpy as np

import TidalPy
from TidalPy.tides.potential import tidal_potential_simple
from TidalPy.constants import G
from TidalPy.toolbox.conversions import orbital_motion2semi_a

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array =  np.linspace(0., 1.e6, 11)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i+1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex)
host_mass = 10. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.01

def test_float():
    potential, potential_partial_theta, potential_partial_phi, \
        potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_simple(radius_array[-1], longitude=0.1, colatitude=0.1, orbital_frequency=orbital_freq,
                               eccentricity=eccentricity, time=1000.)

    assert type(potential) in [float, np.float, np.float64]
    assert type(potential_partial_theta) in [float, np.float, np.float64]
    assert type(potential_partial_phi) in [float, np.float, np.float64]
    assert type(potential_partial2_theta2) in [float, np.float, np.float64]
    assert type(potential_partial2_phi2) in [float, np.float, np.float64]
    assert type(potential_partial2_theta_phi) in [float, np.float, np.float64]


def test_array():
    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = \
        tidal_potential_simple(radius_array[-1], longitude=np.linspace(0.1, 0.8, 10),
                               colatitude=np.linspace(0.1, 0.8, 10), orbital_frequency=orbital_freq,
                               eccentricity=eccentricity, time=1000.)

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert len(potential) == 10
    assert len(potential_partial_theta) == 10
    assert len(potential_partial_phi) == 10
    assert len(potential_partial2_theta2) == 10
    assert len(potential_partial2_phi2) == 10
    assert len(potential_partial2_theta_phi) == 10