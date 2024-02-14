""" Tests for extracting useful information out of a multilayer tidal propagation
"""
import numpy as np

import TidalPy


from TidalPy.constants import G
from TidalPy.tides.potential import tidal_potential_nsr, tidal_potential_obliquity_nsr, \
    tidal_potential_gen_obliquity_nsr, tidal_potential_simple
from TidalPy.utilities.conversions import orbital_motion2semi_a

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
longitude_array = np.radians(np.linspace(0., 360., 12))
colat_array = np.radians(np.linspace(0.5, 179.5, 13))

volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex128)
host_mass = 50000. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.3
obliquity = np.radians(15.)

time_array = np.linspace(0., 2. * np.pi / orbital_freq, 5)

long_mtx, colat_mtx, time_mtx = np.meshgrid(longitude_array, colat_array, time_array)


def test_tidal_potential_simple():
    """ Test the basic tidal potential assuming low eccentricity, no obliquity, and synchronous rotation """
    # Test with floats
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_simple(
            radius_array[-1], longitude=0.1, colatitude=0.1, time=1000.,
            orbital_frequency=orbital_freq,
            eccentricity=eccentricity,
            host_mass=host_mass, semi_major_axis=semi_major_axis
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) in [float, np.float64]
    assert type(potential_partial_theta) in [float, np.float64]
    assert type(potential_partial_phi) in [float, np.float64]
    assert type(potential_partial2_theta2) in [float, np.float64]
    assert type(potential_partial2_phi2) in [float, np.float64]
    assert type(potential_partial2_theta_phi) in [float, np.float64]

    # Test with matrix
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_simple(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq,
            eccentricity=eccentricity,
            host_mass=host_mass, semi_major_axis=semi_major_axis
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == long_mtx.shape
    assert potential_partial_theta.shape == long_mtx.shape
    assert potential_partial_phi.shape == long_mtx.shape
    assert potential_partial2_theta2.shape == long_mtx.shape
    assert potential_partial2_phi2.shape == long_mtx.shape
    assert potential_partial2_theta_phi.shape == long_mtx.shape

    assert type(potential[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_theta[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_phi[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_phi2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta_phi[0, 0, 0]) in [float, np.float64]


def test_tidal_potential_nsr():
    """ Test the tidal potential equation assuming moderate eccentricity, no obliquity, and non-synchronous rotation """
    # Test with floats; static=False
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=0.1, colatitude=0.1, time=1000.,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == tuple()
    assert potential_partial_theta.shape == tuple()
    assert potential_partial_phi.shape == tuple()
    assert potential_partial2_theta2.shape == tuple()
    assert potential_partial2_phi2.shape == tuple()
    assert potential_partial2_theta_phi.shape == tuple()

    assert potential.dtype in [float, np.float64]
    assert potential_partial_theta.dtype in [float, np.float64]
    assert potential_partial_phi.dtype in [float, np.float64]
    assert potential_partial2_theta2.dtype in [float, np.float64]
    assert potential_partial2_phi2.dtype in [float, np.float64]
    assert potential_partial2_theta_phi.dtype in [float, np.float64]

    # Test with floats; static=True
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=0.1, colatitude=0.1, time=1000.,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=True
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == tuple()
    assert potential_partial_theta.shape == tuple()
    assert potential_partial_phi.shape == tuple()
    assert potential_partial2_theta2.shape == tuple()
    assert potential_partial2_phi2.shape == tuple()
    assert potential_partial2_theta_phi.shape == tuple()

    assert potential.dtype in [float, np.float64]
    assert potential_partial_theta.dtype in [float, np.float64]
    assert potential_partial_phi.dtype in [float, np.float64]
    assert potential_partial2_theta2.dtype in [float, np.float64]
    assert potential_partial2_phi2.dtype in [float, np.float64]
    assert potential_partial2_theta_phi.dtype in [float, np.float64]

    # Test with matrix; static=False
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == long_mtx.shape
    assert potential_partial_theta.shape == long_mtx.shape
    assert potential_partial_phi.shape == long_mtx.shape
    assert potential_partial2_theta2.shape == long_mtx.shape
    assert potential_partial2_phi2.shape == long_mtx.shape
    assert potential_partial2_theta_phi.shape == long_mtx.shape

    assert type(potential[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_theta[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_phi[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_phi2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta_phi[0, 0, 0]) in [float, np.float64]

    # Test with matrix; static=True
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=True
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == long_mtx.shape
    assert potential_partial_theta.shape == long_mtx.shape
    assert potential_partial_phi.shape == long_mtx.shape
    assert potential_partial2_theta2.shape == long_mtx.shape
    assert potential_partial2_phi2.shape == long_mtx.shape
    assert potential_partial2_theta_phi.shape == long_mtx.shape

    assert type(potential[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_theta[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_phi[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_phi2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta_phi[0, 0, 0]) in [float, np.float64]


def test_tidal_potential_obliquity_nsr():
    """ Test the tidal potential equation assuming moderate eccentricity, medium obliquity, and non-synchronous rotation """
    # Test with floats; static=False
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=0.1, colatitude=0.1, time=1000.,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == tuple()
    assert potential_partial_theta.shape == tuple()
    assert potential_partial_phi.shape == tuple()
    assert potential_partial2_theta2.shape == tuple()
    assert potential_partial2_phi2.shape == tuple()
    assert potential_partial2_theta_phi.shape == tuple()

    assert potential.dtype in [float, np.float64]
    assert potential_partial_theta.dtype in [float, np.float64]
    assert potential_partial_phi.dtype in [float, np.float64]
    assert potential_partial2_theta2.dtype in [float, np.float64]
    assert potential_partial2_phi2.dtype in [float, np.float64]
    assert potential_partial2_theta_phi.dtype in [float, np.float64]

    # Test with floats; static=True
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=0.1, colatitude=0.1, time=1000.,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=True
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == tuple()
    assert potential_partial_theta.shape == tuple()
    assert potential_partial_phi.shape == tuple()
    assert potential_partial2_theta2.shape == tuple()
    assert potential_partial2_phi2.shape == tuple()
    assert potential_partial2_theta_phi.shape == tuple()

    assert potential.dtype in [float, np.float64]
    assert potential_partial_theta.dtype in [float, np.float64]
    assert potential_partial_phi.dtype in [float, np.float64]
    assert potential_partial2_theta2.dtype in [float, np.float64]
    assert potential_partial2_phi2.dtype in [float, np.float64]
    assert potential_partial2_theta_phi.dtype in [float, np.float64]

    # Test with matrix; static=False
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == long_mtx.shape
    assert potential_partial_theta.shape == long_mtx.shape
    assert potential_partial_phi.shape == long_mtx.shape
    assert potential_partial2_theta2.shape == long_mtx.shape
    assert potential_partial2_phi2.shape == long_mtx.shape
    assert potential_partial2_theta_phi.shape == long_mtx.shape

    assert type(potential[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_theta[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_phi[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_phi2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta_phi[0, 0, 0]) in [float, np.float64]

    # Test with matrix; static=True
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=True
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == long_mtx.shape
    assert potential_partial_theta.shape == long_mtx.shape
    assert potential_partial_phi.shape == long_mtx.shape
    assert potential_partial2_theta2.shape == long_mtx.shape
    assert potential_partial2_phi2.shape == long_mtx.shape
    assert potential_partial2_theta_phi.shape == long_mtx.shape

    assert type(potential[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_theta[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_phi[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_phi2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta_phi[0, 0, 0]) in [float, np.float64]


def test_tidal_potential_general_obliquity_nsr():
    """ Test the tidal potential equation assuming moderate eccentricity, general obliquity, and non-synchronous rotation """
    # Test with floats; static=False
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_gen_obliquity_nsr(
            radius_array[-1], longitude=0.1, colatitude=0.1, time=1000.,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == tuple()
    assert potential_partial_theta.shape == tuple()
    assert potential_partial_phi.shape == tuple()
    assert potential_partial2_theta2.shape == tuple()
    assert potential_partial2_phi2.shape == tuple()
    assert potential_partial2_theta_phi.shape == tuple()

    assert potential.dtype in [float, np.float64]
    assert potential_partial_theta.dtype in [float, np.float64]
    assert potential_partial_phi.dtype in [float, np.float64]
    assert potential_partial2_theta2.dtype in [float, np.float64]
    assert potential_partial2_phi2.dtype in [float, np.float64]
    assert potential_partial2_theta_phi.dtype in [float, np.float64]

    # Test with floats; static=True
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_gen_obliquity_nsr(
            radius_array[-1], longitude=0.1, colatitude=0.1, time=1000.,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=True
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == tuple()
    assert potential_partial_theta.shape == tuple()
    assert potential_partial_phi.shape == tuple()
    assert potential_partial2_theta2.shape == tuple()
    assert potential_partial2_phi2.shape == tuple()
    assert potential_partial2_theta_phi.shape == tuple()

    assert potential.dtype in [float, np.float64]
    assert potential_partial_theta.dtype in [float, np.float64]
    assert potential_partial_phi.dtype in [float, np.float64]
    assert potential_partial2_theta2.dtype in [float, np.float64]
    assert potential_partial2_phi2.dtype in [float, np.float64]
    assert potential_partial2_theta_phi.dtype in [float, np.float64]

    # Test with matrix; static=False
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_gen_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == long_mtx.shape
    assert potential_partial_theta.shape == long_mtx.shape
    assert potential_partial_phi.shape == long_mtx.shape
    assert potential_partial2_theta2.shape == long_mtx.shape
    assert potential_partial2_phi2.shape == long_mtx.shape
    assert potential_partial2_theta_phi.shape == long_mtx.shape

    assert type(potential[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_theta[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_phi[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_phi2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta_phi[0, 0, 0]) in [float, np.float64]

    # Test with matrix; static=True
    frequencies_by_name, modes_by_name, potential_tuple_by_mode = \
        tidal_potential_gen_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=1.5 * orbital_freq,
            eccentricity=eccentricity, obliquity=obliquity,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=True
            )

    assert len(frequencies_by_name) == 1
    assert len(modes_by_name) == 1
    assert len(potential_tuple_by_mode) == 1
    assert 'n' in frequencies_by_name
    assert 'n' in modes_by_name
    assert 'n' in potential_tuple_by_mode

    potential, potential_partial_theta, potential_partial_phi, \
    potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi = potential_tuple_by_mode['n']

    assert type(potential) == np.ndarray
    assert type(potential_partial_theta) == np.ndarray
    assert type(potential_partial_phi) == np.ndarray
    assert type(potential_partial2_theta2) == np.ndarray
    assert type(potential_partial2_phi2) == np.ndarray
    assert type(potential_partial2_theta_phi) == np.ndarray

    assert potential.shape == long_mtx.shape
    assert potential_partial_theta.shape == long_mtx.shape
    assert potential_partial_phi.shape == long_mtx.shape
    assert potential_partial2_theta2.shape == long_mtx.shape
    assert potential_partial2_phi2.shape == long_mtx.shape
    assert potential_partial2_theta_phi.shape == long_mtx.shape

    assert type(potential[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_theta[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial_phi[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_phi2[0, 0, 0]) in [float, np.float64]
    assert type(potential_partial2_theta_phi[0, 0, 0]) in [float, np.float64]

def test_tidal_potential_simple_vs_nsr():
    """ Test to see if the simple potential and nsr potential can reproduce one another in certain domains."""

    # The simple and nsr potentials should equal to one another for a spin-synchronous world at low eccentricity and
    #   ignoring the static portion of the potential
    eccen = 0.0001
    spin_rate = orbital_freq
    freqs, modes, potentials_simple_by_mode = \
        tidal_potential_simple(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis
            )
    freqs, modes, potentials_nsr_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potentials_simple = potentials_simple_by_mode['n']
    potentials_nsr = potentials_nsr_by_mode['n']

    # Uncomment below to show plot of the two potentials.
    # import matplotlib.pyplot as plt
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_simple[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()
    #
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_nsr[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()

    # Test the potential and its derivatives.
    for i in range(6):
        simple = potentials_simple[i]
        nsr = potentials_nsr[i]
        # There can be differences near 0 that are not well picked up by the percent difference.
        # TODO: what causes these? for now lets set the near 0 values to zero for testing purposes
        simple[np.abs(simple) < 1.e-4] = 0.
        nsr[np.abs(nsr) < 1.e-4] = 0.
        # Have a large rtol due to same issue above.
        assert np.allclose(simple, nsr, rtol=1e-2)

    # The two arrays will not be close if we use a large eccentricity and/or a nsr
    # Test large eccentricity
    eccen = 0.3
    spin_rate = orbital_freq
    freqs, modes, potentials_simple_by_mode = \
        tidal_potential_simple(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis
            )
    freqs, modes, potentials_nsr_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potentials_simple = potentials_simple_by_mode['n']
    potentials_nsr = potentials_nsr_by_mode['n']

    # Uncomment below to show plot of the two potentials.
    # import matplotlib.pyplot as plt
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_simple[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()
    #
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_nsr[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()

    # Test the potential and its derivatives.
    for i in range(6):
        simple = potentials_simple[i]
        nsr = potentials_nsr[i]
        # There can be differences near 0 that are not well picked up by the percent difference.
        # TODO: what causes these? for now lets set the near 0 values to zero for testing purposes
        simple[np.abs(simple) < 1.e-4] = 0.
        nsr[np.abs(nsr) < 1.e-4] = 0.
        # Have a large rtol due to same issue above.
        assert not np.allclose(simple, nsr, rtol=1e-2)

    # Test large nsr
    eccen = 0.0001
    spin_rate = orbital_freq * 1.5
    freqs, modes, potentials_simple_by_mode = \
        tidal_potential_simple(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis
            )
    freqs, modes, potentials_nsr_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potentials_simple = potentials_simple_by_mode['n']
    potentials_nsr = potentials_nsr_by_mode['n']

    # Uncomment below to show plot of the two potentials.
    # import matplotlib.pyplot as plt
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_simple[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()
    #
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_nsr[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()

    # Test the potential and its derivatives.
    for i in range(6):
        simple = potentials_simple[i]
        nsr = potentials_nsr[i]
        # There can be differences near 0 that are not well picked up by the percent difference.
        # TODO: what causes these? for now lets set the near 0 values to zero for testing purposes
        simple[np.abs(simple) < 1.e-4] = 0.
        nsr[np.abs(nsr) < 1.e-4] = 0.
        # Have a large rtol due to same issue above.
        assert not np.allclose(simple, nsr, rtol=1e-2)


def test_tidal_potential_nsr_vs_obliquity():
    """ Test to see if the nsr potential and nsr+obliquity potential can reproduce one another in certain domains."""

    # The nsr and obliquity potentials should equal to one another for obliquity = 0.
    eccen = 0.1
    obliq = np.radians(0.)
    spin_rate = orbital_freq
    freqs, modes, potentials_nsr_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )
    freqs, modes, potentials_obli_by_mode = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen, obliquity=obliq,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potentials_obliquity = potentials_obli_by_mode['n']
    potentials_nsr = potentials_nsr_by_mode['n']

    # Uncomment below to show plot of the two potentials.
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_obliquity[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()
    #
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_nsr[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()

    # Test the potential and its derivatives.
    for i in range(6):
        nsr = potentials_nsr[i]
        obliq_nsr = potentials_obliquity[i]
        assert np.allclose(nsr, obliq_nsr)

    # The two arrays will not be close if we use a non-zero obliquity
    eccen = 0.1
    obliq = np.radians(20.)
    spin_rate = orbital_freq
    freqs, modes, potentials_nsr_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )
    freqs, modes, potentials_obli_by_mode = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen, obliquity=obliq,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potentials_obliquity = potentials_obli_by_mode['n']
    potentials_nsr = potentials_nsr_by_mode['n']

    # Uncomment below to show plot of the two potentials.
    # import matplotlib.pyplot as plt
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_obliquity[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()
    #
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_nsr[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()

    # Test the potential and its derivatives.
    for i in range(6):
        nsr = potentials_nsr[i]
        obliq_nsr = potentials_obliquity[i]
        assert not np.allclose(nsr, obliq_nsr)

    # The two arrays should still be close if obliquity = 0 but there is a NSR
    eccen = 0.1
    obliq = np.radians(0.)
    spin_rate = orbital_freq * 2.5
    freqs, modes, potentials_nsr_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )
    freqs, modes, potentials_obli_by_mode = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen, obliquity=obliq,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potentials_obliquity = potentials_obli_by_mode['n']
    potentials_nsr = potentials_nsr_by_mode['n']

    # Uncomment below to show plot of the two potentials.
    # import matplotlib.pyplot as plt
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_obliquity[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()
    #
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_nsr[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()

    # Test the potential and its derivatives.
    for i in range(6):
        nsr = potentials_nsr[i]
        obliq_nsr = potentials_obliquity[i]
        assert np.allclose(nsr, obliq_nsr)

    # But they will break if NSR and non-zero obliquity
    eccen = 0.1
    obliq = np.radians(20.)
    spin_rate = orbital_freq * 2.5
    freqs, modes, potentials_nsr_by_mode = \
        tidal_potential_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )
    freqs, modes, potentials_obli_by_mode = \
        tidal_potential_obliquity_nsr(
            radius_array[-1], longitude=long_mtx, colatitude=colat_mtx, time=time_mtx,
            orbital_frequency=orbital_freq, rotation_frequency=spin_rate,
            eccentricity=eccen, obliquity=obliq,
            host_mass=host_mass, semi_major_axis=semi_major_axis,
            use_static=False
            )

    potentials_obliquity = potentials_obli_by_mode['n']
    potentials_nsr = potentials_nsr_by_mode['n']

    # Uncomment below to show plot of the two potentials.
    # import matplotlib.pyplot as plt
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_obliquity[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()
    #
    # cbdata = plt.contourf(longitude_array, colat_array, potentials_nsr[0][:, :, 3])
    # plt.colorbar(cbdata)
    # plt.show()

    # Test the potential and its derivatives.
    for i in range(6):
        nsr = potentials_nsr[i]
        obliq_nsr = potentials_obliquity[i]
        assert not np.allclose(nsr, obliq_nsr)
