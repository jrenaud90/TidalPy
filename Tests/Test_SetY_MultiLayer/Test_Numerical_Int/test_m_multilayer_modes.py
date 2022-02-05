""" Tests for the multilayer mode calculator module

"""

import numpy as np

import TidalPy
from TidalPy.constants import G, mass_trap1
from TidalPy.toolbox.conversions import orbital_motion2semi_a
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.toolbox.multilayer import calculate_homogen_solid, calculate_sls, calculate_ssls
from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes
from TidalPy.tides.modes.multilayer_modes import collapse_multilayer_modes

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
R = radius_array[-1]
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
host_mass = mass_trap1
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.float64)
viscosity_array = 1.0e19 * np.ones(10, dtype=np.float64)
bulk_array = 10.e10 * np.ones(10, dtype=np.float64)
radius_array = radius_array[1:]
orbital_frequency = 2. * np.pi / (86400. * 1.)
spin_frequency = 1.5 * orbital_frequency
semi_major_axis = orbital_motion2semi_a(orbital_frequency, host_mass, planet_mass)
eccentricity = 0.2


# Build other domains and the final multidimensional arrays
colatitude = np.radians(np.linspace(0.1, 179.9, 5))
longitude = np.radians(np.linspace(0., 360., 7))
time = np.linspace(0., 2. * np.pi / orbital_frequency, 12)

# Figure out volume for each voxel
voxel_volumes = calculate_voxel_volumes(radius_array, longitude, colatitude)

radius_matrix, longitude_matrix, colatitude_matrix, time_matrix = \
    np.meshgrid(radius_array, longitude, colatitude, time, indexing='ij')


def test_collapse_multilayer_modes_OrbAvg_homogen():
    """ Test multilayer mode collapse function for a homogeneous world with orbit averaging on. """

    interface_properties = {'use_static': False, 'use_kamata': True, }
    integration_parameters = {'use_julia': False}

    heating, volumetric_heating, volumetric_heating_by_mode, strains, stresses = \
        collapse_multilayer_modes(shear_array, viscosity_array, bulk_array, density_array, gravity_array,
                              maxwell, radius_matrix, longitude_matrix, colatitude_matrix, time_matrix, voxel_volumes,
                              orbital_frequency, spin_frequency, semi_major_axis, eccentricity, host_mass,
                              calculate_homogen_solid, interface_properties=interface_properties, order_l=2,
                              complex_compliance_input=None,
                              orbit_average_results=True, integration_parameters=integration_parameters)

    # Check types
    assert type(heating) == np.ndarray
    assert type(volumetric_heating) == np.ndarray
    assert type(volumetric_heating_by_mode) == dict
    assert type(strains) == np.ndarray
    assert type(stresses) == np.ndarray
    assert type(heating[0, 0, 0]) in (np.float64, np.float, float)
    assert type(volumetric_heating[0, 0, 0]) in (np.float64, np.float, float)

    # Check shape
    assert heating.shape == (10, 7, 5)
    assert volumetric_heating.shape == (10, 7, 5)
    assert strains.shape == (6, 10, 7, 5)
    assert stresses.shape == (6, 10, 7, 5)

    # Look at stress/strain shapes
    for i in range(6):
        stress_ = stresses[i, :, :, :]
        strain_ = strains[i, :, :, :]

        # Check type
        assert type(stress_) == np.ndarray
        assert type(strain_) == np.ndarray
        assert type(stress_[0, 0, 0]) in (np.complex128, np.complex, complex)
        assert type(strain_[0, 0, 0]) in (np.complex128, np.complex, complex)

        # Check shape
        assert stress_.shape == (10, 7, 5)
        assert strain_.shape == (10, 7, 5)

    # Look at heating per mode
    for mode_name, heating_at_mode in volumetric_heating_by_mode.items():
        # Check types
        assert type(mode_name) == str
        assert type(heating_at_mode) == np.ndarray
        assert type(heating_at_mode[0, 0, 0]) in (np.float64, np.float, float)

        # Check shape
        assert heating_at_mode.shape == (10, 7, 5)


def test_collapse_multilayer_modes_NoOrbAvg_homogen():
    """ Test multilayer mode collapse function for a homogeneous world with orbit averaging off. """

    interface_properties = {'use_static': False, 'use_kamata': True, }
    integration_parameters = {'use_julia': False}

    heating, volumetric_heating, volumetric_heating_by_mode, strains, stresses = \
        collapse_multilayer_modes(shear_array, viscosity_array, bulk_array, density_array, gravity_array,
                              maxwell, radius_matrix, longitude_matrix, colatitude_matrix, time_matrix, voxel_volumes,
                              orbital_frequency, spin_frequency, semi_major_axis, eccentricity, host_mass,
                              calculate_homogen_solid, interface_properties=interface_properties, order_l=2,
                              complex_compliance_input=None,
                              orbit_average_results=False, integration_parameters=integration_parameters)

    # Check types
    assert type(heating) == np.ndarray
    assert type(volumetric_heating) == np.ndarray
    assert type(volumetric_heating_by_mode) == dict
    assert type(strains) == np.ndarray
    assert type(stresses) == np.ndarray
    assert type(heating[0, 0, 0, 0]) in (np.float64, np.float, float)
    assert type(volumetric_heating[0, 0, 0, 0]) in (np.float64, np.float, float)

    # Check shape
    assert heating.shape == (10, 7, 5, 12)
    assert volumetric_heating.shape == (10, 7, 5, 12)
    assert strains.shape == (6, 10, 7, 5, 12)
    assert stresses.shape == (6, 10, 7, 5, 12)

    # Look at stress/strain shapes
    for i in range(6):
        stress_ = stresses[i, :, :, :, :]
        strain_ = strains[i, :, :, :, :]

        # Check type
        assert type(stress_) == np.ndarray
        assert type(strain_) == np.ndarray
        assert type(stress_[0, 0, 0, 0]) in (np.complex128, np.complex, complex)
        assert type(strain_[0, 0, 0, 0]) in (np.complex128, np.complex, complex)

        # Check shape
        assert stress_.shape == (10, 7, 5, 12)
        assert strain_.shape == (10, 7, 5, 12)

    # Look at heating per mode
    for mode_name, heating_at_mode in volumetric_heating_by_mode.items():
        # Check types
        assert type(mode_name) == str
        assert type(heating_at_mode) == np.ndarray
        assert type(heating_at_mode[0, 0, 0, 0]) in (np.float64, np.float, float)

        # Check shape
        assert heating_at_mode.shape == (10, 7, 5, 12)