""" Tests for the multilayer mode calculator module

"""
import numpy as np

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G, mass_trap1
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.tides.modes.multilayer_modes import collapse_multilayer_modes
from TidalPy.toolbox.conversions import orbital_motion2semi_a
from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes

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
planet_bulk_density = planet_mass / np.sum(volume_array)
spin_frequency = 1. * orbital_frequency
semi_major_axis = orbital_motion2semi_a(orbital_frequency, host_mass, planet_mass)
eccentricity = 0.05
obliquity = None
layer_i = radius_array > 0.

# Build other domains and the final multidimensional arrays
colatitude = np.radians(np.linspace(0.1, 179.9, 5))
longitude = np.radians(np.linspace(0., 360., 7))
time = np.linspace(0., 2. * np.pi / orbital_frequency, 12)

# Figure out volume for each voxel
voxel_volumes = calculate_voxel_volumes(radius_array, longitude, colatitude)

longitude_matrix, colatitude_matrix, time_matrix = \
    np.meshgrid(longitude, colatitude, time, indexing='ij')

tidal_y_int_kwargs = {
    'use_kamata'         : False,
    'planet_bulk_density': planet_bulk_density,
    'nondimensionalize'  : False,
    'integrator'         : 'scipy',
    'integration_method' : None,
    'integration_rtol'   : 1.0e-6,
    'integration_atol'   : 1.0e-8
    }

input_kwargs = {
    'orbital_frequency'          : orbital_frequency, 'spin_frequency': spin_frequency,
    'semi_major_axis'            : semi_major_axis, 'eccentricity': eccentricity, 'host_mass': host_mass,
    'radius_array'               : radius_array, 'shear_array': shear_array, 'bulk_array': bulk_array,
    'viscosity_array'            : viscosity_array, 'density_array': density_array, 'gravity_array': gravity_array,
    'longitude_matrix'           : longitude_matrix, 'colatitude_matrix': colatitude_matrix,
    'time_matrix'                : time_matrix, 'voxel_volume': voxel_volumes, 'complex_compliance_function': maxwell,
    'is_solid_by_layer'          : [True], 'is_static_by_layer': [False], 'indices_by_layer': [layer_i],
    'obliquity'                  : obliquity,
    'surface_boundary_conditions': None, 'solve_load_numbers': False,
    'complex_compliance_input'   : None, 'force_mode_calculation': False,
    'order_l'                    : 2,
    'use_modes'                  : False, 'use_static_potential': False, 'use_simple_potential': False,
    'orbit_average_results'      : False,
    **tidal_y_int_kwargs
    }


def test_collapse_multilayer_modes_homogen_modesoff_noorbitavg():
    """ Test the multilayer multimode calculator and collapse function.
    For this function we assume a homogeneous planet with only 1 mode and no orbit averaging. """

    # Test non-orbit averaged
    input_kwargs_to_use = {**input_kwargs}
    heating, volumetric_heating, strains, stresses, \
    total_potential, tidal_potential, complex_shears_avg, tidal_y_avg, \
    (love_k_by_mode, love_h_by_mode, love_l_by_mode), tidal_modes, modes_skipped = \
        collapse_multilayer_modes(**input_kwargs_to_use)

    expected_shape = (*radius_array.shape, *longitude_matrix.shape)
    assert heating.shape == expected_shape
    assert volumetric_heating.shape == expected_shape
    assert total_potential.shape == expected_shape
    assert tidal_potential.shape == colatitude_matrix.shape
    assert strains.shape == (6, *expected_shape)
    assert stresses.shape == (6, *expected_shape)
    assert complex_shears_avg.shape == radius_array.shape
    assert tidal_y_avg.shape == (6, *radius_array.shape)

    assert heating.dtype in [float, np.float64]
    assert volumetric_heating.dtype in [float, np.float64]
    assert total_potential.dtype in [complex, np.complex128]
    assert tidal_potential.dtype in [complex, np.complex128]
    assert strains.dtype in [complex, np.complex128]
    assert stresses.dtype in [complex, np.complex128]
    assert complex_shears_avg.dtype in [complex, np.complex128]
    assert tidal_y_avg.dtype in [complex, np.complex128]

    assert type(love_k_by_mode) == dict
    assert type(love_h_by_mode) == dict
    assert type(love_l_by_mode) == dict

    # For modes off we expect there top be only one mode
    assert len(love_k_by_mode) == 1
    assert len(love_h_by_mode) == 1
    assert len(love_l_by_mode) == 1

    # and it should be 'n'
    assert 'n' in love_k_by_mode
    assert 'n' in love_h_by_mode
    assert 'n' in love_l_by_mode
    assert 'n' in tidal_modes

    assert type(love_k_by_mode['n']) in [complex, np.complex128]
    assert type(love_h_by_mode['n']) in [complex, np.complex128]
    assert type(love_l_by_mode['n']) in [complex, np.complex128]
    assert type(tidal_modes['n']) in [float, np.float64]

    # No modes should have been skipped
    assert type(modes_skipped) == dict
    assert len(modes_skipped) == 0


def test_collapse_multilayer_modes_homogen_modesoff():
    """ Test the multilayer multimode calculator and collapse function.
    For this function we assume a homogeneous planet with only 1 mode.
    Does include: orbit averaging. """

    # Test non-orbit averaged
    input_kwargs_to_use = {**input_kwargs}
    input_kwargs_to_use['orbit_average_results'] = True
    heating, volumetric_heating, strains, stresses, \
    total_potential, tidal_potential, complex_shears_avg, tidal_y_avg, \
    (love_k_by_mode, love_h_by_mode, love_l_by_mode), tidal_modes, modes_skipped = \
        collapse_multilayer_modes(**input_kwargs_to_use)

    expected_shape = (*radius_array.shape, *longitude_matrix.shape[:2])
    assert heating.shape == expected_shape
    assert volumetric_heating.shape == expected_shape
    assert total_potential.shape == expected_shape
    assert tidal_potential.shape == colatitude_matrix.shape[:2]
    assert strains.shape == (6, *expected_shape)
    assert stresses.shape == (6, *expected_shape)
    assert complex_shears_avg.shape == radius_array.shape
    assert tidal_y_avg.shape == (6, *radius_array.shape)

    assert heating.dtype in [float, np.float64]
    assert volumetric_heating.dtype in [float, np.float64]
    assert total_potential.dtype in [complex, np.complex128]
    assert tidal_potential.dtype in [complex, np.complex128]
    assert strains.dtype in [complex, np.complex128]
    assert stresses.dtype in [complex, np.complex128]
    assert complex_shears_avg.dtype in [complex, np.complex128]
    assert tidal_y_avg.dtype in [complex, np.complex128]

    assert type(love_k_by_mode) == dict
    assert type(love_h_by_mode) == dict
    assert type(love_l_by_mode) == dict

    # For modes off we expect there top be only one mode
    assert len(love_k_by_mode) == 1
    assert len(love_h_by_mode) == 1
    assert len(love_l_by_mode) == 1

    # and it should be 'n'
    assert 'n' in love_k_by_mode
    assert 'n' in love_h_by_mode
    assert 'n' in love_l_by_mode
    assert 'n' in tidal_modes

    assert type(love_k_by_mode['n']) in [complex, np.complex128]
    assert type(love_h_by_mode['n']) in [complex, np.complex128]
    assert type(love_l_by_mode['n']) in [complex, np.complex128]
    assert type(tidal_modes['n']) in [float, np.float64]

    # No modes should have been skipped
    assert type(modes_skipped) == dict
    assert len(modes_skipped) == 0


def test_collapse_multilayer_modes_homogen_modeson():
    """ Test the multilayer multimode calculator and collapse function.
    For this function we assume a homogeneous planet.
    Does include: orbit averaging, multiple modes"""

    # Test non-orbit averaged
    input_kwargs_to_use = {**input_kwargs}
    input_kwargs_to_use['orbit_average_results'] = True
    input_kwargs_to_use['use_modes'] = True
    # Assume spin-sync
    input_kwargs_to_use['spin_frequency'] = input_kwargs_to_use['orbital_frequency']

    heating, volumetric_heating, strains, stresses, \
    total_potential, tidal_potential, complex_shears_avg, tidal_y_avg, \
    (love_k_by_mode, love_h_by_mode, love_l_by_mode), tidal_modes, modes_skipped = \
        collapse_multilayer_modes(**input_kwargs_to_use)

    expected_shape = (*radius_array.shape, *longitude_matrix.shape[:2])
    assert heating.shape == expected_shape
    assert volumetric_heating.shape == expected_shape
    assert total_potential.shape == expected_shape
    assert tidal_potential.shape == colatitude_matrix.shape[:2]
    assert strains.shape == (6, *expected_shape)
    assert stresses.shape == (6, *expected_shape)
    assert complex_shears_avg.shape == radius_array.shape
    assert tidal_y_avg.shape == (6, *radius_array.shape)

    assert heating.dtype in [float, np.float64]
    assert volumetric_heating.dtype in [float, np.float64]
    assert total_potential.dtype in [complex, np.complex128]
    assert tidal_potential.dtype in [complex, np.complex128]
    assert strains.dtype in [complex, np.complex128]
    assert stresses.dtype in [complex, np.complex128]
    assert complex_shears_avg.dtype in [complex, np.complex128]
    assert tidal_y_avg.dtype in [complex, np.complex128]

    assert type(love_k_by_mode) == dict
    assert type(love_h_by_mode) == dict
    assert type(love_l_by_mode) == dict

    # For modes off we expect there top be only one mode
    modes_at_these_assumptions = ('n', '2n', '3n', '2o+n', '2o-n', '2o-2n', '2o-3n', '2o-4n', '2o-5n')
    assert len(love_k_by_mode) == len(modes_at_these_assumptions)
    assert len(love_h_by_mode) == len(modes_at_these_assumptions)
    assert len(love_l_by_mode) == len(modes_at_these_assumptions)

    # and it should be 'n'
    for mode in modes_at_these_assumptions:
        assert mode in love_k_by_mode
        assert mode in love_h_by_mode
        assert mode in love_l_by_mode
        assert mode in tidal_modes

        assert type(love_k_by_mode[mode]) in [complex, np.complex128]
        assert type(love_h_by_mode[mode]) in [complex, np.complex128]
        assert type(love_l_by_mode[mode]) in [complex, np.complex128]
        assert type(tidal_modes[mode]) in [float, np.float64]

    # For spin synch, one mode should have been skipped
    assert type(modes_skipped) == dict
    assert len(modes_skipped) == 1


def test_collapse_multilayer_modes_homogen_modeson_nonsor():
    """ Test the multilayer multimode calculator and collapse function.
    For this function we assume a homogeneous planet.
    Does include: orbit averaging, multiple modes, planet is not in a spin-orbit resonance"""

    # Test non-orbit averaged
    input_kwargs_to_use = {**input_kwargs}
    input_kwargs_to_use['orbit_average_results'] = True
    input_kwargs_to_use['use_modes'] = True
    # Assume spin-sync
    input_kwargs_to_use['spin_frequency'] = 0.7 * input_kwargs_to_use['orbital_frequency']

    heating, volumetric_heating, strains, stresses, \
    total_potential, tidal_potential, complex_shears_avg, tidal_y_avg, \
    (love_k_by_mode, love_h_by_mode, love_l_by_mode), tidal_modes, modes_skipped = \
        collapse_multilayer_modes(**input_kwargs_to_use)

    expected_shape = (*radius_array.shape, *longitude_matrix.shape[:2])
    assert heating.shape == expected_shape
    assert volumetric_heating.shape == expected_shape
    assert total_potential.shape == expected_shape
    assert tidal_potential.shape == colatitude_matrix.shape[:2]
    assert strains.shape == (6, *expected_shape)
    assert stresses.shape == (6, *expected_shape)
    assert complex_shears_avg.shape == radius_array.shape
    assert tidal_y_avg.shape == (6, *radius_array.shape)

    assert heating.dtype in [float, np.float64]
    assert volumetric_heating.dtype in [float, np.float64]
    assert total_potential.dtype in [complex, np.complex128]
    assert tidal_potential.dtype in [complex, np.complex128]
    assert strains.dtype in [complex, np.complex128]
    assert stresses.dtype in [complex, np.complex128]
    assert complex_shears_avg.dtype in [complex, np.complex128]
    assert tidal_y_avg.dtype in [complex, np.complex128]

    assert type(love_k_by_mode) == dict
    assert type(love_h_by_mode) == dict
    assert type(love_l_by_mode) == dict

    # For modes off we expect there top be only one mode
    modes_at_these_assumptions = ('n', '2n', '3n', '2o+n', '2o-n', '2o-2n', '2o-3n', '2o-4n', '2o-5n')
    assert len(love_k_by_mode) == len(modes_at_these_assumptions)
    assert len(love_h_by_mode) == len(modes_at_these_assumptions)
    assert len(love_l_by_mode) == len(modes_at_these_assumptions)

    # and it should be 'n'
    for mode in modes_at_these_assumptions:
        assert mode in love_k_by_mode
        assert mode in love_h_by_mode
        assert mode in love_l_by_mode
        assert mode in tidal_modes

        assert type(love_k_by_mode[mode]) in [complex, np.complex128]
        assert type(love_h_by_mode[mode]) in [complex, np.complex128]
        assert type(love_l_by_mode[mode]) in [complex, np.complex128]
        assert type(tidal_modes[mode]) in [float, np.float64]

    # For non-SOR, no modes should have been skipped
    assert type(modes_skipped) == dict
    assert len(modes_skipped) == 0


def test_collapse_multilayer_modes_homogen_modeson_negative_spin():
    """ Test the multilayer multimode calculator and collapse function.
    For this function we assume a homogeneous planet.
    Does include: orbit averaging, multiple modes, planet has negative spin """

    # Test non-orbit averaged
    input_kwargs_to_use = {**input_kwargs}
    input_kwargs_to_use['orbit_average_results'] = True
    input_kwargs_to_use['use_modes'] = True
    # Assume spin-sync
    input_kwargs_to_use['spin_frequency'] = -1. * input_kwargs_to_use['orbital_frequency']

    heating, volumetric_heating, strains, stresses, \
    total_potential, tidal_potential, complex_shears_avg, tidal_y_avg, \
    (love_k_by_mode, love_h_by_mode, love_l_by_mode), tidal_modes, modes_skipped = \
        collapse_multilayer_modes(**input_kwargs_to_use)

    expected_shape = (*radius_array.shape, *longitude_matrix.shape[:2])
    assert heating.shape == expected_shape
    assert volumetric_heating.shape == expected_shape
    assert total_potential.shape == expected_shape
    assert tidal_potential.shape == colatitude_matrix.shape[:2]
    assert strains.shape == (6, *expected_shape)
    assert stresses.shape == (6, *expected_shape)
    assert complex_shears_avg.shape == radius_array.shape
    assert tidal_y_avg.shape == (6, *radius_array.shape)

    assert heating.dtype in [float, np.float64]
    assert volumetric_heating.dtype in [float, np.float64]
    assert total_potential.dtype in [complex, np.complex128]
    assert tidal_potential.dtype in [complex, np.complex128]
    assert strains.dtype in [complex, np.complex128]
    assert stresses.dtype in [complex, np.complex128]
    assert complex_shears_avg.dtype in [complex, np.complex128]
    assert tidal_y_avg.dtype in [complex, np.complex128]

    assert type(love_k_by_mode) == dict
    assert type(love_h_by_mode) == dict
    assert type(love_l_by_mode) == dict

    # For modes off we expect there top be only one mode
    modes_at_these_assumptions = ('n', '2n', '3n', '2o+n', '2o-n', '2o-2n', '2o-3n', '2o-4n', '2o-5n')
    assert len(love_k_by_mode) == len(modes_at_these_assumptions)
    assert len(love_h_by_mode) == len(modes_at_these_assumptions)
    assert len(love_l_by_mode) == len(modes_at_these_assumptions)

    # and it should be 'n'
    for mode in modes_at_these_assumptions:
        assert mode in love_k_by_mode
        assert mode in love_h_by_mode
        assert mode in love_l_by_mode
        assert mode in tidal_modes

        assert type(love_k_by_mode[mode]) in [complex, np.complex128]
        assert type(love_h_by_mode[mode]) in [complex, np.complex128]
        assert type(love_l_by_mode[mode]) in [complex, np.complex128]
        assert type(tidal_modes[mode]) in [float, np.float64]

    # For O=-n, no modes should have been skipped
    assert type(modes_skipped) == dict
    assert len(modes_skipped) == 0


def test_collapse_multilayer_modes_homogen_modeson_negative_n():
    """ Test the multilayer multimode calculator and collapse function.
    For this function we assume a homogeneous planet.
    Does include: orbit averaging, multiple modes, planet has negative orbital motion """

    # Test non-orbit averaged
    input_kwargs_to_use = {**input_kwargs}
    input_kwargs_to_use['orbit_average_results'] = True
    input_kwargs_to_use['use_modes'] = True
    # Assume spin-sync
    input_kwargs_to_use['orbital_frequency'] *= -1.
    input_kwargs_to_use['spin_frequency'] = input_kwargs_to_use['orbital_frequency']

    heating, volumetric_heating, strains, stresses, \
    total_potential, tidal_potential, complex_shears_avg, tidal_y_avg, \
    (love_k_by_mode, love_h_by_mode, love_l_by_mode), tidal_modes, modes_skipped = \
        collapse_multilayer_modes(**input_kwargs_to_use)

    expected_shape = (*radius_array.shape, *longitude_matrix.shape[:2])
    assert heating.shape == expected_shape
    assert volumetric_heating.shape == expected_shape
    assert total_potential.shape == expected_shape
    assert tidal_potential.shape == colatitude_matrix.shape[:2]
    assert strains.shape == (6, *expected_shape)
    assert stresses.shape == (6, *expected_shape)
    assert complex_shears_avg.shape == radius_array.shape
    assert tidal_y_avg.shape == (6, *radius_array.shape)

    assert heating.dtype in [float, np.float64]
    assert volumetric_heating.dtype in [float, np.float64]
    assert total_potential.dtype in [complex, np.complex128]
    assert tidal_potential.dtype in [complex, np.complex128]
    assert strains.dtype in [complex, np.complex128]
    assert stresses.dtype in [complex, np.complex128]
    assert complex_shears_avg.dtype in [complex, np.complex128]
    assert tidal_y_avg.dtype in [complex, np.complex128]

    assert type(love_k_by_mode) == dict
    assert type(love_h_by_mode) == dict
    assert type(love_l_by_mode) == dict

    # For modes off we expect there top be only one mode
    modes_at_these_assumptions = ('n', '2n', '3n', '2o+n', '2o-n', '2o-2n', '2o-3n', '2o-4n', '2o-5n')
    assert len(love_k_by_mode) == len(modes_at_these_assumptions)
    assert len(love_h_by_mode) == len(modes_at_these_assumptions)
    assert len(love_l_by_mode) == len(modes_at_these_assumptions)

    # and it should be 'n'
    for mode in modes_at_these_assumptions:
        assert mode in love_k_by_mode
        assert mode in love_h_by_mode
        assert mode in love_l_by_mode
        assert mode in tidal_modes

        assert type(love_k_by_mode[mode]) in [complex, np.complex128]
        assert type(love_h_by_mode[mode]) in [complex, np.complex128]
        assert type(love_l_by_mode[mode]) in [complex, np.complex128]
        assert type(tidal_modes[mode]) in [float, np.float64]

    # This is technically spin synch just with negative spin and orb motion, one mode should have been skipped
    assert type(modes_skipped) == dict
    assert len(modes_skipped) == 1


def test_collapse_multilayer_modes_homogen_modeson_obliquity():
    """ Test the multilayer multimode calculator and collapse function.
    For this function we assume a homogeneous planet.
    Does include: orbit averaging, multiple modes, non-zero obliquity"""

    # Test non-orbit averaged
    input_kwargs_to_use = {**input_kwargs}
    input_kwargs_to_use['orbit_average_results'] = True
    input_kwargs_to_use['use_modes'] = True
    # Assume spin-sync
    input_kwargs_to_use['spin_frequency'] = input_kwargs_to_use['orbital_frequency']
    # Assume non-zero obliquity
    input_kwargs_to_use['obliquity'] = np.radians(45.)

    heating, volumetric_heating, strains, stresses, \
    total_potential, tidal_potential, complex_shears_avg, tidal_y_avg, \
    (love_k_by_mode, love_h_by_mode, love_l_by_mode), tidal_modes, modes_skipped = \
        collapse_multilayer_modes(**input_kwargs_to_use)

    expected_shape = (*radius_array.shape, *longitude_matrix.shape[:2])
    assert heating.shape == expected_shape
    assert volumetric_heating.shape == expected_shape
    assert total_potential.shape == expected_shape
    assert tidal_potential.shape == colatitude_matrix.shape[:2]
    assert strains.shape == (6, *expected_shape)
    assert stresses.shape == (6, *expected_shape)
    assert complex_shears_avg.shape == radius_array.shape
    assert tidal_y_avg.shape == (6, *radius_array.shape)

    assert heating.dtype in [float, np.float64]
    assert volumetric_heating.dtype in [float, np.float64]
    assert total_potential.dtype in [complex, np.complex128]
    assert tidal_potential.dtype in [complex, np.complex128]
    assert strains.dtype in [complex, np.complex128]
    assert stresses.dtype in [complex, np.complex128]
    assert complex_shears_avg.dtype in [complex, np.complex128]
    assert tidal_y_avg.dtype in [complex, np.complex128]

    assert type(love_k_by_mode) == dict
    assert type(love_h_by_mode) == dict
    assert type(love_l_by_mode) == dict

    # For modes off we expect there top be only one mode
    modes_at_these_assumptions = (
        'n', '2n', '3n', '2o+n', '2o-n', '2o-2n', '2o-3n', '2o-4n', '2o-5n', 'o', '2o', 'o+n', 'o+2n', 'o-n', 'o-2n',
        'o-3n', 'o-4n')
    assert len(love_k_by_mode) == len(modes_at_these_assumptions)
    assert len(love_h_by_mode) == len(modes_at_these_assumptions)
    assert len(love_l_by_mode) == len(modes_at_these_assumptions)

    # and it should be 'n'
    for mode in modes_at_these_assumptions:
        assert mode in love_k_by_mode
        assert mode in love_h_by_mode
        assert mode in love_l_by_mode
        assert mode in tidal_modes

        assert type(love_k_by_mode[mode]) in [complex, np.complex128]
        assert type(love_h_by_mode[mode]) in [complex, np.complex128]
        assert type(love_l_by_mode[mode]) in [complex, np.complex128]
        assert type(tidal_modes[mode]) in [float, np.float64]

    # For spin synch, two modes should have been skipped
    assert type(modes_skipped) == dict
    assert len(modes_skipped) == 2


def test_collapse_multilayer_modes_liquid_solid():
    """ Test the multilayer multimode calculator and collapse function.
    For this function we assume a planet with a liquid core and solid mantle
    Does include: orbit averaging, multiple modes"""

    # Test non-orbit averaged
    input_kwargs_to_use = {**input_kwargs}
    input_kwargs_to_use['orbit_average_results'] = True
    input_kwargs_to_use['use_modes'] = True
    # Assume spin-sync
    input_kwargs_to_use['spin_frequency'] = input_kwargs_to_use['orbital_frequency']
    # Create new structure
    r_core = radius_array[-1] / 2.
    core_index = radius_array <= r_core
    mantle_index = radius_array > r_core
    input_kwargs_to_use['shear_array'][core_index] = 0.
    input_kwargs_to_use['is_solid_by_layer'] = [False, True]
    input_kwargs_to_use['is_static_by_layer'] = [True, False]
    input_kwargs_to_use['indices_by_layer'] = [core_index, mantle_index]

    heating, volumetric_heating, strains, stresses, \
    total_potential, tidal_potential, complex_shears_avg, tidal_y_avg, \
    (love_k_by_mode, love_h_by_mode, love_l_by_mode), tidal_modes, modes_skipped = \
        collapse_multilayer_modes(**input_kwargs_to_use)

    expected_shape = (*radius_array.shape, *longitude_matrix.shape[:2])
    assert heating.shape == expected_shape
    assert volumetric_heating.shape == expected_shape
    assert total_potential.shape == expected_shape
    assert tidal_potential.shape == colatitude_matrix.shape[:2]
    assert strains.shape == (6, *expected_shape)
    assert stresses.shape == (6, *expected_shape)
    assert complex_shears_avg.shape == radius_array.shape
    assert tidal_y_avg.shape == (6, *radius_array.shape)

    assert heating.dtype in [float, np.float64]
    assert volumetric_heating.dtype in [float, np.float64]
    assert total_potential.dtype in [complex, np.complex128]
    assert tidal_potential.dtype in [complex, np.complex128]
    assert strains.dtype in [complex, np.complex128]
    assert stresses.dtype in [complex, np.complex128]
    assert complex_shears_avg.dtype in [complex, np.complex128]
    assert tidal_y_avg.dtype in [complex, np.complex128]

    assert type(love_k_by_mode) == dict
    assert type(love_h_by_mode) == dict
    assert type(love_l_by_mode) == dict

    # For modes off we expect there top be only one mode
    modes_at_these_assumptions = ('n', '2n', '3n', '2o+n', '2o-n', '2o-2n', '2o-3n', '2o-4n', '2o-5n')
    assert len(love_k_by_mode) == len(modes_at_these_assumptions)
    assert len(love_h_by_mode) == len(modes_at_these_assumptions)
    assert len(love_l_by_mode) == len(modes_at_these_assumptions)

    # and it should be 'n'
    for mode in modes_at_these_assumptions:
        assert mode in love_k_by_mode
        assert mode in love_h_by_mode
        assert mode in love_l_by_mode
        assert mode in tidal_modes

        assert type(love_k_by_mode[mode]) in [complex, np.complex128]
        assert type(love_h_by_mode[mode]) in [complex, np.complex128]
        assert type(love_l_by_mode[mode]) in [complex, np.complex128]
        assert type(tidal_modes[mode]) in [float, np.float64]

    # For spin synch, one mode should have been skipped
    assert type(modes_skipped) == dict
    assert len(modes_skipped) == 1
