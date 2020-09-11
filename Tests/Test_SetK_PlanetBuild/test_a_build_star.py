import os
import numpy as np
from scipy.constants import G

import TidalPy

cancri_star_config = {
    "name": "55-Cancri",
    "type": "star",
    "radius": 6.679e8,
    "mass": 1.91e30,
    "luminosity": 2.266e26
}

volume = (4. / 3.) * np.pi * cancri_star_config['radius']**3
surf_area = 4. * np.pi * cancri_star_config['radius']**2
surf_gravity = G * cancri_star_config['mass'] / cancri_star_config['radius']**2

def value_check(built_world, non_numerical: bool = False):

    assert built_world.name == '55-Cancri'
    if not non_numerical:
        np.testing.assert_approx_equal(built_world.mass, cancri_star_config['mass'])
        np.testing.assert_approx_equal(built_world.radius, cancri_star_config['radius'])
        np.testing.assert_approx_equal(built_world.volume, volume)
        np.testing.assert_approx_equal(built_world.density_bulk, cancri_star_config['mass'] / volume)
        np.testing.assert_approx_equal(built_world.surface_area, surf_area)
        np.testing.assert_approx_equal(built_world.radius_inner, 0.)
        np.testing.assert_approx_equal(built_world.thickness, cancri_star_config['radius'])
        np.testing.assert_approx_equal(built_world.gravity_inner, 0.)
        np.testing.assert_approx_equal(built_world.gravity_surf, surf_gravity)

    # Check aliased names
    assert built_world.M is built_world.mass
    assert built_world.R is built_world.radius
    assert built_world.V is built_world.volume
    assert built_world.radius_outer is built_world.radius
    assert built_world.dx is built_world.thickness
    assert built_world.gravity_surface is built_world.gravity_surf
    assert built_world.gravity_outer is built_world.gravity_surf

    return True


def test_build_star_from_manual_config():

    # Test loading a star from a user-provided configuration dictionary
    cancri_star = TidalPy.build_world('55-Cancri', cancri_star_config)

    # Test that its attributes match expectations
    assert value_check(cancri_star)

def test_build_star_from_file_loaded_config():

    # Test loading a star from a user-provided configuration file
    worlds_dir = TidalPy.structures.worldConfigs.WORLD_CONFIG_LOC
    cancri_filepath = os.path.join(worlds_dir, '55cnc.json')
    with open(cancri_filepath, 'r') as cancri_config_file:

        cancri_star = TidalPy.build_world('55-Cancri', cancri_config_file)

    # Test that its attributes match expectations
    assert value_check(cancri_star)

def test_build_star_from_prebuilt():

    # Test loading a star from a pre-built configuration file
    cancri_star = TidalPy.build_world('55-Cancri')

    # The pre-built config may not have the same values as the ones used in this file so we will only perform
    #    non-numerical checks.
    assert value_check(cancri_star, non_numerical=True)