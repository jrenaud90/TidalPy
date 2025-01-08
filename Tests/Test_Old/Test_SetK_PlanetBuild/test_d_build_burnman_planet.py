import os

import pytest
import toml
import numpy as np

from TidalPy.structures import build_world, scale_from_world
from TidalPy.Extending.burnman import burnman_installed

io_config = {
    "name"           : "Io",
    "type"           : "burnman",
    "radius"         : 1821.49e3,
    "orbital_period" : 1.769,
    "eccentricity"   : 0.0041,
    "spin_period"    : 1.769,
    "albedo"         : 0.63,
    "force_spin_sync": True,
    "layers"         : {
        "Core"  : {
            "type"              : "iron",
            "is_tidal"          : False,
            "radius"            : 810.0e3,
            "material"          : ["Pyrite", "Fe_Dewaele"],
            "material_source"   : ["TidalPy", "other"],
            "material_fractions": [0.5, 0.5],
            "temperature_mode"  : "user-defined",
            "temperature_fixed" : 1800.0
            },
        "Mantle": {
            "type"               : "rock",
            "is_tidal"           : True,
            "radius"             : 1821.49e3,
            "material"           : ["forsterite", "mg_perovskite"],
            "material_source"    : ["SLB_2011", "SLB_2011"],
            "material_fractions" : [0.65, 0.35],
            "temperature_mode"   : "adiabatic",
            "temperature_top"    : 1800.0,
            "surface_temperature": 100.0
            }
        }
    }


def planet_asserts(planet):
    assert planet.name == 'Io'
    assert planet.albedo == 0.63
    assert planet.radius == 1821.49e3
    assert planet.mantle.radius == 1821.49e3
    assert planet.core.radius == 810.0e3

    return True


def test_build_burnman_world_from_dict():

    if not burnman_installed:
        pytest.skip('Burnman not installed. Skipping test.')

    result = build_world('Io', world_config=io_config)

    assert planet_asserts(result)


def test_build_burnman_world_from_toml():

    if not burnman_installed:
        pytest.skip('Burnman not installed. Skipping test.')

    with open('temp_io.toml', 'w') as planet_file:
        toml.dump(io_config, planet_file)

    with open('temp_io.toml', 'r') as planet_file:
        result = build_world('Io', world_config=planet_file)

    os.remove('temp_io.toml')

    assert planet_asserts(result)


def test_build_from_burnman_world_scale_radius():
    """ This will test building a secondary BurnmanWorld from an already built one using a radius scaling technique. """

    if not burnman_installed:
        pytest.skip('Burnman not installed. Skipping test.')

    # Build a regular layered Io first.
    io = build_world('Io')

    # Now build a second version that only has a name and eccentricity change.
    io_2 = scale_from_world(io, new_name='small_io', radius_scale=0.5)

    # Check name
    assert io_2.name == 'small_io'

    # Check values
    np.testing.assert_approx_equal(io_2.radius, io.radius / 2.)
    for layer, layer_2 in zip(io, io_2):
        np.testing.assert_approx_equal(layer_2.radius, layer.radius / 2.)

    # Check for scaling a larger radius
    io_3 = scale_from_world(io, new_name='large_io', radius_scale=2.)

    # Check name
    assert io_3.name == 'large_io'

    # Check values
    np.testing.assert_approx_equal(io_3.radius, io.radius * 2.)
    for layer, layer_3 in zip(io, io_3):
        np.testing.assert_approx_equal(layer_3.radius, layer.radius * 2.)


def test_build_burnman_world_from_prebuilt():

    if not burnman_installed:
        pytest.skip('Burnman not installed. Skipping test.')

    result = build_world('Io')

    assert planet_asserts(result)


def test_paint_burnman():

    if not burnman_installed:
        pytest.skip('Burnman not installed. Skipping test.')
        
    # This will test the slicing features of the world_types, as well as the painting graphics tool
    io = build_world('Io')
    assert io.paint(auto_show=False)