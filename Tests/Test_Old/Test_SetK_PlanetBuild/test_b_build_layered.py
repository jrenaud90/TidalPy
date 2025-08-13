import numpy as np

import TidalPy



from TidalPy.structures import build_world, build_from_world, scale_from_world

io_config = {
    "name"           : "Io",
    "type"           : "layered",
    "radius"         : 1821.49e3,
    "orbital_period" : 1.769,
    "eccentricity"   : 0.0041,
    "spin_period"    : 1.769,
    "albedo"         : 0.63,
    "force_spin_sync": True,
    "layers"         : {
        "Core"  : {
            "type"    : "iron",
            "is_tidal": False,
            "radius"  : 810.0e3,
            "density" : 5200.
            },
        "Mantle": {
            "type"               : "rock",
            "is_tidal"           : True,
            "radius"             : 1821.49e3,
            "surface_temperature": 100.0,
            "density"            : 3200.
            }
        }
    }


def value_check(world, config_to_compare: dict, check_name: bool = True):
    """ This is a helper function to test the values of a LayeredWorld. Used in LayeredWorld builder tests. """

    if check_name:
        assert world.name == config_to_compare['name']
    assert world.num_layers == 2
    assert len(world.layers) == 2
    # Note: <world>.layers_by_name will usually contain more than the expected number of layers due to allowing both
    #    title and lowercase keys (number of layer instances will not change).

    layer_mass_combined = 0
    for layer, layer_name_expected in zip(world, config_to_compare['layers']):
        layer_dict = config_to_compare['layers'][layer_name_expected]

        # Check that layer is in the world_types dict
        assert layer_name_expected in world.__dict__
        assert layer_name_expected.title() in world.__dict__
        assert layer_name_expected in world.layers_by_name
        assert layer_name_expected.title() in world.layers_by_name

        # Check that the layers were built correctly
        assert layer.name == layer_name_expected
        np.testing.assert_approx_equal(layer.radius, layer_dict['radius'])

        layer_mass_combined += layer.mass

    # Check that the sum of layer's mass is the same as the world_types.
    np.testing.assert_approx_equal(world.mass, layer_mass_combined)
    return True

def test_build_layered_from_manual_config():
    """ This will test loading a LayeredWorld from a user-provided configuration dictionary. """

    # Test loading a star from a user-provided configuration dictionary
    io = build_world('Io', io_config)

    # Test that its attributes match expectations
    assert value_check(io, io_config)


def test_build_layered_from_prebuilt_config():
    """ This will test loading a LayeredWorld from a pre-built configuration file. """

    #    Note: the "Io_Simple" config that comes with TidalPy is not designed for BurnMan whereas "Io" is. For this
    #        test we want the former.
    io = build_world('Io_Simple')

    # The pre-built config may not have the same values as the ones used in this file so we will only perform
    #    non-numerical checks.
    assert value_check(io, io.config)


def test_paint_layered():
    """ This will test the slicing features of a LayeredWorld, as well as the painting graphics tool. """

    io = build_world('Io_Simple')
    assert io.paint(auto_show=False)


def test_build_from_layered_world():
    """ This will test building a secondary LayeredWorld from an already built one. """

    # Build a regular layered Io first.
    io = build_world('Io_Simple')

    # Now build a second version that only has a name and eccentricity change.
    new_config = {'name': 'Io_Simple_2', 'eccentricity': 0.5}
    io_2 = build_from_world(io, new_config)

    # io_2 should be almost identical to the original io, so the regular value check should pass
    assert value_check(io_2, io_config, check_name=False)
    assert io.name == 'Io_Simple'
    assert io_2.name == 'Io_Simple_2'
    assert io.config['eccentricity'] == 0.0041
    assert io_2.config['eccentricity'] == 0.5


def test_build_from_layered_world_scale_radius():
    """ This will test building a secondary LayeredWorld from an already built one using a radius scaling technique. """

    # Build a regular layered Io first.
    io = build_world('Io_Simple')

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


def test_build_complex_layered_world():
    """ For this test we will build a world with more than two layers. """
    earth_simple = build_world('earth_simple')

    # Make sure that there are four layers in the world
    assert len(earth_simple.layers) == 4

def test_planet_array_edit():
    """ Make sure that planet.paint does not change the underlying arrays. """
    
    from TidalPy.structures import build_world

    io_simple = build_world('io_simple')

    # Current work around; make a copy of any relevant arrays.
    radii_2 = np.copy(io_simple.radii)

    # Make a call to something that could change the underlying array.
    _ = io_simple.paint(auto_show=False) # Comment/uncomment this line to check the result.
    
    assert np.allclose(io_simple.radii, radii_2)
    
    # While the values in these arrays should be the same; the arrays themselves should not be.
    assert np.allclose(io_simple._radii, io_simple.radii)
    assert not (io_simple._radii is io_simple.radii)
