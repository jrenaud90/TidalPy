import numpy as np

import TidalPy

io_dict = {
    "name": "Io",
    "type": "layered",
    "radius": 1821.49e3,
    "orbital_period": 1.769,
    "eccentricity": 0.0041,
    "spin_period": 1.769,
    "albedo": 0.63,
    "force_spin_sync": True,
    "layers": {
        "Core": {
            "type": "iron",
            "is_tidal": False,
            "radius": 810.0e3,
            "density": 5200.
        },
        "Mantle": {
            "type": "rock",
            "is_tidal": True,
            "radius": 1821.49e3,
            "surface_temperature": 100.0,
            "density": 3200.
        }
    }
}

def value_check(world, config_to_compare: dict):

    assert world.name == config_to_compare['name']
    assert world.num_layers == 2
    assert len(world.layers) == 2
    # Note: <world>.layers_by_name will usually contain more than the expected number of layers due to allowing both
    #    title and lowercase keys (number of layer instances will not change).

    layer_mass_combined = 0
    for layer, layer_name_expected in zip(world, config_to_compare['layers']):
        layer_dict = config_to_compare['layers'][layer_name_expected]

        # Check that layer is in the worlds dict
        assert layer_name_expected in world.__dict__
        assert layer_name_expected.title() in world.__dict__
        assert layer_name_expected in world.layers_by_name
        assert layer_name_expected.title() in world.layers_by_name

        # Check that the layers were built correctly
        assert layer.name == layer_name_expected
        np.testing.assert_approx_equal(layer.radius, layer_dict['radius'])

        layer_mass_combined += layer.mass

    # Check that the sum of layer's mass is the same as the worlds.
    np.testing.assert_approx_equal(world.mass, layer_mass_combined)
    return True


def test_build_layered_from_manual_config():

    # Test loading a star from a user-provided configuration dictionary
    io = TidalPy.build_world('Io', io_dict)

    # Test that its attributes match expectations
    assert value_check(io, io_dict)


def test_build_layered_from_prebuilt_config():
    # Test loading a star from a pre-built configuration file
    #    Note: the "Io_Simple" config that comes with TidalPy is not designed for BurnMan whereas "Io" is. For this
    #        test we want the former.
    io = TidalPy.build_world('Io_Simple')

    # The pre-built config may not have the same values as the ones used in this file so we will only perform
    #    non-numerical checks.
    assert value_check(io, io.config)


def test_paint_layered():
    # This will test the slicing features of the worlds, as well as the painting graphics tool
    io = TidalPy.build_world('Io_Simple')
    assert io.paint(auto_show=False)