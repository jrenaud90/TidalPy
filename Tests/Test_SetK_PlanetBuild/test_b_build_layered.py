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

def value_check(world, check_numerical: bool = True):

    assert world.name == io_dict['name']
    assert len(world.layers) == 2
    assert len(world.layers_by_name) == 2

    layer_mass_combined = 0
    for layer, layer_name_expected in zip(world, io_dict['layers']):
        layer_dict = io_dict['layers'][layer_name_expected]

        # Check that layer is in the worlds dict
        assert layer_name_expected in world.__dict__

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
    assert value_check(io)