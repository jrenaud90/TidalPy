"""The world builder refuses numbers that do not describe a possible world, and names the one at fault."""
import copy
import math

import pytest

from TidalPy.Structures import available_worlds, build_world
from TidalPy.Structures.configs import load_toml, resolve_world_path, validate_world_config
from TidalPy.Structures.configs.toml_loader import validate_physical_values


_DELETE = object()


def _config():
    return {
        "name": "checked", "type": "terrestrial", "radius_m": 2.0e6, "mass_kg": 1.0e23,
        "layers": {
            "core":   {"class": "base", "type": "iron", "radius_fraction": 0.5},
            "mantle": {"class": "base", "type": "mantle_rock", "radius_fraction": 1.0},
        },
    }


def _edited(*edits):
    """The sound config with each ``(path, value)`` applied; ``_DELETE`` removes the key."""
    config = _config()
    for path, value in edits:
        target = config
        for key in path[:-1]:
            target = target[key]
        if value is _DELETE:
            del target[path[-1]]
        else:
            target[path[-1]] = value
    return config


def test_a_sound_world_passes():
    validate_world_config(_config())
    build_world(_config())


@pytest.mark.parametrize("world_name", available_worlds())
def test_every_bundled_world_passes(world_name):
    config = load_toml(resolve_world_path(world_name))
    if "data_file" in config or "data" in config:
        # A profile world is expanded into layers before it is validated, so the build covers it.
        build_world(world_name)
    else:
        validate_physical_values(config)


@pytest.mark.parametrize("path, value, fragment", [
    (("radius_m",), -2.0e6, "radius_m"),
    (("radius_m",), 0.0, "radius_m"),
    (("radius_m",), math.nan, "finite"),
    (("radius_m",), "2000 km", "must be a number"),
    (("mass_kg",), 0.0, "mass_kg"),
    (("mass_kg",), math.inf, "finite"),
    (("albedo",), 1.2, "albedo"),
    (("albedo",), -0.1, "albedo"),
    (("emissivity",), 0.0, "emissivity"),
    (("spin_frequency_rad_s",), math.nan, "spin_frequency_rad_s"),
    (("layers", "core", "radius_fraction"), 0.0, "radius_fraction"),
    (("layers", "core", "radius_fraction"), -0.5, "radius_fraction"),
    (("layers", "mantle", "radius_fraction"), 1.5, "radius_fraction"),
    (("layers", "mantle", "mass_kg"), -1.0, "mass_kg"),
    (("layers", "mantle", "tidal_scale"), math.nan, "tidal_scale"),
    (("layers", "mantle", "temperature_k"), -5.0, "temperature_k"),
    (("layers", "core", "layer_index"), -1, "layer_index"),
    (("layers", "core", "layer_index"), 0.5, "layer_index"),
])
def test_an_impossible_value_is_named(path, value, fragment):
    with pytest.raises(ValueError, match=fragment):
        validate_world_config(_edited((path, value)))
    with pytest.raises(ValueError, match=fragment):
        build_world(_edited((path, value)))


@pytest.mark.parametrize("edits, fragment", [
    pytest.param(
        [(("layers", "core", "radius_fraction"), 0.6), (("layers", "mantle", "radius_fraction"), 0.6)],
        "not above where it starts",
        id="zero-thickness"),
    pytest.param([(("layers", "mantle", "radius_fraction"), 0.9)], "short of the world's radius", id="short"),
    pytest.param(
        [(("layers", "mantle", "radius_fraction"), _DELETE), (("layers", "mantle", "radius_outer_m"), 2.5e6)],
        "above the world's radius",
        id="past-the-surface"),
])
def test_the_layers_must_stack_to_the_surface(edits, fragment):
    with pytest.raises(ValueError, match=fragment):
        validate_world_config(_edited(*edits))


def test_volume_fractions_are_stacked_like_the_builder_stacks_them():
    config = _config()
    for name, fraction in (("core", 0.2), ("mantle", 0.8)):
        del config["layers"][name]["radius_fraction"]
        config["layers"][name]["volume_fraction"] = fraction
    validate_world_config(config)
    config["layers"]["mantle"]["volume_fraction"] = 0.5
    with pytest.raises(ValueError, match="short of the world's radius"):
        validate_world_config(config)


def test_two_layers_cannot_share_an_index():
    config = _config()
    # The mantle takes index 1 from its place in the table.
    config["layers"]["core"]["layer_index"] = 1
    with pytest.raises(ValueError, match="both resolve to layer_index 1"):
        validate_world_config(config)
    # Swapped outright is a valid, if odd, way to order them.
    config["layers"]["mantle"]["layer_index"] = 0
    config["layers"]["core"]["radius_fraction"] = 1.0
    config["layers"]["mantle"]["radius_fraction"] = 0.5
    validate_world_config(config)


def test_a_star_is_checked_too():
    star = {"name": "s", "type": "star", "radius_m": 7.0e8, "mass_kg": 2.0e30, "luminosity_w": -1.0}
    with pytest.raises(ValueError, match="luminosity_w"):
        validate_world_config(star)
    # Zero asks for the luminosity to be derived from the temperature.
    star["luminosity_w"] = 0.0
    validate_world_config(star)


def test_the_input_is_left_alone():
    config = _config()
    snapshot = copy.deepcopy(config)
    validate_world_config(config)
    assert config == snapshot
