"""``get_config_dict`` of a layer, world, or system rebuilds the same object through ``build_*_from_dict``."""
import copy
import math

import pytest

from TidalPy.Dynamics.spin import Spin
from TidalPy.Material import load_material
from TidalPy.Stellar.luminosity import make_luminosity
from TidalPy.Structures import (
    available_worlds,
    build_layer_from_dict,
    build_system,
    build_system_from_dict,
    build_world,
    build_world_from_dict,
)
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.system.system import System


# =====================================================================================================================
# Layers
# =====================================================================================================================
def _base_layer():
    return Layer(
        "shell",
        1,
        1.0e6,
        2.0e6,
        3.0e22,
        material=load_material("simple_rock"),
        use_tides=False,
        tidal_scale=0.25,
    )


def _layer_with_state():
    return Layer(
        "mantle",
        2,
        1.0e6,
        2.0e6,
        3.0e22,
        material=load_material("peridotite"),
        state="liquid",
        is_volume_fixed=False,
        is_static=False,
        is_incompressible=True,
        temperature=1500.0,
        use_thermal_expansion=True,
        use_melting=True,
        use_pressure_melting=True,
        use_melt_density=True,
        use_heating=True,
    )


def _world_layer(world_name, layer_name):
    """A layer of a bundled world, with every model the builder attaches."""
    world = build_world(world_name)
    for layer in world:
        if layer.name == layer_name:
            return layer
    raise KeyError(layer_name)


@pytest.mark.parametrize("make_layer", [
    pytest.param(_base_layer, id="base"),
    pytest.param(_layer_with_state, id="base-with-state"),
    pytest.param(lambda: _world_layer("io", "mantle"), id="io-mantle"),
    pytest.param(lambda: _world_layer("io", "core"), id="io-core"),
    pytest.param(lambda: _world_layer("jupiter_simple", "envelope"), id="jupiter-envelope"),
])
def test_layer_rebuilds_from_its_config_dict(make_layer):
    layer = make_layer()
    config = layer.get_config_dict()
    snapshot = copy.deepcopy(config)

    rebuilt = build_layer_from_dict(config)
    assert config == snapshot, "the input dictionary must not be modified"
    assert type(rebuilt) is type(layer) is Layer
    assert rebuilt.get_config_dict() == config
    assert rebuilt.is_liquid == layer.is_liquid


@pytest.mark.parametrize("missing", ["name", "radius_inner_m", "radius_outer_m"])
def test_layer_config_needs_its_standalone_keys(missing):
    config = _layer_with_state().get_config_dict()
    del config[missing]
    with pytest.raises(ValueError, match=missing):
        build_layer_from_dict(config)


@pytest.mark.parametrize("builder, source, match", [
    pytest.param(build_layer_from_dict, "mantle", None, id="layer"),
    pytest.param(build_world_from_dict, "io", "build_world", id="world"),
    pytest.param(build_system_from_dict, "sol_system", "build_system", id="system"),
])
def test_builder_rejects_a_non_dict(builder, source, match):
    with pytest.raises(TypeError, match=match):
        builder(source)


# =====================================================================================================================
# Worlds
# =====================================================================================================================
@pytest.mark.parametrize("world_name", available_worlds())
def test_world_rebuilds_from_its_config_dict(world_name):
    world = build_world(world_name)
    config = world.get_config_dict()
    snapshot = copy.deepcopy(config)

    rebuilt = build_world_from_dict(config)
    assert config == snapshot, "the input dictionary must not be modified"
    assert type(rebuilt) is type(world)
    assert rebuilt.get_config_dict() == config
    assert rebuilt.config is not config


def test_world_changes_made_after_the_build_are_in_the_dict():
    """The dictionary is the world as it stands, not the file it came from."""
    world = build_world("io")
    world.set_spin_frequency(1.234e-5)
    world.set_obliquity(0.05)
    world.set_spin_model(Spin(moment_of_inertia_factor=0.3769))
    for layer in world:
        if layer.name == "mantle":
            layer.temperature = 1777.0
            layer.is_incompressible = True

    rebuilt = build_world_from_dict(world.get_config_dict())
    assert rebuilt.spin_frequency == 1.234e-5
    assert rebuilt.obliquity == 0.05
    assert rebuilt.get_config_dict()["moment_of_inertia_factor"] == 0.3769
    # Before an EOS solve the moment of inertia is the spin model's estimate.
    assert rebuilt.get_moment_of_inertia() == pytest.approx(0.3769 * rebuilt.mass * rebuilt.radius**2)
    mantle = [layer for layer in rebuilt if layer.name == "mantle"][0]
    assert mantle.temperature == 1777.0
    assert mantle.is_incompressible is True


def test_rebuilt_world_solves_to_the_same_numbers():
    world = build_world("io")
    rebuilt = build_world_from_dict(world.get_config_dict())
    for each in (world, rebuilt):
        each.solve_eos()
        each.solve_love_numbers(frequency=4.11e-5)
    assert rebuilt.planet_mass_eos == world.planet_mass_eos
    assert rebuilt.love_number_k == world.love_number_k


def test_star_luminosity_model_survives_the_round_trip():
    star = build_world("sol")
    star.set_luminosity_model(make_luminosity("power_law", {"power_law_coeff": 1.2, "power_law_exponent": 3.8}))
    config = star.get_config_dict()
    assert config["luminosity"]["model"] == "power_law"

    rebuilt = build_world_from_dict(config)
    assert rebuilt.luminosity_model_set
    assert rebuilt.get_config_dict() == config
    assert rebuilt.calc_luminosity_from_mass() == star.calc_luminosity_from_mass()
    # Attaching the model does not move the stored luminosity.
    assert rebuilt.luminosity == star.luminosity


@pytest.mark.parametrize("key, value", [
    pytest.param("luminosity", {"model": "fixed"}, id="luminosity-on-a-planet"),
    pytest.param("moment_of_inertia_factor", 0.9, id="moi-factor-out-of-range"),
])
def test_world_builder_rejects_a_bad_key(key, value):
    config = build_world("io").get_config_dict()
    config[key] = value
    with pytest.raises(ValueError, match=key):
        build_world_from_dict(config)


# =====================================================================================================================
# Systems
# =====================================================================================================================
def test_bundled_system_rebuilds_from_its_config_dict():
    system = build_system("sol_system")
    config = system.get_config_dict()
    snapshot = copy.deepcopy(config)

    rebuilt = build_system_from_dict(config)
    assert config == snapshot, "the input dictionary must not be modified"
    assert isinstance(rebuilt, System)
    assert rebuilt.get_config_dict() == config
    assert [world.name for world in rebuilt] == [world.name for world in system]
    assert rebuilt.get_tidal_host("earth").name == "sun"
    assert rebuilt.star.name == "sun"


def test_system_dict_carries_the_live_state_of_its_worlds():
    system = build_system("sol_system")
    system["earth"].set_spin_frequency(9.9e-5)
    system.set_eccentricity("earth", 0.05)

    rebuilt = build_system_from_dict(system.get_config_dict())
    assert rebuilt["earth"].spin_frequency == 9.9e-5
    assert rebuilt.get_eccentricity("earth") == 0.05
    assert math.isclose(rebuilt.get_semi_major_axis("earth"), system.get_semi_major_axis("earth"))


def test_mutual_pair_survives_the_round_trip():
    system = System("pair")
    earth = build_world("earth_simple")
    earth.name = "earth"
    moon = build_world("luna")
    moon.name = "moon"
    system.add_world(earth)
    system.add_world(moon, tidal_host=earth, semi_major_axis=3.844e8, eccentricity=0.0549)
    system.set_tidal_host(earth, moon)

    rebuilt = build_system_from_dict(system.get_config_dict())
    assert rebuilt.is_mutual_pair("earth") and rebuilt.is_mutual_pair("moon")
    assert rebuilt.get_semi_major_axis("earth") == 3.844e8
    assert rebuilt.get_config_dict() == system.get_config_dict()
