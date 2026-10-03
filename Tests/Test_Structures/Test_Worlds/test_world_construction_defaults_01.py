"""A world built by hand takes the same defaults as one the builder makes, a layer without a temperature is named
once, and the tide, luminosity, and radiogenics models a world or layer is given are copied in."""
import math

import pytest

import TidalPy
from TidalPy.Radiogenics.radiogenics import make_radiogenics
from TidalPy.Stellar.luminosity import make_luminosity
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.worlds.gasgiant import GasGiantWorld
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Structures.worlds.terrestrial import TerrestrialWorld
from TidalPy.Tides.classes.tide import make_tide


_SUN_RADIUS = 6.957e8   # [m]
_SUN_MASS = 1.988e30    # [kg]
_COLD_WARNING = "has no temperature of its own"


def _scalars(world):
    config = world.get_config_dict()
    scalars = {key: config[key] for key in ("albedo", "emissivity", "obliquity_rad", "spin_frequency_rad_s")}
    # The configured moment-of-inertia factor is not written back, so compare the spin model's value.
    scalars["moment_of_inertia"] = world.get_moment_of_inertia()
    return scalars


# =====================================================================================================================
# Manual worlds take the [worlds] defaults
# =====================================================================================================================
@pytest.mark.parametrize("world_class, world_type", [
    (BaseWorld, "layered"), (TerrestrialWorld, "terrestrial"), (GasGiantWorld, "gasgiant"), (StarWorld, "star")])
def test_a_manual_world_matches_a_built_one(world_class, world_type):
    manual = world_class("manual", _SUN_RADIUS, _SUN_MASS, world_type=world_type)
    config = {"schema_version": "0.2.0", "name": "built", "type": world_type, "radius_m": _SUN_RADIUS,
              "mass_kg": _SUN_MASS}
    if world_type != "star":
        config["layers"] = {"body": {"radius_fraction": 1.0, "material": "simple_rock"}}
    assert _scalars(manual) == _scalars(build_world(config))


def test_a_manual_star_takes_the_star_defaults():
    star = StarWorld("sun", _SUN_RADIUS, _SUN_MASS)
    defaults = TidalPy.config["worlds"]["star"]
    # The [worlds.star] default is the configuration's to supply, so the star's config leaves it out.
    assert "moment_of_inertia_factor" not in star.get_config_dict()
    assert star.albedo == defaults["albedo"]
    assert star.effective_temperature == defaults["effective_temperature_k"]
    assert star.get_moment_of_inertia() == pytest.approx(
        defaults["moment_of_inertia_factor"] * _SUN_MASS * _SUN_RADIUS**2, rel=1e-14)


def test_given_values_win():
    world = TerrestrialWorld("planet", 6.4e6, 6.0e24, albedo=0.1, emissivity=0.9, obliquity=0.2, spin_frequency=1e-5)
    assert (world.albedo, world.emissivity, world.obliquity, world.spin_frequency) == (0.1, 0.9, 0.2, 1e-5)


# =====================================================================================================================
# A layer without a temperature
# =====================================================================================================================
def _rock_world(temperature=None):
    radius = 1.0e6
    kwargs = {} if temperature is None else {"temperature": temperature}
    layer = Layer("rock", 0, 0.0, radius, material="peridotite", **kwargs)
    world = BaseWorld("rock", radius, 4.0 / 3.0 * math.pi * 3300.0 * radius**3)
    world.add_layer(layer)
    return world


def test_a_rigid_layer_without_a_temperature_is_named_once(spdlog_text):
    world = _rock_world()
    world.solve_eos()
    world.solve_eos()
    text = spdlog_text()
    assert text.count(_COLD_WARNING) == 1
    assert "'rock'" in text
    # Setting the temperature again (here, still to 0 K) re-arms it.
    world.rock.temperature = 0.0
    world.solve_eos()
    assert spdlog_text().count(_COLD_WARNING) == 2


def test_a_layer_with_a_temperature_is_not_named(spdlog_text):
    _rock_world(1600.0).solve_eos()
    assert _COLD_WARNING not in spdlog_text()


# =====================================================================================================================
# Models are copied in
# =====================================================================================================================
def test_one_tide_model_serves_several_worlds():
    tide = make_tide("fixed_q", {"fixed_k": [0.3], "fixed_q": [50.0]})
    first = BaseWorld("first", 1.0e6, 1.0e22)
    second = BaseWorld("second", 2.0e6, 2.0e22)
    first.set_tide_model(tide)
    second.set_tide_model(tide)
    assert first.tide_model_set and second.tide_model_set
    assert first.get_config_dict()["tides"] == second.get_config_dict()["tides"]
    assert tide.get_config_dict()["fixed_q"][0] == 50.0


def test_one_luminosity_model_serves_several_stars():
    model = make_luminosity("power_law", {"power_law_coeff": 1.2, "power_law_exponent": 3.9})
    stars = [StarWorld(name, _SUN_RADIUS, _SUN_MASS) for name in ("a", "b")]
    for star in stars:
        star.set_luminosity_model(model)
    assert stars[0].calc_luminosity_from_mass() == stars[1].calc_luminosity_from_mass()
    assert model.get_config_dict()["power_law_exponent"] == 3.9


def test_one_radiogenics_model_serves_several_layers():
    model = make_radiogenics("fixed", {"fixed_heat_production_w_kg": 1.0e-11})
    layers = [Layer(name, 0, 0.0, 1.0e6, material="simple_rock") for name in ("a", "b")]
    for layer in layers:
        layer.radiogenics = model
    assert layers[0].calc_radiogenic_heating(0.0, 1.0e20) == layers[1].calc_radiogenic_heating(0.0, 1.0e20)
    assert layers[0].calc_radiogenic_heating(0.0, 1.0e20) == pytest.approx(1.0e9)
    assert model.get_config_dict()["fixed_heat_production_w_kg"] == 1.0e-11
