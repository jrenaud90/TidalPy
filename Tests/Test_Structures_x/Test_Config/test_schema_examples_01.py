"""The schema example files (``Documentation/structures_x/config/examples``) build, solve, round-trip, and between them
use every key the schema accepts, so a key added without an example shows up here."""
import math
import os

import pytest
import toml

import TidalPy.schema_x as schema
from TidalPy.constants import G, mass_jupiter
from TidalPy.Material_x.eos.material_eos import MATERIAL_EOS_CONFIG_KEYS
from TidalPy.cooling_x.cooling import COOLING_CONFIG_KEYS
from TidalPy.partial_melt_x.partial_melt import PARTIAL_MELT_CONFIG_KEYS
from TidalPy.radiogenics_x.radiogenics import RADIOGENICS_CONFIG_KEYS
from TidalPy.rheology_x.rheology import RHEOLOGY_CONFIG_KEYS
from TidalPy.stellar_x.luminosity import LUMINOSITY_CONFIG_KEYS
from TidalPy.structures_x import build_system, build_world, build_world_from_dict
from TidalPy.structures_x.configs.toml_loader import SYSTEM_WORLD_KEYS
from TidalPy.Tides_x.classes.tide import TIDE_CONFIG_KEYS
from TidalPy.viscosity_x.viscosity import VISCOSITY_CONFIG_KEYS

EXAMPLES = os.path.normpath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "Documentation", "structures_x", "config", "examples"))
WORLD_FILES = ("example_world", "example_gasgiant", "example_star", "example_profile_world")
ALL_FILES = WORLD_FILES + ("example_system",)


def _path(name):
    return os.path.join(EXAMPLES, f"{name}.toml")


def _load(name):
    with open(_path(name), encoding="utf-8") as file:
        return toml.load(file)


@pytest.fixture(scope="module")
def examples():
    assert os.path.isdir(EXAMPLES), EXAMPLES
    return {name: _load(name) for name in ALL_FILES}


# =====================================================================================================================
# Every file builds, and the worlds solve
# =====================================================================================================================
@pytest.mark.parametrize("name", WORLD_FILES)
def test_example_world_builds_and_round_trips(name):
    world = build_world(_path(name))
    config = world.get_config_dict()
    rebuilt = build_world_from_dict(config)
    assert rebuilt.get_config_dict() == config
    if world.world_type != "star":
        assert rebuilt.get_solver_defaults() == world.get_solver_defaults()


@pytest.mark.parametrize("name", ("example_world", "example_gasgiant", "example_profile_world"))
def test_example_world_solves(name):
    world = build_world(_path(name))
    result = world.solve_eos(G_to_use=G)
    assert result["success"], result["message"]
    assert world.planet_mass_eos == pytest.approx(world.mass, rel=5.0e-3)
    love = world.solve_love_numbers(frequency=1.5e-5, degree_l=2)
    assert love["success"], love["message"]
    assert 0.0 < world.love_number_k.real < 1.5 * (1.0 + 1.0e-9)
    world.calc_tides(1.5e-5, world.spin_frequency, 0.02, 0.01, 4.0e8, mass_jupiter)
    assert math.isfinite(world.get_tidal_heating()) and world.get_tidal_heating() > 0.0


def test_example_world_pins_every_solver_key():
    world = build_world(_path("example_world"))
    pinned = world.get_solver_defaults()
    assert set(pinned["eos_solver"]) == schema.EOS_SOLVER_KEYS
    assert set(pinned["radial_solver"]) == schema.RADIAL_SOLVER_KEYS


def test_example_star_carries_its_luminosity_model():
    star = build_world(_path("example_star"))
    config = star.get_config_dict()
    assert config["luminosity"]["model"] == "power_law"
    assert config["tides"]["global_tidal_model"] == "fixed_q"


def test_example_system_builds(monkeypatch):
    # The planet's `world` is a relative path, which resolves against the working directory.
    monkeypatch.chdir(EXAMPLES)
    system = build_system(_path("example_system"))
    assert [world.name for world in system] == ["star", "planet", "moon"]
    assert system.get_semi_major_axis("moon") == pytest.approx(4.0e8)


# =====================================================================================================================
# Coverage: every accepted key appears in some example
# =====================================================================================================================
def _layer_tables(examples):
    """Every layer table of every example, including the inline world of the system."""
    for name in WORLD_FILES:
        for layer in (examples[name].get("layers", {}) or {}).values():
            yield layer
    for member in examples["example_system"]["worlds"].values():
        if isinstance(member["world"], dict):
            for layer in member["world"].get("layers", {}).values():
                yield layer


def _model_tables(examples, section):
    for layer in _layer_tables(examples):
        table = layer.get(section)
        if isinstance(table, dict):
            yield table
        if section != "material" and isinstance(layer.get("material"), dict):
            nested = layer["material"].get(section)
            if isinstance(nested, dict):
                yield nested


def _keys(tables):
    found = set()
    for table in tables:
        found |= set(table)
    found.discard("model")
    return found


def test_every_world_level_key_has_an_example(examples):
    layered = {key for name in WORLD_FILES for key in examples[name] if examples[name]["type"] != "star"}
    star_keys = set(examples["example_star"])
    for world_type, allowed in schema.ALLOWED_WORLD_SCALAR_KEYS.items():
        used = star_keys if world_type == "star" else layered
        assert allowed <= used, (world_type, allowed - used)
    assert set(schema.WORLD_MODEL_SECTIONS) <= star_keys
    assert {"data_file", "eos_solver", "radial_solver", "tides"} <= layered


def test_every_tides_key_has_an_example(examples):
    used = _keys(examples[name]["tides"] for name in WORLD_FILES if "tides" in examples[name])
    assert schema.ALLOWED_TIDES_KEYS <= used, schema.ALLOWED_TIDES_KEYS - used


def test_every_layer_key_has_an_example(examples):
    used = _keys(_layer_tables(examples))
    for layer_class, allowed in schema.ALLOWED_LAYER_SCALAR_KEYS.items():
        assert allowed <= used, (layer_class, allowed - used)
    assert set(schema.LAYER_GEOMETRY_SPEC_KEYS) <= used
    assert set(schema.LAYER_MODEL_SECTIONS) <= used
    classes = {layer.get("class") for layer in _layer_tables(examples)}
    assert classes >= {"physics", "solidliquid", "gas"}


@pytest.mark.parametrize("section, accepted", [
    ("material", MATERIAL_EOS_CONFIG_KEYS),
    ("shear_viscosity", VISCOSITY_CONFIG_KEYS),
    ("partial_melt", PARTIAL_MELT_CONFIG_KEYS),
    ("shear_rheology", RHEOLOGY_CONFIG_KEYS),
    ("cooling", COOLING_CONFIG_KEYS),
    ("radiogenics", RADIOGENICS_CONFIG_KEYS),
])
def test_every_model_key_has_an_example(examples, section, accepted):
    used = _keys(_model_tables(examples, section))
    if section == "shear_viscosity":
        used |= _keys(_model_tables(examples, "bulk_viscosity"))
    if section == "shear_rheology":
        used |= _keys(_model_tables(examples, "bulk_rheology"))
    assert accepted <= used, (section, accepted - used)


def test_every_luminosity_tide_and_system_key_has_an_example(examples):
    assert LUMINOSITY_CONFIG_KEYS - {"luminosity_w"} <= set(examples["example_star"]["luminosity"])
    tides = _keys(examples[name]["tides"] for name in WORLD_FILES if "tides" in examples[name])
    assert TIDE_CONFIG_KEYS <= tides
    members = _keys(examples["example_system"]["worlds"].values())
    assert set(SYSTEM_WORLD_KEYS) <= members, set(SYSTEM_WORLD_KEYS) - members


def test_every_model_of_every_family_is_named_somewhere(examples):
    """Each family's models appear across the layer tables (aliases aside)."""
    named = {}
    for section in ("material", "shear_viscosity", "bulk_viscosity", "partial_melt", "shear_rheology", "bulk_rheology",
                    "cooling", "radiogenics"):
        named[section] = {table["model"] for table in _model_tables(examples, section) if "model" in table}
    assert named["material"] >= {"constant", "bm", "vinet", "interpolate"}
    assert named["shear_viscosity"] | named["bulk_viscosity"] >= {"constant", "reference", "arrhenius"}
    assert named["partial_melt"] >= {"off", "spohn", "henning"}
    rheologies = named["shear_rheology"] | named["bulk_rheology"]
    assert rheologies >= {"elastic", "maxwell", "andrade", "sundberg", "voigt", "burgers"}
    assert named["cooling"] >= {"off", "convection", "conduction"}
    assert named["radiogenics"] >= {"off", "isotope", "fixed"}
