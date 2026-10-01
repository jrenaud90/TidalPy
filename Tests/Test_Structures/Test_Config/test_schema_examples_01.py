"""The schema example files build, solve, round-trip, and between them use every key the schema accepts."""
import math
import os

import pytest
import toml

import TidalPy.schema as schema
from TidalPy.constants import G, mass_jupiter
from TidalPy.Cooling.cooling import COOLING_CONFIG_KEYS, make_cooling
from TidalPy.Material import available_materials
from TidalPy.Radiogenics.radiogenics import RADIOGENICS_CONFIG_KEYS, make_radiogenics
from TidalPy.Rheology.rheology import RHEOLOGY_CONFIG_KEYS
from TidalPy.Stellar.luminosity import LUMINOSITY_CONFIG_KEYS
from TidalPy.Structures import build_system, build_world, build_world_from_dict
from TidalPy.Structures.configs.toml_loader import SYSTEM_WORLD_KEYS
from TidalPy.Tides.classes.tide import TIDE_CONFIG_KEYS
from TidalPy.Utilities.classes.families import get_family

EXAMPLES = os.path.normpath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "Documentation", "Structures", "config", "examples"))
WORLD_FILES = ("example_world", "example_gasgiant", "example_star", "example_profile_world",
               "example_profile_q_world")
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


@pytest.mark.parametrize(
    "name", ("example_world", "example_gasgiant", "example_profile_world", "example_profile_q_world"))
def test_example_world_solves(name):
    world = build_world(_path(name))
    result = world.solve_eos(G_to_use=G)
    assert result["success"], result["message"]
    assert world.planet_mass_eos == pytest.approx(world.mass, rel=5.0e-3)
    love = world.solve_love_numbers(frequency=1.5e-5, degree_l=2)
    assert love["success"], love["message"]
    assert 0.0 < world.love_number_k.real < 1.5 * (1.0 + 1.0e-9)
    world.calc_tides(
        1.5e-5,
        world.spin_frequency,
        0.02,
        0.01,
        4.0e8,
        mass_jupiter,
    )
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
    """The tables of a layer-held model family: the layer's own, and its material phases' for a rheology."""
    for layer in _layer_tables(examples):
        table = layer.get(section)
        if isinstance(table, dict):
            yield table
        material = layer.get("material")
        if isinstance(material, dict):
            for slot in ("solid", "liquid"):
                nested = (material.get(slot) or {}).get(section)
                if isinstance(nested, dict):
                    yield nested


def _material_tables(examples):
    """Every layer's material, as written: a MatPack name or a table."""
    for layer in _layer_tables(examples):
        if "material" in layer:
            yield layer["material"]


# The law slots of a material table, each with the model family it holds: a phase's laws, and the melting laws.
_PHASE_LAW_FAMILIES = {
    "eos":             "equation of state",
    "shear_modulus":   "shear modulus",
    "shear_viscosity": "viscosity",
    "bulk_viscosity":  "viscosity",
    "shear_rheology":  "rheology",
    "bulk_rheology":   "rheology",
}
_MELTING_LAW_FAMILIES = {
    "solidus":               "melting curve",
    "liquidus":              "melting curve",
    "weakening":             "melt weakening",
    "bulk_modulus_mixing":   "bulk-modulus mixing",
    "bulk_viscosity_mixing": "bulk-viscosity mixing",
}


def _law_tables(examples, family):
    """The law tables of one family in the examples' material tables (a preset's own laws are not written there)."""
    for material in _material_tables(examples):
        if not isinstance(material, dict):
            continue
        for slot in ("solid", "liquid"):
            phase = material.get(slot) or {}
            for law, law_family in _PHASE_LAW_FAMILIES.items():
                table = phase.get(law)
                if law_family == family and isinstance(table, dict):
                    yield table
                    # A composite viscosity's mechanisms are viscosity laws too.
                    for mechanism in table.get("mechanisms", ()):
                        yield mechanism
        melting = material.get("melting") or {}
        for law, law_family in _MELTING_LAW_FAMILIES.items():
            if law_family == family and isinstance(melting.get(law), dict):
                yield melting[law]


def _phase_tables(examples):
    for material in _material_tables(examples):
        if isinstance(material, dict):
            for slot in ("solid", "liquid"):
                if isinstance(material.get(slot), dict):
                    yield material[slot]


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
    assert schema.LAYER_SCALAR_KEYS <= used, schema.LAYER_SCALAR_KEYS - used
    assert set(schema.LAYER_GEOMETRY_SPEC_KEYS) <= used
    assert set(schema.LAYER_MODEL_SECTIONS) <= used
    assert "layer_index" in used


def test_every_material_form_has_an_example(examples):
    """A material is given as a MatPack name, as a preset with overrides, and as a full table."""
    materials = list(_material_tables(examples))
    names = [material for material in materials if isinstance(material, str)]
    presets = [material for material in materials if isinstance(material, dict) and "preset" in material]
    full = [material for material in materials if isinstance(material, dict) and "preset" not in material
            and ("solid" in material or "liquid" in material)]
    assert names and presets and full
    assert set(names) | {material["preset"] for material in presets} <= set(available_materials())
    # A full table of each kind: one with a solid phase and one with only a liquid phase.
    assert any("solid" in material for material in full)
    assert any("liquid" in material and "solid" not in material for material in full)


def test_every_material_and_phase_key_has_an_example(examples):
    material_keys = _keys(material for material in _material_tables(examples) if isinstance(material, dict))
    assert get_family("material").config_keys_of("material") <= material_keys
    phase_keys = _keys(_phase_tables(examples))
    assert get_family("phase").config_keys_of("phase") <= phase_keys


@pytest.mark.parametrize("sections, accepted", [
    (("shear_rheology", "bulk_rheology"), RHEOLOGY_CONFIG_KEYS),
    (("cooling",), COOLING_CONFIG_KEYS),
    (("radiogenics",), RADIOGENICS_CONFIG_KEYS),
], ids=lambda value: value[0] if isinstance(value, tuple) else None)
def test_every_model_key_has_an_example(examples, sections, accepted):
    used = set()
    for section in sections:
        used |= _keys(_model_tables(examples, section))
    assert accepted <= used, (sections, accepted - used)


def test_every_luminosity_tide_and_system_key_has_an_example(examples):
    assert LUMINOSITY_CONFIG_KEYS - {"luminosity_w"} <= set(examples["example_star"]["luminosity"])
    tides = _keys(examples[name]["tides"] for name in WORLD_FILES if "tides" in examples[name])
    assert TIDE_CONFIG_KEYS <= tides
    members = _keys(examples["example_system"]["worlds"].values())
    assert set(SYSTEM_WORLD_KEYS) <= members, set(SYSTEM_WORLD_KEYS) - members


def _canonical(family_name, model_name):
    # Cooling and radiogenics have no Python name lookup; a model built by name reports its canonical name.
    if family_name == "cooling":
        return make_cooling(model_name, {}).model_name
    if family_name == "radiogenics":
        return make_radiogenics(model_name, {}).model_name
    return get_family(family_name).canonical_name(model_name)


@pytest.mark.parametrize("sections, family_name, expected_models", [
    pytest.param(("shear_rheology", "bulk_rheology"), "rheology", set(get_family("rheology").model_names()),
                 id="rheology"),
    pytest.param(("cooling",), "cooling", {"off", "convection", "conduction"}, id="cooling"),
    pytest.param(("radiogenics",), "radiogenics", {"off", "isotope", "fixed"}, id="radiogenics"),
])
def test_every_model_of_every_layer_family_is_named_somewhere(examples, sections, family_name, expected_models):
    """Each layer-held family's models appear across the layer tables (aliases aside)."""
    named = set()
    for section in sections:
        named |= {_canonical(family_name, table["model"])
                  for table in _model_tables(examples, section) if "model" in table}
    expected = {_canonical(family_name, model) for model in expected_models}
    assert named >= expected, expected - named


# The models the examples' material tables name, by law family. The examples do not name every model of these
# families; each model they name is checked to show every key it reads.
_MATERIAL_LAW_MODELS = {
    "equation of state":     {"constant", "birch_murnaghan", "vinet", "interpolate"},
    "shear modulus":         {"linear", "interpolate"},
    "viscosity":             {"arrhenius", "reference", "constant", "interpolate", "composite"},
    "melting curve":         {"simon_glatzel_2"},
    "melt weakening":        {"spohn"},
    "bulk-modulus mixing":   {"hashin_shtrikman"},
    "bulk-viscosity mixing": {"compaction"},
}


@pytest.mark.parametrize("family_name", sorted(_MATERIAL_LAW_MODELS))
def test_every_material_law_model_and_its_keys_have_an_example(examples, family_name):
    """The material tables name these models of each law family, and between them show every key those models read."""
    family = get_family(family_name)
    tables = list(_law_tables(examples, family_name))
    named = {family.canonical_name(table["model"]) for table in tables}
    assert named >= _MATERIAL_LAW_MODELS[family_name], _MATERIAL_LAW_MODELS[family_name] - named
    needed = set().union(*(family.config_keys_of(model) for model in named))
    used = _keys(tables)
    assert needed <= used, needed - used
