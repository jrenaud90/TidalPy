"""Builder conveniences: binary files through load_world, load_system, and the builders; overrides; the schema-version
rule for dicts; closest-name suggestions; and the TidalPy.Structures exports."""
import copy
import pathlib
import warnings

import pytest

import TidalPy.Structures as structures
from TidalPy.Structures import (
    BaseWorld,
    GasGiantWorld,
    Layer,
    StarWorld,
    System,
    TerrestrialWorld,
    available_worlds,
    build_system,
    build_world,
    load_system,
    load_world,
    make_tide,
)
from TidalPy.Utilities.binary import binary_file_class


# =====================================================================================================================
# Binary files
# =====================================================================================================================
@pytest.mark.parametrize("name, world_class", [
    ("io", TerrestrialWorld),
    ("jupiter_simple", GasGiantWorld),
    ("sol", StarWorld),
])
def test_load_world_needs_no_placeholder(tmp_path, name, world_class):
    world = build_world(name)
    path = pathlib.Path(tmp_path) / f"{name}.tpyb"
    world.save_binary(path)
    assert binary_file_class(path) == world_class.__name__
    for loaded in (load_world(path), load_world(str(path)), build_world(path)):
        assert type(loaded) is world_class
        assert loaded.get_config_dict() == world.get_config_dict()


def test_load_system_needs_no_placeholder(tmp_path):
    system = build_system("sol_system")
    path = pathlib.Path(tmp_path) / "sol_system.tpyb"
    system.save_binary(path)
    assert binary_file_class(path) == "System"
    for loaded in (load_system(path), build_system(path)):
        assert isinstance(loaded, System)
        assert loaded.get_config_dict() == system.get_config_dict()


def test_loaders_name_what_a_file_holds(tmp_path):
    world_path = pathlib.Path(tmp_path) / "io.tpyb"
    build_world("io").save_binary(world_path)
    system_path = pathlib.Path(tmp_path) / "sol_system.tpyb"
    build_system("sol_system").save_binary(system_path)
    toml_path = pathlib.Path(tmp_path) / "io.toml"
    build_world("io").save_to_toml(toml_path)
    assert binary_file_class(toml_path) is None
    with pytest.raises(IOError, match="a System file; load it with load_system"):
        load_world(system_path)
    with pytest.raises(IOError, match="a TerrestrialWorld file; load it with load_world"):
        load_system(world_path)
    with pytest.raises(IOError, match="not a TidalPy binary file"):
        load_world(toml_path)
    with pytest.raises(FileNotFoundError):
        load_world(pathlib.Path(tmp_path) / "missing.tpyb")
    with pytest.raises(ValueError, match="binary file"):
        build_world(world_path, overrides={"name": "other"})


def test_a_corrupt_world_file_names_the_file(tmp_path):
    path = pathlib.Path(tmp_path) / "io.tpyb"
    build_world("io").save_binary(path)
    path.write_bytes(path.read_bytes()[:-5])
    with pytest.raises(IOError, match="io.tpyb"):
        load_world(path)


# =====================================================================================================================
# Overrides
# =====================================================================================================================
def test_build_world_overrides_merge_table_by_table():
    overrides = {"tides": {"eccentricity_trunc_lvl": 20}, "layers": {"mantle": {"temperature_k": 1700.0}}}
    snapshot = copy.deepcopy(overrides)
    world = build_world("io", overrides)
    reference = build_world("io")
    assert overrides == snapshot
    assert world.get_tide_config()["eccentricity_trunc_lvl"] == 20
    assert world["mantle"].temperature == 1700.0
    # Everything the overrides leave out is as the file has it.
    assert world["mantle"].radius_outer == reference["mantle"].radius_outer
    assert world.get_tide_config()["max_degree_l"] == reference.get_tide_config()["max_degree_l"]


def test_build_system_overrides():
    # The bundled file gives the Earth's one orbit about the Sun in both element sets, so both change.
    system = build_system(
        "sol_system", overrides={"worlds": {"earth": {"eccentricity": 0.05, "stellar_eccentricity": 0.05}}})
    assert system.get_eccentricity("earth") == 0.05
    assert system.get_eccentricity("jupiter") == 0.0489


def test_overrides_must_be_a_dict():
    with pytest.raises(TypeError, match="overrides"):
        build_world("io", overrides=[("name", "x")])


# =====================================================================================================================
# Schema version
# =====================================================================================================================
def test_a_dict_without_a_schema_version_builds_silently():
    config = build_world("io").get_config_dict()
    del config["schema_version"]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        world = build_world(config)
    assert world.source_config["schema_version"] == structures.SCHEMA_VERSION


def test_a_file_without_a_schema_version_still_warns(tmp_path):
    path = pathlib.Path(tmp_path) / "io.toml"
    build_world("io").save_to_toml(path)
    text = path.read_text(encoding="utf-8")
    path.write_text("\n".join(line for line in text.splitlines() if not line.startswith("schema_version")),
                    encoding="utf-8")
    with pytest.warns(UserWarning, match="schema_version"):
        build_world(path)


# =====================================================================================================================
# Closest-name suggestions
# =====================================================================================================================
def test_an_unknown_bundled_world_names_the_closest():
    with pytest.raises(FileNotFoundError, match="did you mean 'io'") as error:
        build_world("ioo")
    assert "earth_simple" in str(error.value)


def test_an_unknown_bundled_system_names_the_closest():
    with pytest.raises(FileNotFoundError, match="did you mean 'sol_system'") as error:
        build_system("sol_sistem")
    assert "Bundled systems: sol_system" in str(error.value)


@pytest.mark.parametrize("edit, match", [
    pytest.param(lambda config: config["layers"]["mantle"].update(temprature_k=1600.0),
                 "did you mean 'temperature_k'", id="layer-key"),
    pytest.param(lambda config: config["tides"].update(max_degre_l=3), "did you mean 'max_degree_l'", id="tides-key"),
    pytest.param(lambda config: config.update(albedoo=0.3), "did you mean 'albedo'", id="world-key"),
])
def test_an_unknown_key_names_the_closest(edit, match):
    config = build_world("io").get_config_dict()
    edit(config)
    with pytest.raises(ValueError, match=match):
        build_world(config)


def test_an_unknown_system_key_names_the_closest():
    config = build_system("sol_system").get_config_dict()
    config["worlds"]["earth"]["eccentricty"] = 0.1
    with pytest.raises(ValueError, match="did you mean 'eccentricity'"):
        build_system(config)


# =====================================================================================================================
# Exports
# =====================================================================================================================
def test_structures_exports_the_classes_and_loaders():
    assert {BaseWorld, TerrestrialWorld, GasGiantWorld, StarWorld, Layer, System} <= {
        getattr(structures, name) for name in structures.__all__}
    assert make_tide("fixed_q").model_name == "fixed_q"
    assert "io" in available_worlds()
    for removed in ("build_world_from_dict", "build_system_from_dict", "construct_world", "construct_layer"):
        assert not hasattr(structures, removed)
