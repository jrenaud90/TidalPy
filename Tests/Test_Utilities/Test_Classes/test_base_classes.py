"""Tests for TidalPyBaseClass, StructureBase, and PhysicsBase in ``TidalPy.Utilities.classes.classes``."""
import math

import pytest

from TidalPy.Utilities.classes import classes

_EARTH_RADIUS_M = 6.371e6
_EARTH_MASS_KG  = 5.972e24
_G              = 6.674e-11


def test_tidalpy_base_not_directly_instantiable():
    """TidalPyBaseClass raises RuntimeError on a method call when it wraps no C++ object."""
    base = classes.TidalPyBaseClass()
    with pytest.raises(RuntimeError):
        base.get_schema_version_str()


@pytest.mark.parametrize(
    "make_object",
    [lambda: classes.StructureBase(1.0, 1.0), lambda: classes.PhysicsBase("elastic")],
    ids=["StructureBase", "PhysicsBase"])
def test_schema_version(make_object):
    """get_schema_version_str returns '0.2.0'."""
    assert make_object().get_schema_version_str() == "0.2.0"


def test_structure_base_construction():
    """StructureBase stores radius and mass."""
    structure = classes.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    assert structure.radius == pytest.approx(_EARTH_RADIUS_M)
    assert structure.mass   == pytest.approx(_EARTH_MASS_KG)


def test_structure_base_get_config_dict():
    """get_config_dict returns the radius and mass keys."""
    structure = classes.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    cfg = structure.get_config_dict()
    assert "radius_m" in cfg
    assert "mass_kg"  in cfg
    assert cfg["radius_m"] == pytest.approx(_EARTH_RADIUS_M)
    assert cfg["mass_kg"]  == pytest.approx(_EARTH_MASS_KG)


@pytest.mark.parametrize(
    "structure_args, method, args, expected, rel",
    [
        pytest.param((1.0, 1.0), "calc_surface_area", (1.0,), 4.0 * math.pi, None, id="unit_surface_area"),
        pytest.param((1.0, 1.0), "calc_volume_sphere", (1.0,), 4.0 * math.pi / 3.0, None, id="unit_volume"),
        pytest.param(
            (_EARTH_RADIUS_M, _EARTH_MASS_KG),
            "calc_surface_gravity",
            (_EARTH_MASS_KG, _EARTH_RADIUS_M),
            _G * _EARTH_MASS_KG / _EARTH_RADIUS_M**2,
            1e-3,
            id="earth_surface_gravity"),
        pytest.param(
            (_EARTH_RADIUS_M, _EARTH_MASS_KG),
            "calc_escape_velocity",
            (_EARTH_MASS_KG, _EARTH_RADIUS_M),
            11186.0,
            0.01,
            id="earth_escape_velocity"),
        pytest.param((0.0, 0.0), "calc_surface_area", (0.0,), 0.0, None, id="zero_surface_area"),
    ])
def test_structure_geometry(structure_args, method, args, expected, rel):
    """StructureBase geometry helpers return the expected values."""
    structure = classes.StructureBase(*structure_args)
    assert getattr(structure, method)(*args) == pytest.approx(expected, rel=rel)


def test_calc_volume_shell():
    """Shell volume equals the outer sphere volume minus the inner sphere volume."""
    structure = classes.StructureBase(1.0, 1.0)
    expected = structure.calc_volume_sphere(2.0) - structure.calc_volume_sphere(1.0)
    assert structure.calc_volume_shell(2.0, 1.0) == pytest.approx(expected)


def test_calc_mean_density():
    """Earth's mean density is about 5515 kg/m^3."""
    structure = classes.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    volume = structure.calc_volume_sphere(_EARTH_RADIUS_M)
    assert structure.calc_mean_density(_EARTH_MASS_KG, volume) == pytest.approx(5515.0, rel=0.01)


def test_calc_surface_gravity_zero_radius():
    """calc_surface_gravity returns exactly 0 for a zero radius."""
    structure = classes.StructureBase(0.0, 1.0)
    assert structure.calc_surface_gravity(1.0, 0.0) == 0.0


def test_structure_base_binary_roundtrip(tmp_path):
    """save_binary then load_binary preserves radius and mass."""
    path = str(tmp_path / "structure.tpyb")
    classes.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG).save_binary(path)
    reloaded = classes.StructureBase(0.0, 0.0)
    reloaded.load_binary(path)
    assert reloaded.radius == pytest.approx(_EARTH_RADIUS_M)
    assert reloaded.mass   == pytest.approx(_EARTH_MASS_KG)


def test_structure_base_load_binary_not_found():
    """load_binary raises FileNotFoundError for a missing path."""
    structure = classes.StructureBase(1.0, 1.0)
    with pytest.raises(FileNotFoundError):
        structure.load_binary("/nonexistent/path/xyz.tpyb")


def test_structure_base_save_config(tmp_path):
    """save_config writes a TOML file holding radius and mass."""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            pytest.skip("Neither tomllib nor tomli available for reading TOML.")

    path = str(tmp_path / "structure.toml")
    classes.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG).save_config(path)
    with open(path, 'rb') as file:
        data = tomllib.load(file)
    assert data["radius_m"] == pytest.approx(_EARTH_RADIUS_M)
    assert data["mass_kg"]  == pytest.approx(_EARTH_MASS_KG)


def test_save_config_writes_the_version_header_with_lf(tmp_path):
    """save_config starts with the version header every saved configuration carries, and writes LF newlines."""
    path = tmp_path / "structure.toml"
    classes.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG).save_config(str(path))
    raw = path.read_bytes()
    assert b"\r\n" not in raw
    text = raw.decode("utf-8")
    assert text.startswith("# ===")
    assert "#  TidalPy StructureBase configuration" in text
    assert "#  TidalPy version:" in text


@pytest.mark.parametrize(
    "name",
    ["elastic", "viscous", "voigt", "maxwell", "burgers", "andrade", "sundberg", "", "a" * 256],
    ids=["elastic", "viscous", "voigt", "maxwell", "burgers", "andrade", "sundberg", "empty", "long"])
def test_physics_base_model_names(name):
    """PhysicsBase stores and returns any model name."""
    assert classes.PhysicsBase(name).model_name == name


def test_physics_base_model_name_is_read_only():
    """model_name cannot be reassigned, since it names what the model computes."""
    physics = classes.PhysicsBase("maxwell")
    with pytest.raises(AttributeError):
        physics.model_name = "andrade"
    assert physics.model_name == "maxwell"


def test_physics_base_get_config_dict():
    """get_config_dict returns the builder's model key."""
    cfg = classes.PhysicsBase("voigt").get_config_dict()
    assert "model" in cfg
    assert cfg["model"] == "voigt"


@pytest.mark.parametrize(
    "name, placeholder",
    [("andrade", "empty"), ("", "placeholder")],
    ids=["andrade", "empty_name"])
def test_physics_base_binary_roundtrip(tmp_path, name, placeholder):
    """save_binary then load_binary preserves model_name."""
    path = str(tmp_path / "physics.tpyb")
    classes.PhysicsBase(name).save_binary(path)
    reloaded = classes.PhysicsBase(placeholder)
    reloaded.load_binary(path)
    assert reloaded.model_name == name
