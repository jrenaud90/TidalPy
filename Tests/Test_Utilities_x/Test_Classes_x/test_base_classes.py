"""
Tests for TidalPy.Utilities_x.classes_x.classes

Covers TidalPyBaseClass, StructureBase, PhysicsBase.
Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import math
import os
import tempfile
import pytest


# =====================================================================================================================
# Helpers
# =====================================================================================================================
def _import_classes():
    try:
        from TidalPy.Utilities_x.classes_x import classes as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.Utilities_x.classes_x.classes not compiled — run uv pip install first."
        )


# Known physical values used across tests (MKS)
_EARTH_RADIUS_M = 6.371e6
_EARTH_MASS_KG  = 5.972e24
_G              = 6.674e-11


# =====================================================================================================================
# TidalPyBaseClass — cannot instantiate directly
# =====================================================================================================================
def test_tidalpy_base_not_directly_instantiable():
    """TidalPyBaseClass raises RuntimeError on any method call when _ptr is NULL."""
    mod = _import_classes()
    base = mod.TidalPyBaseClass()
    with pytest.raises(RuntimeError):
        base.get_schema_version_str()


# =====================================================================================================================
# StructureBase — construction and properties
# =====================================================================================================================
def test_structure_base_construction():
    """StructureBase stores radius and mass correctly."""
    mod = _import_classes()
    sb = mod.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    assert sb.radius == pytest.approx(_EARTH_RADIUS_M)
    assert sb.mass   == pytest.approx(_EARTH_MASS_KG)


def test_structure_base_schema_version():
    """StructureBase.get_schema_version_str returns '0.2.0'."""
    mod = _import_classes()
    sb = mod.StructureBase(1.0, 1.0)
    assert sb.get_schema_version_str() == "0.2.0"


def test_structure_base_get_config_dict():
    """get_config_dict returns correct radius_m and mass_kg keys."""
    mod = _import_classes()
    sb = mod.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    cfg = sb.get_config_dict()
    assert "radius_m" in cfg
    assert "mass_kg"  in cfg
    assert cfg["radius_m"] == pytest.approx(_EARTH_RADIUS_M)
    assert cfg["mass_kg"]  == pytest.approx(_EARTH_MASS_KG)


# =====================================================================================================================
# StructureBase — geometry calculations
# =====================================================================================================================
def test_calc_surface_area():
    """Surface area of a unit sphere is 4*pi."""
    mod = _import_classes()
    sb = mod.StructureBase(1.0, 1.0)
    assert sb.calc_surface_area(1.0) == pytest.approx(4.0 * math.pi)


def test_calc_volume_sphere():
    """Volume of a unit sphere is 4*pi/3."""
    mod = _import_classes()
    sb = mod.StructureBase(1.0, 1.0)
    assert sb.calc_volume_sphere(1.0) == pytest.approx(4.0 * math.pi / 3.0)


def test_calc_volume_shell():
    """Volume of shell = volume of outer sphere minus volume of inner sphere."""
    mod = _import_classes()
    sb = mod.StructureBase(1.0, 1.0)
    expected = sb.calc_volume_sphere(2.0) - sb.calc_volume_sphere(1.0)
    assert sb.calc_volume_shell(2.0, 1.0) == pytest.approx(expected)


def test_calc_surface_gravity_earth():
    """Surface gravity of Earth is approximately 9.8 m/s^2."""
    mod = _import_classes()
    sb = mod.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    g = sb.calc_surface_gravity(_EARTH_MASS_KG, _EARTH_RADIUS_M)
    assert g == pytest.approx(_G * _EARTH_MASS_KG / _EARTH_RADIUS_M**2, rel=1e-3)


def test_calc_mean_density():
    """Mean density of Earth is approximately 5500 kg/m^3."""
    mod = _import_classes()
    sb = mod.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    vol = sb.calc_volume_sphere(_EARTH_RADIUS_M)
    rho = sb.calc_mean_density(_EARTH_MASS_KG, vol)
    assert rho == pytest.approx(5515.0, rel=0.01)


def test_calc_escape_velocity_earth():
    """Escape velocity of Earth is approximately 11.2 km/s."""
    mod = _import_classes()
    sb = mod.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    v_esc = sb.calc_escape_velocity(_EARTH_MASS_KG, _EARTH_RADIUS_M)
    assert v_esc == pytest.approx(11186.0, rel=0.01)


def test_calc_surface_area_zero_radius():
    """calc_surface_area(0) returns 0."""
    mod = _import_classes()
    sb = mod.StructureBase(0.0, 0.0)
    assert sb.calc_surface_area(0.0) == pytest.approx(0.0)


def test_calc_surface_gravity_zero_radius():
    """calc_surface_gravity returns 0 when radius is 0."""
    mod = _import_classes()
    sb = mod.StructureBase(0.0, 1.0)
    assert sb.calc_surface_gravity(1.0, 0.0) == 0.0


# =====================================================================================================================
# StructureBase — binary round-trip
# =====================================================================================================================
def test_structure_base_binary_roundtrip():
    """save_binary + load_binary preserves radius and mass."""
    mod = _import_classes()
    sb1 = mod.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        sb1.save_binary(path)
        sb2 = mod.StructureBase(0.0, 0.0)
        sb2.load_binary(path)
        assert sb2.radius == pytest.approx(_EARTH_RADIUS_M)
        assert sb2.mass   == pytest.approx(_EARTH_MASS_KG)
    finally:
        os.unlink(path)


def test_structure_base_load_binary_not_found():
    """load_binary raises FileNotFoundError for missing path."""
    mod = _import_classes()
    sb = mod.StructureBase(1.0, 1.0)
    with pytest.raises(FileNotFoundError):
        sb.load_binary("/nonexistent/path/xyz.tpyb")


# =====================================================================================================================
# StructureBase — TOML config round-trip
# =====================================================================================================================
def test_structure_base_save_config():
    """save_config writes a valid TOML file containing radius_m and mass_kg."""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            pytest.skip("Neither tomllib nor tomli available for reading TOML.")

    mod = _import_classes()
    sb = mod.StructureBase(_EARTH_RADIUS_M, _EARTH_MASS_KG)
    with tempfile.NamedTemporaryFile(suffix=".toml", delete=False, mode='w') as f:
        path = f.name
    try:
        sb.save_config(path)
        with open(path, 'rb') as f:
            data = tomllib.load(f)
        assert data["radius_m"] == pytest.approx(_EARTH_RADIUS_M)
        assert data["mass_kg"]  == pytest.approx(_EARTH_MASS_KG)
    finally:
        os.unlink(path)


# =====================================================================================================================
# PhysicsBase — construction and properties
# =====================================================================================================================
def test_physics_base_construction():
    """PhysicsBase stores model_name correctly."""
    mod = _import_classes()
    pb = mod.PhysicsBase("maxwell")
    assert pb.model_name == "maxwell"


def test_physics_base_schema_version():
    """PhysicsBase.get_schema_version_str returns '0.2.0'."""
    mod = _import_classes()
    pb = mod.PhysicsBase("elastic")
    assert pb.get_schema_version_str() == "0.2.0"


def test_physics_base_model_name_setter():
    """PhysicsBase.model_name setter updates the stored name."""
    mod = _import_classes()
    pb = mod.PhysicsBase("maxwell")
    pb.model_name = "andrade"
    assert pb.model_name == "andrade"


def test_physics_base_get_config_dict():
    """get_config_dict returns correct model_name key."""
    mod = _import_classes()
    pb = mod.PhysicsBase("voigt")
    cfg = pb.get_config_dict()
    assert "model_name" in cfg
    assert cfg["model_name"] == "voigt"


@pytest.mark.parametrize("name", [
    "elastic", "viscous", "voigt", "maxwell", "burgers", "andrade", "sundberg",
    "",  # empty string is valid
    "a" * 256,  # long name
])
def test_physics_base_model_names(name):
    """PhysicsBase correctly stores and retrieves various model names."""
    mod = _import_classes()
    pb = mod.PhysicsBase(name)
    assert pb.model_name == name


# =====================================================================================================================
# PhysicsBase — binary round-trip
# =====================================================================================================================
def test_physics_base_binary_roundtrip():
    """save_binary + load_binary preserves model_name."""
    mod = _import_classes()
    pb1 = mod.PhysicsBase("andrade")
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        pb1.save_binary(path)
        pb2 = mod.PhysicsBase("empty")
        pb2.load_binary(path)
        assert pb2.model_name == "andrade"
    finally:
        os.unlink(path)


def test_physics_base_binary_roundtrip_empty_name():
    """Binary round-trip works with an empty model name."""
    mod = _import_classes()
    pb1 = mod.PhysicsBase("")
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        pb1.save_binary(path)
        pb2 = mod.PhysicsBase("placeholder")
        pb2.load_binary(path)
        assert pb2.model_name == ""
    finally:
        os.unlink(path)
