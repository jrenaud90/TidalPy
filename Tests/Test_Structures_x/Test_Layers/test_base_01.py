"""
Tests for TidalPy.structures_x.layers.base — BaseLayer (Phase 1).

Covers construction, geometry properties, EOS profile lifecycle, binary
round-trip, and TOML config save.

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
def _import_base():
    try:
        from TidalPy.structures_x.layers import base as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.structures_x.layers.base not compiled."
        )


# Reference values (MKS)
_EARTH_RADIUS_M   = 6.371e6      # [m]
_EARTH_MASS_KG    = 5.972e24     # [kg]
_MANTLE_R_INNER_M = 3.485e6      # [m]  (approx. CMB)
_MANTLE_R_OUTER_M = _EARTH_RADIUS_M
_MANTLE_MASS_KG   = 4.043e24     # [kg]


def _make_mantle():
    """Return a typical Earth-like mantle BaseLayer instance."""
    mod = _import_base()
    return mod.BaseLayer(
        name           = "mantle",
        layer_index    = 1,
        radius_inner_m = _MANTLE_R_INNER_M,
        radius_outer_m = _MANTLE_R_OUTER_M,
        mass_kg        = _MANTLE_MASS_KG,
        material_name  = "perovskite",
        is_tidal       = True,
        tidal_scale    = 1.0,
    )


# =====================================================================================================================
# Construction
# =====================================================================================================================

def test_base_layer_construction_basic():
    """BaseLayer stores all config values at construction."""
    mod = _import_base()
    bl = mod.BaseLayer("core", 0, 0.0, 3.485e6, 1.932e24)
    assert bl.name          == "core"
    assert bl.layer_index   == 0
    assert bl.radius_inner  == pytest.approx(0.0)
    assert bl.radius_outer  == pytest.approx(3.485e6)
    assert bl.mass          == pytest.approx(1.932e24)
    assert bl.material_name == ""
    assert bl.is_tidal      is True
    assert bl.tidal_scale   == pytest.approx(1.0)


def test_base_layer_construction_full_kwargs():
    """BaseLayer stores all optional fields when provided."""
    bl = _make_mantle()
    assert bl.name          == "mantle"
    assert bl.layer_index   == 1
    assert bl.radius_inner  == pytest.approx(_MANTLE_R_INNER_M)
    assert bl.radius_outer  == pytest.approx(_MANTLE_R_OUTER_M)
    assert bl.mass          == pytest.approx(_MANTLE_MASS_KG)
    assert bl.material_name == "perovskite"
    assert bl.is_tidal      is True
    assert bl.tidal_scale   == pytest.approx(1.0)


def test_base_layer_is_tidal_false():
    """is_tidal=False is stored correctly."""
    mod = _import_base()
    bl = mod.BaseLayer("gas", 2, 6e6, 7e7, 1e26, is_tidal=False, tidal_scale=0.5)
    assert bl.is_tidal    is False
    assert bl.tidal_scale == pytest.approx(0.5)


# =====================================================================================================================
# Derived geometry properties
# =====================================================================================================================
def test_base_layer_radius_alias():
    """radius property returns the outer radius."""
    bl = _make_mantle()
    assert bl.radius == pytest.approx(bl.radius_outer)


def test_base_layer_thickness():
    """Thickness equals outer minus inner radius."""
    bl = _make_mantle()
    expected = _MANTLE_R_OUTER_M - _MANTLE_R_INNER_M
    assert bl.thickness == pytest.approx(expected)


def test_base_layer_volume():
    """Volume equals the spherical shell formula."""
    bl = _make_mantle()
    v_outer = (4.0 / 3.0) * math.pi * _MANTLE_R_OUTER_M**3
    v_inner = (4.0 / 3.0) * math.pi * _MANTLE_R_INNER_M**3
    assert bl.volume == pytest.approx(v_outer - v_inner, rel=1e-9)


def test_base_layer_surface_area_outer():
    """Outer surface area equals 4*pi*r_outer^2."""
    bl = _make_mantle()
    expected = 4.0 * math.pi * _MANTLE_R_OUTER_M**2
    assert bl.surface_area_outer == pytest.approx(expected, rel=1e-9)


def test_base_layer_surface_area_inner():
    """Inner surface area equals 4*pi*r_inner^2."""
    bl = _make_mantle()
    expected = 4.0 * math.pi * _MANTLE_R_INNER_M**2
    assert bl.surface_area_inner == pytest.approx(expected, rel=1e-9)


def test_base_layer_zero_inner_radius():
    """A full-sphere layer (inner radius = 0) has zero inner surface area."""
    mod = _import_base()
    bl = mod.BaseLayer("core", 0, 0.0, 3.485e6, 1.932e24)
    assert bl.surface_area_inner == pytest.approx(0.0)
    assert bl.thickness == pytest.approx(bl.radius_outer)


# =====================================================================================================================
# Schema version (inherited from TidalPyBaseClass)
# =====================================================================================================================
def test_base_layer_schema_version():
    """BaseLayer.get_schema_version_str returns '0.2.0'."""
    bl = _make_mantle()
    assert bl.get_schema_version_str() == "0.2.0"


# =====================================================================================================================
# EOS profile lifecycle
# =====================================================================================================================
def test_eos_unpopulated_initially():
    """eos_data_populated is False on a freshly constructed layer."""
    bl = _make_mantle()
    assert bl.eos_data_populated is False


def test_eos_returns_nan_when_unpopulated():
    """get_density/gravity/pressure return NaN before EOS data is loaded."""
    import math as _math
    bl = _make_mantle()
    assert _math.isnan(bl.get_density(_MANTLE_R_INNER_M))
    assert _math.isnan(bl.get_gravity(_MANTLE_R_INNER_M))
    assert _math.isnan(bl.get_pressure(_MANTLE_R_INNER_M))


def test_eos_populated_after_update():
    """eos_data_populated is True after update_eos_data is called."""
    bl = _make_mantle()
    radii    = [_MANTLE_R_INNER_M, _MANTLE_R_OUTER_M]
    density  = [5000.0, 4000.0]
    gravity  = [10.0, 9.8]
    pressure = [1e11, 0.0]
    bl.update_eos_data(radii, density, gravity, pressure)
    assert bl.eos_data_populated is True


def test_eos_boundary_values():
    """Linear interpolation clamps to boundary values outside range."""
    bl = _make_mantle()
    radii    = [_MANTLE_R_INNER_M, _MANTLE_R_OUTER_M]
    density  = [5000.0, 4000.0]
    gravity  = [10.0, 9.8]
    pressure = [1e11, 0.0]
    bl.update_eos_data(radii, density, gravity, pressure)

    # At exact boundaries
    assert bl.get_density(_MANTLE_R_INNER_M) == pytest.approx(5000.0)
    assert bl.get_density(_MANTLE_R_OUTER_M) == pytest.approx(4000.0)

    # Below lower boundary -> clamped to first value
    assert bl.get_density(_MANTLE_R_INNER_M - 1.0) == pytest.approx(5000.0)

    # Above upper boundary -> clamped to last value
    assert bl.get_density(_MANTLE_R_OUTER_M + 1.0) == pytest.approx(4000.0)


def test_eos_midpoint_interpolation():
    """Linear interpolation returns the midpoint value at the midpoint radius."""
    bl = _make_mantle()
    r0, r1   = _MANTLE_R_INNER_M, _MANTLE_R_OUTER_M
    rho0, rho1 = 5000.0, 3000.0
    bl.update_eos_data([r0, r1], [rho0, rho1], [10.0, 9.8], [1e11, 0.0])
    r_mid  = 0.5 * (r0 + r1)
    rho_mid = 0.5 * (rho0 + rho1)
    assert bl.get_density(r_mid) == pytest.approx(rho_mid, rel=1e-12)


@pytest.mark.parametrize("n_pts", [3, 10, 100])
def test_eos_multi_point_profiles(n_pts):
    """EOS data with multiple sample points interpolates correctly."""
    import numpy as _np
    bl = _make_mantle()
    r   = list(_np.linspace(_MANTLE_R_INNER_M, _MANTLE_R_OUTER_M, n_pts))
    rho = [5000.0 - 1000.0 * (ri - _MANTLE_R_INNER_M) / (_MANTLE_R_OUTER_M - _MANTLE_R_INNER_M)
           for ri in r]
    g   = [10.0] * n_pts
    p   = [1e11 * (1.0 - (ri - _MANTLE_R_INNER_M) / (_MANTLE_R_OUTER_M - _MANTLE_R_INNER_M))
           for ri in r]
    bl.update_eos_data(r, rho, g, p)
    assert bl.eos_data_populated is True
    assert bl.get_density(r[0]) == pytest.approx(rho[0], rel=1e-9)
    assert bl.get_density(r[-1]) == pytest.approx(rho[-1], rel=1e-9)


# =====================================================================================================================
# StructureBase geometry calc methods (inherited; pure-function, not stored state)
# =====================================================================================================================
def test_calc_surface_area_via_base_layer():
    """calc_surface_area inherited from StructureBase works on BaseLayer instances."""
    bl = _make_mantle()
    r  = 1.0
    assert bl.calc_surface_area(r) == pytest.approx(4.0 * math.pi * r * r)


def test_calc_volume_sphere_via_base_layer():
    """calc_volume_sphere inherited from StructureBase works on BaseLayer instances."""
    bl = _make_mantle()
    r  = 1.0
    assert bl.calc_volume_sphere(r) == pytest.approx(4.0 * math.pi * r**3 / 3.0)


# =====================================================================================================================
# get_config_dict
# =====================================================================================================================
def test_get_config_dict_keys():
    """get_config_dict returns all expected keys."""
    bl  = _make_mantle()
    cfg = bl.get_config_dict()
    for key in ("name", "layer_index", "radius_inner_m", "radius_outer_m",
                "mass_kg", "material_name", "is_tidal", "tidal_scale"):
        assert key in cfg, f"Missing key: {key}"


def test_get_config_dict_values():
    """get_config_dict values match constructor arguments."""
    bl  = _make_mantle()
    cfg = bl.get_config_dict()
    assert cfg["name"]           == "mantle"
    assert cfg["layer_index"]    == 1
    assert cfg["radius_inner_m"] == pytest.approx(_MANTLE_R_INNER_M)
    assert cfg["radius_outer_m"] == pytest.approx(_MANTLE_R_OUTER_M)
    assert cfg["mass_kg"]        == pytest.approx(_MANTLE_MASS_KG)
    assert cfg["material_name"]  == "perovskite"
    assert cfg["is_tidal"]       is True
    assert cfg["tidal_scale"]    == pytest.approx(1.0)


# =====================================================================================================================
# save_config (TOML)
# =====================================================================================================================
def test_save_config_writes_toml():
    """save_config writes a TOML file containing all geometry keys."""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            pytest.skip("Neither tomllib nor tomli available.")

    bl = _make_mantle()
    with tempfile.NamedTemporaryFile(suffix=".toml", delete=False, mode="w") as f:
        path = f.name
    try:
        bl.save_config(path)
        with open(path, "rb") as f:
            data = tomllib.load(f)
        assert data["name"]           == "mantle"
        assert data["radius_inner_m"] == pytest.approx(_MANTLE_R_INNER_M)
        assert data["radius_outer_m"] == pytest.approx(_MANTLE_R_OUTER_M)
        assert data["mass_kg"]        == pytest.approx(_MANTLE_MASS_KG)
    finally:
        os.unlink(path)


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
def test_binary_roundtrip():
    """save_binary + load_binary preserves all geometry fields."""
    mod = _import_base()
    bl1 = _make_mantle()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        bl1.save_binary(path)

        # Create a fresh instance with different values, then load into it.
        bl2 = mod.BaseLayer("placeholder", 99, 0.0, 1.0, 1.0)
        bl2.load_binary(path)

        assert bl2.name           == "mantle"
        assert bl2.layer_index    == 1
        assert bl2.radius_inner   == pytest.approx(_MANTLE_R_INNER_M)
        assert bl2.radius_outer   == pytest.approx(_MANTLE_R_OUTER_M)
        assert bl2.mass           == pytest.approx(_MANTLE_MASS_KG)
        assert bl2.material_name  == "perovskite"
        assert bl2.is_tidal       is True
        assert bl2.tidal_scale    == pytest.approx(1.0)
    finally:
        os.unlink(path)


def test_binary_roundtrip_derived_fields():
    """After load_binary, derived geometry fields (thickness, volume, areas) are correct."""
    mod = _import_base()
    bl1 = _make_mantle()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        bl1.save_binary(path)
        bl2 = mod.BaseLayer("placeholder", 0, 0.0, 1.0, 1.0)
        bl2.load_binary(path)

        expected_thickness = _MANTLE_R_OUTER_M - _MANTLE_R_INNER_M
        expected_volume    = (4.0 / 3.0) * math.pi * (
            _MANTLE_R_OUTER_M**3 - _MANTLE_R_INNER_M**3)

        assert bl2.thickness == pytest.approx(expected_thickness, rel=1e-9)
        assert bl2.volume    == pytest.approx(expected_volume, rel=1e-9)
    finally:
        os.unlink(path)


def test_binary_roundtrip_eos_not_preserved():
    """EOS data is NOT serialized; eos_data_populated is False after load_binary."""
    mod = _import_base()
    bl1 = _make_mantle()
    bl1.update_eos_data(
        [_MANTLE_R_INNER_M, _MANTLE_R_OUTER_M],
        [5000.0, 4000.0],
        [10.0, 9.8],
        [1e11, 0.0],
    )
    assert bl1.eos_data_populated is True

    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        bl1.save_binary(path)
        bl2 = mod.BaseLayer("placeholder", 0, 0.0, 1.0, 1.0)
        bl2.load_binary(path)
        # EOS data should NOT be restored from binary.
        assert bl2.eos_data_populated is False
    finally:
        os.unlink(path)


def test_binary_load_file_not_found():
    """load_binary raises FileNotFoundError for a missing path."""
    bl = _make_mantle()
    with pytest.raises(FileNotFoundError):
        bl.load_binary("/nonexistent/path/xyz.tpyb")


# =====================================================================================================================
# isinstance checks
# =====================================================================================================================
def test_base_layer_is_structure_base():
    """BaseLayer is an instance of StructureBase."""
    from TidalPy.Utilities_x.classes_x.classes import StructureBase
    bl = _make_mantle()
    assert isinstance(bl, StructureBase)


def test_base_layer_is_tidalpy_base():
    """BaseLayer is an instance of TidalPyBaseClass."""
    from TidalPy.Utilities_x.classes_x.classes import TidalPyBaseClass
    bl = _make_mantle()
    assert isinstance(bl, TidalPyBaseClass)
