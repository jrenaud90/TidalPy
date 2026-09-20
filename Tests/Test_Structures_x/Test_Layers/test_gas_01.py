"""
Tests for TidalPy.structures_x.layers.gas — GasLayer.

Covers construction, geometry/EOS/tidal inheritance (via _layer_ptr), gas
property getters, binary round-trip, TOML config save, and isinstance checks.

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

def _import_gas():
    try:
        from TidalPy.structures_x.layers import gas as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.structures_x.layers.gas not compiled — run uv pip install first."
        )


# Reference values (jovian atmosphere-ish, MKS)
_R_INNER_M   = 0.0
_R_OUTER_M   = 7.0e7     # [m]
_MASS_KG     = 1.0e27    # [kg]
_MW          = 2.0e-3    # [kg/mol] hydrogen
_GAMMA       = 1.4       # [dimensionless]
_T_REF_K     = 300.0     # [K]
_RHO_REF     = 1.2       # [kg/m³]


def _make_layer(**kw):
    mod = _import_gas()
    defaults = dict(
        name         = "atmosphere",
        layer_index  = 0,
        radius_inner = _R_INNER_M,
        radius_outer = _R_OUTER_M,
        mass         = _MASS_KG,
        material_name = "hydrogen",
        is_tidal    = False,
        tidal_scale = 1.0,
        mean_molecular_weight = _MW,
        adiabatic_index       = _GAMMA,
        reference_temperature = _T_REF_K,
        reference_density     = _RHO_REF,
    )
    defaults.update(kw)
    return mod.GasLayer(**defaults)


# =====================================================================================================================
# Construction
# =====================================================================================================================
def test_gas_construction_basic():
    """GasLayer stores all config values at construction."""
    gl = _make_layer()
    assert gl.name                         == "atmosphere"
    assert gl.layer_index                  == 0
    assert gl.mean_molecular_weight        == pytest.approx(_MW)
    assert gl.adiabatic_index              == pytest.approx(_GAMMA)
    assert gl.reference_temperature        == pytest.approx(_T_REF_K)
    assert gl.reference_density            == pytest.approx(_RHO_REF)


def test_gas_defaults():
    """GasLayer optional parameters default to expected values."""
    mod = _import_gas()
    gl  = mod.GasLayer("test", 0, 0.0, 1e7, 1e27)
    assert gl.mean_molecular_weight == pytest.approx(2.0e-3)
    assert gl.adiabatic_index       == pytest.approx(1.4)
    assert gl.reference_temperature == pytest.approx(300.0)
    assert gl.reference_density     == pytest.approx(1.0)
    # A gas carries no shear stress, so the radial solver treats it as a static liquid by default.
    assert (gl.is_solid, gl.is_static, gl.is_incompressible) == (False, True, False)
    assert gl.get_config_dict()["is_solid"] is False


def test_gas_layer_assumption_flags_from_constructor():
    gl = _make_layer(is_solid=True, is_static=False, is_incompressible=True)
    assert (gl.is_solid, gl.is_static, gl.is_incompressible) == (True, False, True)


# =====================================================================================================================
# Inheritance checks (via _layer_ptr — critical regression)
# =====================================================================================================================
def test_gas_inherits_geometry():
    """Thickness, volume, surface areas resolved correctly via _layer_ptr."""
    gl = _make_layer(radius_inner=1e7, radius_outer=7e7)
    expected_t = 7e7 - 1e7
    expected_v = (4.0 / 3.0) * math.pi * (7e7**3 - 1e7**3)
    assert gl.thickness == pytest.approx(expected_t)
    assert gl.volume    == pytest.approx(expected_v, rel=1e-9)
    assert gl.radius    == pytest.approx(7e7)


def test_gas_inherits_eos():
    """EOS operations inherited from BaseLayer work on GasLayer."""
    gl = _make_layer()
    assert gl.eos_data_populated is False
    assert math.isnan(gl.get_density(_R_OUTER_M))
    gl.update_eos_data(
        [0.0, _R_OUTER_M],
        [5.0, 1.0],
        [0.0, 25.0],
        [1e5, 0.0],
    )
    assert gl.eos_data_populated is True
    assert gl.get_density(0.0) == pytest.approx(5.0)


# =====================================================================================================================
# get_config_dict
# =====================================================================================================================
_ALL_KEYS = (
    "name", "layer_index", "radius_inner_m", "radius_outer_m", "mass_kg",
    "material_name", "is_tidal", "tidal_scale",
    "shear_modulus_static_pa", "bulk_modulus_static_pa",
    "shear_viscosity_static_pas", "bulk_viscosity_static_pas",
    "love_number_k_re", "love_number_k_im",
    "love_number_h_re", "love_number_h_im",
    "love_number_l_re", "love_number_l_im",
    "mean_molecular_weight_kg_mol", "adiabatic_index",
    "reference_temperature_k", "reference_density_kg_m3",
)


# =====================================================================================================================
def test_get_config_dict_has_all_keys():
    """get_config_dict contains all expected keys."""
    gl  = _make_layer()
    cfg = gl.get_config_dict()
    for key in _ALL_KEYS:
        assert key in cfg, f"Missing key: {key}"


def test_get_config_dict_values():
    """get_config_dict values match constructor arguments."""
    gl  = _make_layer()
    cfg = gl.get_config_dict()
    assert cfg["name"]                        == "atmosphere"
    assert cfg["mean_molecular_weight_kg_mol"] == pytest.approx(_MW)
    assert cfg["adiabatic_index"]              == pytest.approx(_GAMMA)
    assert cfg["reference_temperature_k"]      == pytest.approx(_T_REF_K)
    assert cfg["reference_density_kg_m3"]      == pytest.approx(_RHO_REF)


# =====================================================================================================================
# save_config (TOML)
# =====================================================================================================================
def test_save_config_gas():
    """save_config writes a TOML file with all GasLayer keys."""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            pytest.skip("Neither tomllib nor tomli available.")

    gl = _make_layer()
    with tempfile.NamedTemporaryFile(suffix=".toml", delete=False, mode="w") as f:
        path = f.name
    try:
        gl.save_config(path)
        with open(path, "rb") as f:
            data = tomllib.load(f)
        assert data["name"]                        == "atmosphere"
        assert data["mean_molecular_weight_kg_mol"] == pytest.approx(_MW)
        assert data["adiabatic_index"]              == pytest.approx(_GAMMA)
    finally:
        os.unlink(path)


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
def test_binary_roundtrip_gas():
    """save_binary + load_binary preserves all fields."""
    mod = _import_gas()
    gl1 = _make_layer()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        gl1.save_binary(path)
        gl2 = mod.GasLayer("placeholder", 0, 0.0, 1.0, 1.0)
        gl2.load_binary(path)
        assert gl2.name                        == "atmosphere"
        assert gl2.layer_index                 == 0
        assert gl2.radius_outer                == pytest.approx(_R_OUTER_M)
        assert gl2.mass                        == pytest.approx(_MASS_KG)
        assert gl2.mean_molecular_weight       == pytest.approx(_MW)
        assert gl2.adiabatic_index             == pytest.approx(_GAMMA)
        assert gl2.reference_temperature       == pytest.approx(_T_REF_K)
        assert gl2.reference_density           == pytest.approx(_RHO_REF)
    finally:
        os.unlink(path)


def test_binary_roundtrip_derived_fields():
    """After load_binary, derived geometry fields are recomputed correctly."""
    mod = _import_gas()
    gl1 = _make_layer(radius_inner=1e7, radius_outer=7e7)
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        gl1.save_binary(path)
        gl2 = mod.GasLayer("placeholder", 0, 0.0, 1.0, 1.0)
        gl2.load_binary(path)
        assert gl2.thickness == pytest.approx(7e7 - 1e7)
    finally:
        os.unlink(path)


def test_binary_load_file_not_found():
    """load_binary raises FileNotFoundError for a missing path."""
    gl = _make_layer()
    with pytest.raises(FileNotFoundError):
        gl.load_binary("/nonexistent/path/xyz.tpyb")


# =====================================================================================================================
# Rheology attachment and recursive binary serialization (inherited from PhysicsLayer)
# =====================================================================================================================

def _import_rheology():
    from TidalPy.rheology_x import rheology as _mod
    return _mod


def test_gas_attach_rheology_sets_flag():
    """GasLayer inherits set_shear_rheology from PhysicsLayer."""
    rheo = _import_rheology()
    gl = _make_layer(shear_modulus_static=1.0e9, shear_viscosity_static=1.0e18)
    assert gl.shear_rheology_set is False
    gl.set_shear_rheology(rheo.Maxwell())
    assert gl.shear_rheology_set is True


def test_gas_binary_roundtrip_with_rheology():
    """save_binary + load_binary restores a rheology model attached to a GasLayer."""
    mod  = _import_gas()
    rheo = _import_rheology()
    freq = 1.0e-5
    gl1  = _make_layer(shear_modulus_static=1.0e9, shear_viscosity_static=1.0e18)
    gl1.set_shear_rheology(rheo.Maxwell())
    mu_before = gl1.calc_complex_shear_modulus(freq)

    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        gl1.save_binary(path)
        gl2 = _make_layer(name="placeholder")
        gl2.load_binary(path)
        assert gl2.shear_rheology_set is True
        mu_after = gl2.calc_complex_shear_modulus(freq)
        assert mu_after.real == pytest.approx(mu_before.real, rel=1e-12)
        assert mu_after.imag == pytest.approx(mu_before.imag, rel=1e-12)
    finally:
        os.unlink(path)


# =====================================================================================================================
# isinstance checks
# =====================================================================================================================
def test_gas_is_physics_layer():
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    assert isinstance(_make_layer(), PhysicsLayer)


def test_gas_is_base_layer():
    from TidalPy.structures_x.layers.base import BaseLayer
    assert isinstance(_make_layer(), BaseLayer)


def test_gas_is_structure_base():
    from TidalPy.Utilities_x.classes_x.classes import StructureBase
    assert isinstance(_make_layer(), StructureBase)


def test_gas_is_tidalpy_base():
    from TidalPy.Utilities_x.classes_x.classes import TidalPyBaseClass
    assert isinstance(_make_layer(), TidalPyBaseClass)


def test_get_config_dict_class_name():
    """A gas layer names the builder class it was built from."""
    assert _make_layer().get_config_dict()["class"] == "gas"
