"""
Tests for TidalPy.structures_x.layers.physics — PhysicsLayer (Phase 2).

Covers construction, geometry inheritance, mechanical property getters,
tidal susceptibility, complex modulus calculations (no rheology), binary
round-trip, TOML config save, and isinstance checks.

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

def _import_physics():
    try:
        from TidalPy.structures_x.layers import physics as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.structures_x.layers.physics not compiled — run uv pip install first."
        )


# Reference values (MKS)
_MANTLE_R_INNER_M = 3.485e6      # [m]  (approx. CMB)
_MANTLE_R_OUTER_M = 6.371e6      # [m]  (Earth surface)
_MANTLE_MASS_KG   = 4.043e24     # [kg]
_SHEAR_MOD_PA     = 1.67e11      # [Pa]  (perovskite, lower mantle)
_BULK_MOD_PA      = 3.57e11      # [Pa]
_VISCOSITY_PAS    = 1.0e21       # [Pa·s]
_G                = 6.674e-11    # [m^3 kg^-1 s^-2]


def _make_mantle(shear=_SHEAR_MOD_PA, bulk=_BULK_MOD_PA, shear_visc=_VISCOSITY_PAS, bulk_visc=_VISCOSITY_PAS):
    mod = _import_physics()
    return mod.PhysicsLayer(
        name                       = "mantle",
        layer_index                = 1,
        radius_inner_m             = _MANTLE_R_INNER_M,
        radius_outer_m             = _MANTLE_R_OUTER_M,
        mass_kg                    = _MANTLE_MASS_KG,
        material_name              = "perovskite",
        is_tidal                   = True,
        tidal_scale                = 1.0,
        shear_modulus_static_pa    = shear,
        bulk_modulus_static_pa     = bulk,
        shear_viscosity_static_pas = shear_visc,
        bulk_viscosity_static_pas  = bulk_visc,
        love_number_k              = 0.0 + 0.0j,
        love_number_h              = 0.0 + 0.0j,
        love_number_l              = 0.0 + 0.0j,
    )


# =====================================================================================================================
# Construction
# =====================================================================================================================
def test_physics_layer_construction_basic():
    """PhysicsLayer stores all config values at construction."""
    mod = _import_physics()
    pl = mod.PhysicsLayer("core", 0, 0.0, 3.485e6, 1.932e24,
                          shear_modulus_static_pa=5e10, bulk_modulus_static_pa=2e11,
                          shear_viscosity_static_pas=1e20, bulk_viscosity_static_pas=1e22)
    assert pl.name                   == "core"
    assert pl.layer_index            == 0
    assert pl.radius_inner           == pytest.approx(0.0)
    assert pl.radius_outer           == pytest.approx(3.485e6)
    assert pl.mass                   == pytest.approx(1.932e24)
    assert pl.shear_modulus_static   == pytest.approx(5e10)
    assert pl.bulk_modulus_static    == pytest.approx(2e11)
    assert pl.shear_viscosity_static == pytest.approx(1e20)
    assert pl.bulk_viscosity_static  == pytest.approx(1e22)
    assert pl.love_number_k          == pytest.approx(0.0 + 0.0j)
    assert pl.love_number_h          == pytest.approx(0.0 + 0.0j)
    assert pl.love_number_l          == pytest.approx(0.0 + 0.0j)


def test_physics_layer_defaults():
    """PhysicsLayer optional mechanical parameters default to 0.0."""
    mod = _import_physics()
    pl = mod.PhysicsLayer("test", 0, 0.0, 1e6, 1e20)
    assert pl.shear_modulus_static   == pytest.approx(0.0)
    assert pl.bulk_modulus_static    == pytest.approx(0.0)
    assert pl.shear_viscosity_static == pytest.approx(0.0)
    assert pl.bulk_viscosity_static  == pytest.approx(0.0)
    assert pl.love_number_k          == pytest.approx(0.0 + 0.0j)
    assert pl.love_number_h          == pytest.approx(0.0 + 0.0j)
    assert pl.love_number_l          == pytest.approx(0.0 + 0.0j)


# =====================================================================================================================
# Inherited geometry properties (via _layer_ptr refactor — critical regression check)
# =====================================================================================================================
def test_physics_layer_inherits_name():
    """Geometry properties inherited from BaseLayer work on PhysicsLayer instances."""
    pl = _make_mantle()
    assert pl.name          == "mantle"
    assert pl.layer_index   == 1
    assert pl.material_name == "perovskite"
    assert pl.is_tidal      is True
    assert pl.tidal_scale   == pytest.approx(1.0)


def test_physics_layer_inherits_geometry():
    """Thickness, volume, and surface areas computed correctly via BaseLayer getters."""
    pl = _make_mantle()
    expected_thickness = _MANTLE_R_OUTER_M - _MANTLE_R_INNER_M
    expected_volume    = (4.0 / 3.0) * math.pi * (_MANTLE_R_OUTER_M**3 - _MANTLE_R_INNER_M**3)
    expected_sa_outer  = 4.0 * math.pi * _MANTLE_R_OUTER_M**2
    expected_sa_inner  = 4.0 * math.pi * _MANTLE_R_INNER_M**2

    assert pl.radius           == pytest.approx(_MANTLE_R_OUTER_M)
    assert pl.mass             == pytest.approx(_MANTLE_MASS_KG)
    assert pl.thickness        == pytest.approx(expected_thickness)
    assert pl.volume           == pytest.approx(expected_volume, rel=1e-9)
    assert pl.surface_area_outer == pytest.approx(expected_sa_outer, rel=1e-9)
    assert pl.surface_area_inner == pytest.approx(expected_sa_inner, rel=1e-9)


def test_physics_layer_eos_inherited():
    """EOS profile operations (inherited from BaseLayer) work correctly."""
    import math as _math
    pl = _make_mantle()
    assert pl.eos_data_populated is False
    assert _math.isnan(pl.get_density(_MANTLE_R_INNER_M))

    pl.update_eos_data(
        [_MANTLE_R_INNER_M, _MANTLE_R_OUTER_M],
        [5000.0, 4000.0],
        [10.0,   9.81],
        [1e11,   0.0],
    )
    assert pl.eos_data_populated is True
    assert pl.get_density(_MANTLE_R_INNER_M) == pytest.approx(5000.0)


# =====================================================================================================================
# Tidal susceptibility
# =====================================================================================================================
def test_tidal_susceptibility_formula():
    """calc_tidal_susceptibility matches (3/2) * r^5 / (G * m^2) analytically."""
    pl = _make_mantle()
    r  = _MANTLE_R_OUTER_M
    m  = _MANTLE_MASS_KG
    expected = (1.5 * r**5) / (_G * m**2)
    # Allow 0.1 % relative tolerance to account for the G value in TidalPy config.
    assert pl.calc_tidal_susceptibility() == pytest.approx(expected, rel=1e-3)


def test_tidal_susceptibility_zero_mass():
    """calc_tidal_susceptibility returns 0.0 when mass is zero."""
    mod = _import_physics()
    pl  = mod.PhysicsLayer("test", 0, 0.0, 1e6, 0.0)
    assert pl.calc_tidal_susceptibility() == pytest.approx(0.0)


# =====================================================================================================================
# Complex moduli (no rheology attached)
# =====================================================================================================================

def test_complex_shear_no_rheology_real_only():
    """Without rheology, calc_complex_shear_modulus returns purely real value."""
    pl  = _make_mantle()
    mu  = pl.calc_complex_shear_modulus(1e-5)
    assert mu.real == pytest.approx(_SHEAR_MOD_PA)
    assert mu.imag == pytest.approx(0.0)


def test_complex_bulk_no_rheology_real_only():
    """Without rheology, calc_complex_bulk_modulus returns purely real value."""
    pl = _make_mantle()
    K  = pl.calc_complex_bulk_modulus(1e-5)
    assert K.real == pytest.approx(_BULK_MOD_PA)
    assert K.imag == pytest.approx(0.0)


@pytest.mark.parametrize("freq", [1e-10, 1e-5, 1.0, 1e5])
def test_complex_shear_frequency_independent_without_rheology(freq):
    """Without rheology, complex shear modulus is frequency-independent."""
    pl = _make_mantle()
    mu = pl.calc_complex_shear_modulus(freq)
    assert mu.real == pytest.approx(_SHEAR_MOD_PA)
    assert mu.imag == pytest.approx(0.0)


def test_complex_modulus_zero_modulus():
    """With zero static modulus and no rheology, complex modulus is zero."""
    mod = _import_physics()
    pl  = mod.PhysicsLayer("test", 0, 0.0, 1e6, 1e20,
                           shear_modulus_static_pa=0.0, bulk_modulus_static_pa=0.0)
    mu = pl.calc_complex_shear_modulus(1e-5)
    K  = pl.calc_complex_bulk_modulus(1e-5)
    assert mu.real == pytest.approx(0.0)
    assert K.real  == pytest.approx(0.0)


# =====================================================================================================================
# Rheology set flags
# =====================================================================================================================

def test_rheology_not_set_initially():
    """shear_rheology_set and bulk_rheology_set are False on a new layer."""
    pl = _make_mantle()
    assert pl.shear_rheology_set is False
    assert pl.bulk_rheology_set  is False


# =====================================================================================================================
# get_config_dict
# =====================================================================================================================

def test_get_config_dict_has_all_keys():
    """get_config_dict returns all BaseLayer keys plus all physics keys."""
    pl  = _make_mantle()
    cfg = pl.get_config_dict()
    for key in ("name", "layer_index", "radius_inner_m", "radius_outer_m",
                "mass_kg", "material_name", "is_tidal", "tidal_scale",
                "shear_modulus_static_pa", "bulk_modulus_static_pa",
                "shear_viscosity_static_pas", "bulk_viscosity_static_pas",
                "love_number_k_re", "love_number_k_im",
                "love_number_h_re", "love_number_h_im",
                "love_number_l_re", "love_number_l_im"):
        assert key in cfg, f"Missing key: {key}"


def test_get_config_dict_values():
    """get_config_dict values match constructor arguments."""
    pl  = _make_mantle()
    cfg = pl.get_config_dict()
    assert cfg["name"]            == "mantle"
    assert cfg["layer_index"]     == 1
    assert cfg["radius_inner_m"]  == pytest.approx(_MANTLE_R_INNER_M)
    assert cfg["radius_outer_m"]  == pytest.approx(_MANTLE_R_OUTER_M)
    assert cfg["mass_kg"]         == pytest.approx(_MANTLE_MASS_KG)
    assert cfg["shear_modulus_static_pa"]    == pytest.approx(_SHEAR_MOD_PA)
    assert cfg["bulk_modulus_static_pa"]     == pytest.approx(_BULK_MOD_PA)
    assert cfg["shear_viscosity_static_pas"] == pytest.approx(_VISCOSITY_PAS)
    assert cfg["bulk_viscosity_static_pas"]  == pytest.approx(_VISCOSITY_PAS)
    assert cfg["love_number_k_re"]            == pytest.approx(0.0)
    assert cfg["love_number_k_im"]            == pytest.approx(0.0)
    assert cfg["love_number_h_re"]            == pytest.approx(0.0)
    assert cfg["love_number_l_re"]            == pytest.approx(0.0)


# =====================================================================================================================
# save_config (TOML)
# =====================================================================================================================

def test_save_config_physics_layer():
    """save_config writes a TOML file containing all physics keys."""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            pytest.skip("Neither tomllib nor tomli available.")

    pl = _make_mantle()
    with tempfile.NamedTemporaryFile(suffix=".toml", delete=False, mode="w") as f:
        path = f.name
    try:
        pl.save_config(path)
        with open(path, "rb") as f:
            data = tomllib.load(f)
        assert data["name"]                       == "mantle"
        assert data["shear_modulus_static_pa"]    == pytest.approx(_SHEAR_MOD_PA)
        assert data["bulk_modulus_static_pa"]     == pytest.approx(_BULK_MOD_PA)
        assert data["shear_viscosity_static_pas"] == pytest.approx(_VISCOSITY_PAS)
        assert data["bulk_viscosity_static_pas"]  == pytest.approx(_VISCOSITY_PAS)
        assert data["love_number_k_re"]            == pytest.approx(0.0)
        assert data["love_number_k_im"]            == pytest.approx(0.0)
        assert data["love_number_h_re"]            == pytest.approx(0.0)
        assert data["love_number_l_re"]            == pytest.approx(0.0)
    finally:
        os.unlink(path)


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================

def test_binary_roundtrip_physics():
    """save_binary + load_binary preserves all geometry and physics fields."""
    mod = _import_physics()
    pl1 = _make_mantle()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        pl1.save_binary(path)
        pl2 = mod.PhysicsLayer("placeholder", 99, 0.0, 1.0, 1.0)
        pl2.load_binary(path)

        assert pl2.name                   == "mantle"
        assert pl2.layer_index            == 1
        assert pl2.radius_inner           == pytest.approx(_MANTLE_R_INNER_M)
        assert pl2.radius_outer           == pytest.approx(_MANTLE_R_OUTER_M)
        assert pl2.mass                   == pytest.approx(_MANTLE_MASS_KG)
        assert pl2.material_name          == "perovskite"
        assert pl2.shear_modulus_static   == pytest.approx(_SHEAR_MOD_PA)
        assert pl2.bulk_modulus_static    == pytest.approx(_BULK_MOD_PA)
        assert pl2.shear_viscosity_static == pytest.approx(_VISCOSITY_PAS)
        assert pl2.bulk_viscosity_static  == pytest.approx(_VISCOSITY_PAS)
        assert pl2.love_number_k          == pytest.approx(0.0 + 0.0j)
        assert pl2.love_number_h          == pytest.approx(0.0 + 0.0j)
        assert pl2.love_number_l          == pytest.approx(0.0 + 0.0j)
    finally:
        os.unlink(path)


def test_binary_roundtrip_derived_fields():
    """After load_binary, derived geometry fields are recomputed correctly."""
    mod = _import_physics()
    pl1 = _make_mantle()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        pl1.save_binary(path)
        pl2 = mod.PhysicsLayer("placeholder", 0, 0.0, 1.0, 1.0)
        pl2.load_binary(path)

        expected_thickness = _MANTLE_R_OUTER_M - _MANTLE_R_INNER_M
        expected_volume    = (4.0 / 3.0) * math.pi * (
            _MANTLE_R_OUTER_M**3 - _MANTLE_R_INNER_M**3)
        assert pl2.thickness == pytest.approx(expected_thickness, rel=1e-9)
        assert pl2.volume    == pytest.approx(expected_volume,    rel=1e-9)
    finally:
        os.unlink(path)


def test_binary_roundtrip_complex_modulus_preserved():
    """After load_binary, complex modulus uses restored static modulus."""
    mod = _import_physics()
    pl1 = _make_mantle()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        pl1.save_binary(path)
        pl2 = mod.PhysicsLayer("placeholder", 0, 0.0, 1.0, 1.0)
        pl2.load_binary(path)
        mu = pl2.calc_complex_shear_modulus(1e-5)
        assert mu.real == pytest.approx(_SHEAR_MOD_PA)
        assert mu.imag == pytest.approx(0.0)
    finally:
        os.unlink(path)


def test_binary_load_file_not_found():
    """load_binary raises FileNotFoundError for a missing path."""
    pl = _make_mantle()
    with pytest.raises(FileNotFoundError):
        pl.load_binary("/nonexistent/path/xyz.tpyb")


# =====================================================================================================================
# Rheology attachment and recursive binary serialization
# =====================================================================================================================

def _import_rheology():
    try:
        from TidalPy.rheology_x import rheology as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.rheology_x.rheology not compiled — run uv pip install first."
        )


def test_attach_shear_rheology_sets_flag():
    """set_shear_rheology flips the shear_rheology_set flag and consumes the model."""
    rheo = _import_rheology()
    pl   = _make_mantle()
    assert pl.shear_rheology_set is False
    pl.set_shear_rheology(rheo.Maxwell())
    assert pl.shear_rheology_set is True
    assert pl.bulk_rheology_set  is False


def test_attach_shear_rheology_changes_complex_modulus():
    """After attaching Maxwell, calc_complex_shear_modulus matches the model output."""
    rheo = _import_rheology()
    pl   = _make_mantle()
    freq = 1.0e-5
    pl.set_shear_rheology(rheo.Maxwell())
    mu       = pl.calc_complex_shear_modulus(freq)
    expected = rheo.maxwell(freq, _SHEAR_MOD_PA, _VISCOSITY_PAS)
    assert mu.real == pytest.approx(expected.real, rel=1e-9)
    assert mu.imag == pytest.approx(expected.imag, rel=1e-9)
    # A dissipative Maxwell body has a non-zero loss modulus at this frequency.
    assert mu.imag != 0.0


def test_attach_rheology_consumes_wrapper():
    """A rheology model cannot be attached twice (ownership is transferred)."""
    rheo = _import_rheology()
    pl   = _make_mantle()
    model = rheo.Maxwell()
    pl.set_shear_rheology(model)
    with pytest.raises(ValueError):
        pl.set_bulk_rheology(model)


def test_binary_roundtrip_with_rheology():
    """save_binary + load_binary restores attached shear and bulk rheology models."""
    mod  = _import_physics()
    rheo = _import_rheology()
    freq = 1.0e-5

    pl1 = _make_mantle()
    pl1.set_shear_rheology(rheo.Maxwell())
    pl1.set_bulk_rheology(rheo.Andrade(alpha=0.25, zeta=2.0))
    mu_before = pl1.calc_complex_shear_modulus(freq)
    k_before  = pl1.calc_complex_bulk_modulus(freq)

    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        pl1.save_binary(path)
        pl2 = mod.PhysicsLayer("placeholder", 0, 0.0, 1.0, 1.0)
        pl2.load_binary(path)

        assert pl2.shear_rheology_set is True
        assert pl2.bulk_rheology_set  is True
        mu_after = pl2.calc_complex_shear_modulus(freq)
        k_after  = pl2.calc_complex_bulk_modulus(freq)
        assert mu_after.real == pytest.approx(mu_before.real, rel=1e-12)
        assert mu_after.imag == pytest.approx(mu_before.imag, rel=1e-12)
        assert k_after.real  == pytest.approx(k_before.real,  rel=1e-12)
        assert k_after.imag  == pytest.approx(k_before.imag,  rel=1e-12)
    finally:
        os.unlink(path)


def test_binary_roundtrip_without_rheology_flags_false():
    """A layer with no rheology round-trips with both flags still False."""
    mod = _import_physics()
    pl1 = _make_mantle()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        pl1.save_binary(path)
        pl2 = mod.PhysicsLayer("placeholder", 0, 0.0, 1.0, 1.0)
        pl2.load_binary(path)
        assert pl2.shear_rheology_set is False
        assert pl2.bulk_rheology_set  is False
    finally:
        os.unlink(path)


def test_binary_roundtrip_only_shear_rheology():
    """Attaching only the shear model leaves the bulk slot absent after a round-trip."""
    mod  = _import_physics()
    rheo = _import_rheology()
    pl1  = _make_mantle()
    pl1.set_shear_rheology(rheo.Maxwell())
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        pl1.save_binary(path)
        pl2 = mod.PhysicsLayer("placeholder", 0, 0.0, 1.0, 1.0)
        pl2.load_binary(path)
        assert pl2.shear_rheology_set is True
        assert pl2.bulk_rheology_set  is False
    finally:
        os.unlink(path)


# =====================================================================================================================
# isinstance checks
# =====================================================================================================================

def test_physics_layer_is_base_layer():
    """PhysicsLayer is an instance of BaseLayer."""
    from TidalPy.structures_x.layers.base import BaseLayer
    assert isinstance(_make_mantle(), BaseLayer)


def test_physics_layer_is_structure_base():
    """PhysicsLayer is an instance of StructureBase."""
    from TidalPy.Utilities_x.classes_x.classes import StructureBase
    assert isinstance(_make_mantle(), StructureBase)


def test_physics_layer_is_tidalpy_base():
    """PhysicsLayer is an instance of TidalPyBaseClass."""
    from TidalPy.Utilities_x.classes_x.classes import TidalPyBaseClass
    assert isinstance(_make_mantle(), TidalPyBaseClass)
