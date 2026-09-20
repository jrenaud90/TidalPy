"""
Tests for TidalPy.structures_x.layers.solidliquid — SolidLiquidLayer.

Covers construction, geometry/EOS/tidal inheritance (via _layer_ptr), thermal
property getters, thermal conductivity/diffusivity, adiabatic gradient,
conductive heat flux, radiogenic heating (no sub-model), binary round-trip,
TOML config save, and isinstance checks.

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

def _import_solidliquid():
    try:
        from TidalPy.structures_x.layers import solidliquid as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.structures_x.layers.solidliquid not compiled — run uv pip install first."
        )


# Reference values (Earth lower mantle, MKS)
_R_INNER_M      = 3.485e6    # [m]  CMB radius
_R_OUTER_M      = 6.371e6    # [m]  surface radius
_MASS_KG        = 4.043e24   # [kg]
_SHEAR_PA       = 1.67e11    # [Pa]
_SHEAR_VISC_PAS = 1.0e21     # [Pa·s]
_BULK_VISC_PAS  = 2.0e21     # [Pa·s]
_K_COND         = 4.5        # [W/(m·K)]
_ALPHA          = 2.0e-5     # [1/K]
_CP             = 1200.0     # [J/(kg·K)]
_RHO_EOS        = 4000.0     # [kg/m³] the material's density


# Everything about the material goes to the EOS model the layer is given: the static constants and the thermal
# constants (conductivity, expansivity, heat capacity) alike.
_MATERIAL_KEYS = ("shear_modulus_static", "bulk_modulus_static", "shear_viscosity_static",
                  "bulk_viscosity_static", "thermal_conductivity", "thermal_expansion", "heat_capacity",
                  "reference_density")


def _make_layer(**kw):
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    mod = _import_solidliquid()
    defaults = dict(
        name                 = "mantle",
        layer_index          = 1,
        radius_inner         = _R_INNER_M,
        radius_outer         = _R_OUTER_M,
        mass                 = _MASS_KG,
        material_name        = "perovskite",
        is_tidal             = True,
        tidal_scale          = 1.0,
        shear_modulus_static = _SHEAR_PA,
        bulk_modulus_static  = 3.57e11,
        shear_viscosity_static = _SHEAR_VISC_PAS,
        bulk_viscosity_static  = _BULK_VISC_PAS,
        thermal_conductivity = _K_COND,
        thermal_expansion    = _ALPHA,
        heat_capacity        = _CP,
        reference_density    = _RHO_EOS,
    )
    defaults.update(kw)
    material = {key: defaults.pop(key) for key in _MATERIAL_KEYS if key in defaults}
    layer = mod.SolidLiquidLayer(**defaults)
    layer.set_eos(ConstantDensityEOS(**material))
    return layer


# =====================================================================================================================
# Construction
# =====================================================================================================================
def test_solidliquid_construction_basic():
    """SolidLiquidLayer stores all config values at construction."""
    sl = _make_layer()
    assert sl.name                         == "mantle"
    assert sl.layer_index                  == 1
    assert sl.shear_modulus_static         == pytest.approx(_SHEAR_PA)
    assert sl.shear_viscosity_static       == pytest.approx(_SHEAR_VISC_PAS)
    assert sl.bulk_viscosity_static        == pytest.approx(_BULK_VISC_PAS)
    # The thermal constants are read from the layer's material.
    assert sl.thermal_conductivity         == pytest.approx(_K_COND)
    assert sl.thermal_expansion            == pytest.approx(_ALPHA)
    assert sl.heat_capacity                == pytest.approx(_CP)


def test_solidliquid_defaults():
    """Without a material the thermal constants are NaN; a default material supplies k = 4, alpha = 0, c_p = 1200."""
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    mod = _import_solidliquid()
    sl  = mod.SolidLiquidLayer("test", 0, 0.0, 1e6, 1e20)
    assert math.isnan(sl.thermal_conductivity)
    assert math.isnan(sl.shear_viscosity_static)
    sl.set_eos(ConstantDensityEOS())
    assert sl.thermal_conductivity      == pytest.approx(4.0)
    assert sl.thermal_expansion         == 0.0          # one alpha: zero keeps the density law athermal
    assert sl.heat_capacity             == pytest.approx(1200.0)
    assert math.isnan(sl.shear_viscosity_static)
    assert math.isnan(sl.bulk_viscosity_static)
    assert (sl.is_solid, sl.is_static, sl.is_incompressible) == (True, True, False)


def test_solidliquid_layer_assumption_flags_from_constructor():
    sl = _make_layer(is_solid=False, is_static=False, is_incompressible=True)
    assert (sl.is_solid, sl.is_static, sl.is_incompressible) == (False, False, True)
    cfg = sl.get_config_dict()
    assert (cfg["is_solid"], cfg["is_static"], cfg["is_incompressible"]) == (False, False, True)


# =====================================================================================================================
# Inheritance checks (via _layer_ptr — critical regression)
# =====================================================================================================================
def test_solidliquid_inherits_geometry():
    """Thickness, volume, surface areas resolved correctly via _layer_ptr."""
    sl = _make_layer()
    expected_t  = _R_OUTER_M - _R_INNER_M
    expected_v  = (4.0 / 3.0) * math.pi * (_R_OUTER_M**3 - _R_INNER_M**3)
    assert sl.thickness        == pytest.approx(expected_t)
    assert sl.volume           == pytest.approx(expected_v, rel=1e-9)
    assert sl.radius           == pytest.approx(_R_OUTER_M)
    assert sl.radius_inner     == pytest.approx(_R_INNER_M)
    assert sl.mass             == pytest.approx(_MASS_KG)


def test_solidliquid_inherits_eos():
    """EOS operations inherited from BaseLayer work on SolidLiquidLayer."""
    sl = _make_layer()
    assert sl.eos_data_populated is False
    assert math.isnan(sl.get_density(_R_INNER_M))
    sl.update_eos_data(
        [_R_INNER_M, _R_OUTER_M],
        [5000.0, 4000.0],
        [10.0,   9.81],
        [1e11,   0.0],
    )
    assert sl.eos_data_populated is True
    assert sl.get_density(_R_INNER_M) == pytest.approx(5000.0)


def test_solidliquid_inherits_love_numbers():
    """love_numbers property inherited from PhysicsLayer returns LoveNumbers object."""
    mod = _import_solidliquid()
    sl  = mod.SolidLiquidLayer("test", 0, 0.0, 1e6, 1e20,
                               love_number_k=0.3-0.01j,
                               love_number_h=0.6-0.02j,
                               love_number_l=0.1-0.005j)
    assert sl.love_number_k == pytest.approx(0.3 - 0.01j)
    assert sl.love_number_h == pytest.approx(0.6 - 0.02j)
    assert sl.love_number_l == pytest.approx(0.1 - 0.005j)
    ln = sl.love_numbers
    assert ln.k == pytest.approx(0.3 - 0.01j)
    assert ln.h == pytest.approx(0.6 - 0.02j)
    assert ln.l == pytest.approx(0.1 - 0.005j)


# =====================================================================================================================
# Melt fraction
# =====================================================================================================================
def test_thermal_conductivity_constant():
    """calc_thermal_conductivity returns reference value."""
    sl = _make_layer()
    assert sl.calc_thermal_conductivity(1000.0) == pytest.approx(_K_COND)
    assert sl.calc_thermal_conductivity(4000.0) == pytest.approx(_K_COND)


def test_thermal_diffusivity_formula():
    """calc_thermal_diffusivity = k / (rho c_p), with the layer's bulk density (its mass over its volume)."""
    sl       = _make_layer()
    expected = _K_COND / (sl.density_bulk * _CP)
    assert sl.calc_thermal_diffusivity(2000.0) == pytest.approx(expected)


# =====================================================================================================================
# Adiabatic gradient (no EOS)
# =====================================================================================================================
def test_adiabatic_gradient_no_eos_returns_zero():
    """calc_adiabatic_temperature_gradient returns 0 when EOS not populated."""
    sl = _make_layer()
    assert sl.calc_adiabatic_temperature_gradient(2000.0) == pytest.approx(0.0)


# =====================================================================================================================
# Conductive heat flux
# =====================================================================================================================
def test_heat_flux_conductive_formula():
    """calc_heat_flux_conductive = k * dT / thickness."""
    sl        = _make_layer()
    T_base    = 3000.0
    T_top     = 1000.0
    thickness = _R_OUTER_M - _R_INNER_M
    expected  = _K_COND * (T_base - T_top) / thickness
    assert sl.calc_heat_flux_conductive(T_base, T_top) == pytest.approx(expected)


def test_heat_flux_conductive_zero_when_isothermal():
    """calc_heat_flux_conductive is 0 when T_base == T_top."""
    sl = _make_layer()
    assert sl.calc_heat_flux_conductive(2000.0, 2000.0) == pytest.approx(0.0)


def test_heat_flux_conductive_zero_thickness():
    """calc_heat_flux_conductive returns 0 for zero-thickness layer."""
    mod = _import_solidliquid()
    sl  = mod.SolidLiquidLayer("thin", 0, 1e6, 1e6, 1e20)
    assert sl.calc_heat_flux_conductive(3000.0, 1000.0) == pytest.approx(0.0)


# =====================================================================================================================
# Radiogenic heating (no sub-model)
# =====================================================================================================================
def test_radiogenic_heating_no_submodel_is_zero():
    """calc_radiogenic_heating returns 0 before a sub-model is attached."""
    sl = _make_layer()
    assert sl.calc_radiogenic_heating(0.0, _MASS_KG) == pytest.approx(0.0)


def test_radiogenics_not_set_initially():
    """radiogenics_set is False on a new layer."""
    sl = _make_layer()
    assert sl.radiogenics_set is False


def test_cooling_not_set_initially():
    """cooling_set is False on a new layer."""
    sl = _make_layer()
    assert sl.cooling_set is False


# =====================================================================================================================
# get_config_dict
# =====================================================================================================================
_ALL_KEYS = (
    "name", "layer_index", "radius_inner_m", "radius_outer_m", "mass_kg",
    "material_name", "is_tidal", "tidal_scale",
    "love_number_k_re", "love_number_k_im",
    "love_number_h_re", "love_number_h_im",
    "love_number_l_re", "love_number_l_im",
)
# The thermal constants sit in the material table.
_MATERIAL_TABLE_KEYS = ("thermal_conductivity_w_mk", "thermal_expansion_1_k", "heat_capacity_j_kgk")


def test_get_config_dict_has_all_keys():
    """get_config_dict contains all expected keys."""
    sl  = _make_layer()
    cfg = sl.get_config_dict()
    for key in _ALL_KEYS:
        assert key in cfg, f"Missing key: {key}"
    for key in _MATERIAL_TABLE_KEYS:
        assert key in cfg["material"], f"Missing material key: {key}"


def test_get_config_dict_values():
    """get_config_dict values match constructor arguments."""
    sl  = _make_layer()
    cfg = sl.get_config_dict()
    assert cfg["name"]                          == "mantle"
    assert cfg["layer_index"]                   == 1
    assert cfg["material"]["shear_modulus_static_pa"]       == pytest.approx(_SHEAR_PA)
    assert cfg["material"]["shear_viscosity_static_pas"]    == pytest.approx(_SHEAR_VISC_PAS)
    assert cfg["material"]["bulk_viscosity_static_pas"]     == pytest.approx(_BULK_VISC_PAS)
    assert cfg["material"]["thermal_conductivity_w_mk"] == pytest.approx(_K_COND)
    assert cfg["material"]["thermal_expansion_1_k"]     == pytest.approx(_ALPHA)
    assert cfg["material"]["heat_capacity_j_kgk"]       == pytest.approx(_CP)
    assert cfg["material"]["reference_density_kg_m3"]   == pytest.approx(_RHO_EOS)
    for moved in ("thermal_conductivity_ref_w_mk", "reference_density_kg_m3", "reference_temperature_k"):
        assert moved not in cfg
    assert cfg["love_number_k_re"]              == pytest.approx(0.0)
    assert cfg["love_number_k_im"]              == pytest.approx(0.0)


# =====================================================================================================================
# save_config (TOML)
# =====================================================================================================================
def test_save_config_solidliquid():
    """save_config writes a TOML file with all SolidLiquidLayer keys."""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            pytest.skip("Neither tomllib nor tomli available.")

    sl = _make_layer()
    with tempfile.NamedTemporaryFile(suffix=".toml", delete=False, mode="w") as f:
        path = f.name
    try:
        sl.save_config(path)
        with open(path, "rb") as f:
            data = tomllib.load(f)
        assert data["name"]                          == "mantle"
        assert data["material"]["thermal_conductivity_w_mk"] == pytest.approx(_K_COND)
        assert data["material"]["heat_capacity_j_kgk"]       == pytest.approx(_CP)
    finally:
        os.unlink(path)


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
def test_binary_roundtrip_solidliquid():
    """save_binary + load_binary preserves all fields."""
    mod = _import_solidliquid()
    sl1 = _make_layer()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        sl1.save_binary(path)
        sl2 = mod.SolidLiquidLayer("placeholder", 0, 0.0, 1.0, 1.0)
        sl2.load_binary(path)
        assert sl2.name                         == "mantle"
        assert sl2.layer_index                  == 1
        assert sl2.radius_outer                 == pytest.approx(_R_OUTER_M)
        assert sl2.radius_inner                 == pytest.approx(_R_INNER_M)
        assert sl2.mass                         == pytest.approx(_MASS_KG)
        assert sl2.shear_modulus_static         == pytest.approx(_SHEAR_PA)
        assert sl2.shear_viscosity_static       == pytest.approx(_SHEAR_VISC_PAS)
        assert sl2.bulk_viscosity_static        == pytest.approx(_BULK_VISC_PAS)
        assert sl2.thermal_conductivity         == pytest.approx(_K_COND)
        assert sl2.thermal_expansion            == pytest.approx(_ALPHA)
        assert sl2.heat_capacity                == pytest.approx(_CP)
    finally:
        os.unlink(path)


def test_binary_roundtrip_derived_fields():
    """After load_binary, derived geometry fields are recomputed correctly."""
    mod = _import_solidliquid()
    sl1 = _make_layer()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        sl1.save_binary(path)
        sl2 = mod.SolidLiquidLayer("placeholder", 0, 0.0, 1.0, 1.0)
        sl2.load_binary(path)
        expected_thickness = _R_OUTER_M - _R_INNER_M
        assert sl2.thickness == pytest.approx(expected_thickness)
    finally:
        os.unlink(path)


def test_binary_load_file_not_found():
    """load_binary raises FileNotFoundError for a missing path."""
    sl = _make_layer()
    with pytest.raises(FileNotFoundError):
        sl.load_binary("/nonexistent/path/xyz.tpyb")


# =====================================================================================================================
# Sub-model attachment and recursive binary serialization
# =====================================================================================================================

def _import_rheology():
    from TidalPy.rheology_x import rheology as _mod
    return _mod


def _import_cooling():
    from TidalPy.cooling_x import cooling as _mod
    return _mod


def _import_radiogenics():
    from TidalPy.radiogenics_x import radiogenics as _mod
    return _mod


def test_attach_submodels_set_flags():
    """set_cooling / set_radiogenics flip their respective flags."""
    cooling = _import_cooling()
    radio   = _import_radiogenics()
    sl = _make_layer()
    assert sl.cooling_set     is False
    assert sl.radiogenics_set is False
    sl.set_cooling(cooling.ConvectiveCooling())
    sl.set_radiogenics(radio.FixedRadiogenics(fixed_heat_production=1.0e-11))
    assert sl.cooling_set     is True
    assert sl.radiogenics_set is True


def test_attach_radiogenics_changes_heating():
    """After attaching a fixed radiogenics model, calc_radiogenic_heating is non-zero."""
    radio = _import_radiogenics()
    sl = _make_layer()
    sl.set_radiogenics(radio.FixedRadiogenics(fixed_heat_production=2.0e-11))
    q = sl.calc_radiogenic_heating(0.0, _MASS_KG)
    assert q == pytest.approx(2.0e-11 * _MASS_KG, rel=1e-12)


def test_binary_roundtrip_with_all_submodels():
    """save_binary + load_binary restores rheology, cooling, and radiogenics models."""
    mod     = _import_solidliquid()
    rheo    = _import_rheology()
    cooling = _import_cooling()
    radio   = _import_radiogenics()
    freq    = 1.0e-5

    sl1 = _make_layer()
    sl1.set_shear_rheology(rheo.Maxwell())
    sl1.set_cooling(cooling.ConvectiveCooling(convection_alpha=0.9, critical_rayleigh=1200.0))
    sl1.set_radiogenics(radio.FixedRadiogenics(fixed_heat_production=3.0e-11))
    mu_before = sl1.calc_complex_shear_modulus(freq)
    q_before  = sl1.calc_radiogenic_heating(0.0, _MASS_KG)

    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        sl1.save_binary(path)
        sl2 = _make_layer(name="placeholder")
        sl2.load_binary(path)

        assert sl2.shear_rheology_set is True
        assert sl2.cooling_set        is True
        assert sl2.radiogenics_set    is True
        mu_after = sl2.calc_complex_shear_modulus(freq)
        q_after  = sl2.calc_radiogenic_heating(0.0, _MASS_KG)
        assert mu_after.real == pytest.approx(mu_before.real, rel=1e-12)
        assert mu_after.imag == pytest.approx(mu_before.imag, rel=1e-12)
        assert q_after       == pytest.approx(q_before,       rel=1e-12)
    finally:
        os.unlink(path)


def test_binary_roundtrip_isotope_radiogenics():
    """Variable-length isotope radiogenics survives a layer binary round-trip."""
    mod   = _import_solidliquid()
    radio = _import_radiogenics()

    sl1 = _make_layer()
    sl1.set_radiogenics(radio.IsotopeRadiogenics.from_dataset("modern_day_chondritic"))
    q_before = sl1.calc_radiogenic_heating(0.0, _MASS_KG)

    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        sl1.save_binary(path)
        sl2 = _make_layer(name="placeholder")
        sl2.load_binary(path)
        assert sl2.radiogenics_set is True
        q_after = sl2.calc_radiogenic_heating(0.0, _MASS_KG)
        assert q_after == pytest.approx(q_before, rel=1e-12)
    finally:
        os.unlink(path)


def test_binary_roundtrip_no_submodels_flags_false():
    """A layer with no sub-models round-trips with all flags False."""
    mod = _import_solidliquid()
    sl1 = _make_layer()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        sl1.save_binary(path)
        sl2 = _make_layer(name="placeholder")
        sl2.load_binary(path)
        assert sl2.shear_rheology_set is False
        assert sl2.bulk_rheology_set  is False
        assert sl2.cooling_set        is False
        assert sl2.radiogenics_set    is False
    finally:
        os.unlink(path)


# =====================================================================================================================
# isinstance checks
# =====================================================================================================================
def test_solidliquid_is_physics_layer():
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    assert isinstance(_make_layer(), PhysicsLayer)


def test_solidliquid_is_base_layer():
    from TidalPy.structures_x.layers.base import BaseLayer
    assert isinstance(_make_layer(), BaseLayer)


def test_solidliquid_is_structure_base():
    from TidalPy.Utilities_x.classes_x.classes import StructureBase
    assert isinstance(_make_layer(), StructureBase)


def test_solidliquid_is_tidalpy_base():
    from TidalPy.Utilities_x.classes_x.classes import TidalPyBaseClass
    assert isinstance(_make_layer(), TidalPyBaseClass)


def test_get_config_dict_class_and_thermal_model_tables():
    """Cooling and radiogenics models appear as sub-tables once attached."""
    from TidalPy.cooling_x.cooling import ConvectiveCooling
    from TidalPy.radiogenics_x.radiogenics import FixedRadiogenics
    sl = _make_layer()
    cfg = sl.get_config_dict()
    assert cfg["class"] == "solidliquid"
    assert "cooling" not in cfg
    assert "radiogenics" not in cfg
    sl.set_cooling(ConvectiveCooling(convection_alpha=0.9))
    sl.set_radiogenics(FixedRadiogenics(fixed_heat_production=3.0e-11))
    cfg = sl.get_config_dict()
    assert cfg["cooling"]["convection_alpha"] == pytest.approx(0.9)
    assert cfg["radiogenics"]["fixed_heat_production_w_kg"] == pytest.approx(3.0e-11)
