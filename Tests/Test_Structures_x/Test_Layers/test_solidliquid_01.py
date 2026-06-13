"""
Tests for TidalPy.structures_x.layers.solidliquid — SolidLiquidLayer (Phase 3).

Covers construction, geometry/EOS/tidal inheritance (via _layer_ptr), thermal
property getters, melt fraction, Arrhenius viscosity, dynamic shear modulus,
thermal conductivity/diffusivity, adiabatic gradient, conductive heat flux,
radiogenic heating (no sub-model), binary round-trip, TOML config save,
and isinstance checks.

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
_T_SOLIDUS_K    = 3000.0     # [K]
_T_LIQUIDUS_K   = 4000.0     # [K]
_T_REF_K        = 3000.0     # [K]
_K_COND         = 4.5        # [W/(m·K)]
_ALPHA          = 2.0e-5     # [1/K]
_CP             = 1200.0     # [J/(kg·K)]
_E_A            = 300.0e3    # [J/mol]
_V_A            = 5.0e-6     # [m³/mol]
_RHO_REF        = 4000.0     # [kg/m³]
_MELT_EXP       = 1.0
_MELT_RED       = 25.0
_R_GAS          = 8.314462618


def _make_layer(**kw):
    mod = _import_solidliquid()
    defaults = dict(
        name                          = "mantle",
        layer_index                   = 1,
        radius_inner_m                = _R_INNER_M,
        radius_outer_m                = _R_OUTER_M,
        mass_kg                       = _MASS_KG,
        material_name                 = "perovskite",
        is_tidal                      = True,
        tidal_scale                   = 1.0,
        shear_modulus_static_pa       = _SHEAR_PA,
        bulk_modulus_static_pa        = 3.57e11,
        shear_viscosity_static_pas    = _SHEAR_VISC_PAS,
        bulk_viscosity_static_pas     = _BULK_VISC_PAS,
        thermal_conductivity_ref_w_mk = _K_COND,
        thermal_expansion_ref_1_k     = _ALPHA,
        heat_capacity_ref_j_kgk       = _CP,
        activation_energy_j_mol       = _E_A,
        activation_volume_m3_mol      = _V_A,
        solidus_temperature_k         = _T_SOLIDUS_K,
        liquidus_temperature_k        = _T_LIQUIDUS_K,
        melt_fraction_exponent        = _MELT_EXP,
        reference_density_kg_m3       = _RHO_REF,
        reference_temperature_k       = _T_REF_K,
        melt_viscosity_reduction      = _MELT_RED,
    )
    defaults.update(kw)
    return mod.SolidLiquidLayer(**defaults)


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
    assert sl.thermal_conductivity_ref     == pytest.approx(_K_COND)
    assert sl.thermal_expansion_ref        == pytest.approx(_ALPHA)
    assert sl.heat_capacity_ref            == pytest.approx(_CP)
    assert sl.activation_energy            == pytest.approx(_E_A)
    assert sl.activation_volume            == pytest.approx(_V_A)
    assert sl.solidus_temperature          == pytest.approx(_T_SOLIDUS_K)
    assert sl.liquidus_temperature         == pytest.approx(_T_LIQUIDUS_K)
    assert sl.melt_fraction_exponent       == pytest.approx(_MELT_EXP)
    assert sl.reference_density            == pytest.approx(_RHO_REF)
    assert sl.reference_temperature        == pytest.approx(_T_REF_K)
    assert sl.melt_viscosity_reduction     == pytest.approx(_MELT_RED)


def test_solidliquid_defaults():
    """SolidLiquidLayer optional parameters default to expected values."""
    mod = _import_solidliquid()
    sl  = mod.SolidLiquidLayer("test", 0, 0.0, 1e6, 1e20)
    assert sl.thermal_conductivity_ref  == pytest.approx(4.0)
    assert sl.thermal_expansion_ref     == pytest.approx(3.0e-5)
    assert sl.heat_capacity_ref         == pytest.approx(1200.0)
    assert sl.solidus_temperature       == pytest.approx(1600.0)
    assert sl.liquidus_temperature      == pytest.approx(2000.0)
    assert sl.melt_fraction_exponent    == pytest.approx(1.0)
    assert sl.reference_density         == pytest.approx(3500.0)
    assert sl.melt_viscosity_reduction  == pytest.approx(25.0)
    assert sl.shear_viscosity_static    == pytest.approx(0.0)
    assert sl.bulk_viscosity_static     == pytest.approx(0.0)


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


def test_solidliquid_inherits_tidal_susceptibility():
    """calc_tidal_susceptibility inherited from PhysicsLayer works correctly."""
    sl = _make_layer()
    G  = 6.674e-11
    expected = (1.5 * _R_OUTER_M**5) / (G * _MASS_KG**2)
    assert sl.calc_tidal_susceptibility() == pytest.approx(expected, rel=1e-3)


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
def test_melt_fraction_below_solidus():
    """Melt fraction is 0 below solidus."""
    sl = _make_layer()
    assert sl.calc_melt_fraction(_T_SOLIDUS_K - 100.0) == pytest.approx(0.0)


def test_melt_fraction_above_liquidus():
    """Melt fraction is 1 above liquidus."""
    sl = _make_layer()
    assert sl.calc_melt_fraction(_T_LIQUIDUS_K + 100.0) == pytest.approx(1.0)


def test_melt_fraction_at_solidus():
    """Melt fraction is 0 at solidus."""
    sl = _make_layer()
    assert sl.calc_melt_fraction(_T_SOLIDUS_K) == pytest.approx(0.0)


def test_melt_fraction_at_liquidus():
    """Melt fraction is 1 at liquidus."""
    sl = _make_layer()
    assert sl.calc_melt_fraction(_T_LIQUIDUS_K) == pytest.approx(1.0)


def test_melt_fraction_midpoint_linear():
    """Melt fraction is 0.5 at mid-point between solidus and liquidus (exponent=1)."""
    sl   = _make_layer(melt_fraction_exponent=1.0)
    T_mid = (_T_SOLIDUS_K + _T_LIQUIDUS_K) / 2.0
    assert sl.calc_melt_fraction(T_mid) == pytest.approx(0.5)


def test_melt_fraction_exponent_gt1():
    """Melt fraction < 0.5 at midpoint when exponent > 1."""
    sl    = _make_layer(melt_fraction_exponent=2.0)
    T_mid = (_T_SOLIDUS_K + _T_LIQUIDUS_K) / 2.0
    assert sl.calc_melt_fraction(T_mid) == pytest.approx(0.25)


# =====================================================================================================================
# Viscosity
# =====================================================================================================================
def test_viscosity_at_reference_zero_pressure():
    """Viscosity equals shear reference value at reference T, P=0, below solidus."""
    sl  = _make_layer(solidus_temperature_k=1e10)   # push solidus far so phi = 0
    eta = sl.calc_viscosity(_T_REF_K, 0.0)
    assert eta == pytest.approx(_SHEAR_VISC_PAS, rel=1e-9)


def test_viscosity_increases_below_reference_temp():
    """Viscosity is higher at T < T_ref (Arrhenius — colder = more viscous)."""
    sl   = _make_layer(solidus_temperature_k=1e10)
    eta_cold = sl.calc_viscosity(_T_REF_K * 0.5, 0.0)
    eta_ref  = sl.calc_viscosity(_T_REF_K, 0.0)
    assert eta_cold > eta_ref


def test_viscosity_decreases_above_reference_temp():
    """Viscosity is lower at T > T_ref."""
    sl   = _make_layer(solidus_temperature_k=1e10)
    eta_hot = sl.calc_viscosity(_T_REF_K * 2.0, 0.0)
    eta_ref = sl.calc_viscosity(_T_REF_K, 0.0)
    assert eta_hot < eta_ref


def test_viscosity_melt_reduction():
    """Partial melt reduces viscosity."""
    sl = _make_layer()
    T_partial = (_T_SOLIDUS_K + _T_LIQUIDUS_K) / 2.0   # phi = 0.5
    eta_melt  = sl.calc_viscosity(T_partial, 0.0)
    # With phi = 0.5, melt_viscosity_reduction = 25: factor = exp(-12.5)
    T_solid   = _T_SOLIDUS_K - 1.0                      # phi ≈ 0
    eta_solid = sl.calc_viscosity(T_solid, 0.0)
    assert eta_melt < eta_solid


def test_viscosity_zero_temp_guard():
    """calc_viscosity returns reference viscosity rather than crashing at T=0."""
    sl  = _make_layer()
    eta = sl.calc_viscosity(0.0, 0.0)
    assert math.isfinite(eta)


# =====================================================================================================================
# Dynamic shear modulus
# =====================================================================================================================

def test_shear_modulus_below_solidus_equals_static():
    """calc_shear_modulus equals static value when fully solid (phi=0)."""
    sl = _make_layer()
    G  = sl.calc_shear_modulus(_T_SOLIDUS_K - 1.0, 0.0)
    assert G == pytest.approx(_SHEAR_PA)


def test_shear_modulus_above_liquidus_is_zero():
    """calc_shear_modulus is 0 when fully molten (phi=1)."""
    sl = _make_layer()
    G  = sl.calc_shear_modulus(_T_LIQUIDUS_K + 1.0, 0.0)
    assert G == pytest.approx(0.0)


def test_shear_modulus_midpoint():
    """calc_shear_modulus is G_static/2 at midpoint (linear melt, phi=0.5)."""
    sl    = _make_layer(melt_fraction_exponent=1.0)
    T_mid = (_T_SOLIDUS_K + _T_LIQUIDUS_K) / 2.0
    G     = sl.calc_shear_modulus(T_mid, 0.0)
    assert G == pytest.approx(_SHEAR_PA * 0.5)


# =====================================================================================================================
# Thermal conductivity and diffusivity
# =====================================================================================================================

def test_thermal_conductivity_constant():
    """calc_thermal_conductivity returns reference value."""
    sl = _make_layer()
    assert sl.calc_thermal_conductivity(1000.0) == pytest.approx(_K_COND)
    assert sl.calc_thermal_conductivity(4000.0) == pytest.approx(_K_COND)


def test_thermal_diffusivity_formula():
    """calc_thermal_diffusivity = k / (rho_ref * c_p)."""
    sl       = _make_layer()
    expected = _K_COND / (_RHO_REF * _CP)
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
    "shear_modulus_static_pa", "bulk_modulus_static_pa",
    "shear_viscosity_static_pas", "bulk_viscosity_static_pas",
    "love_number_k_re", "love_number_k_im",
    "love_number_h_re", "love_number_h_im",
    "love_number_l_re", "love_number_l_im",
    "thermal_conductivity_ref_w_mk", "thermal_expansion_ref_1_k",
    "heat_capacity_ref_j_kgk", "activation_energy_j_mol", "activation_volume_m3_mol",
    "solidus_temperature_k", "liquidus_temperature_k", "melt_fraction_exponent",
    "reference_density_kg_m3", "reference_temperature_k", "melt_viscosity_reduction",
)


def test_get_config_dict_has_all_keys():
    """get_config_dict contains all expected keys."""
    sl  = _make_layer()
    cfg = sl.get_config_dict()
    for key in _ALL_KEYS:
        assert key in cfg, f"Missing key: {key}"


def test_get_config_dict_values():
    """get_config_dict values match constructor arguments."""
    sl  = _make_layer()
    cfg = sl.get_config_dict()
    assert cfg["name"]                          == "mantle"
    assert cfg["layer_index"]                   == 1
    assert cfg["shear_modulus_static_pa"]       == pytest.approx(_SHEAR_PA)
    assert cfg["shear_viscosity_static_pas"]    == pytest.approx(_SHEAR_VISC_PAS)
    assert cfg["bulk_viscosity_static_pas"]     == pytest.approx(_BULK_VISC_PAS)
    assert cfg["thermal_conductivity_ref_w_mk"] == pytest.approx(_K_COND)
    assert cfg["solidus_temperature_k"]         == pytest.approx(_T_SOLIDUS_K)
    assert cfg["liquidus_temperature_k"]        == pytest.approx(_T_LIQUIDUS_K)
    assert cfg["reference_temperature_k"]       == pytest.approx(_T_REF_K)
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
        assert data["thermal_conductivity_ref_w_mk"] == pytest.approx(_K_COND)
        assert data["solidus_temperature_k"]         == pytest.approx(_T_SOLIDUS_K)
        assert data["liquidus_temperature_k"]        == pytest.approx(_T_LIQUIDUS_K)
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
        assert sl2.thermal_conductivity_ref     == pytest.approx(_K_COND)
        assert sl2.solidus_temperature          == pytest.approx(_T_SOLIDUS_K)
        assert sl2.liquidus_temperature         == pytest.approx(_T_LIQUIDUS_K)
        assert sl2.reference_temperature        == pytest.approx(_T_REF_K)
        assert sl2.melt_viscosity_reduction     == pytest.approx(_MELT_RED)
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


def test_binary_roundtrip_calc_preserves_modulus():
    """After load_binary, calc_shear_modulus works with restored static value."""
    mod = _import_solidliquid()
    sl1 = _make_layer()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        sl1.save_binary(path)
        sl2 = mod.SolidLiquidLayer("placeholder", 0, 0.0, 1.0, 1.0)
        sl2.load_binary(path)
        G = sl2.calc_shear_modulus(_T_SOLIDUS_K - 1.0, 0.0)
        assert G == pytest.approx(_SHEAR_PA)
    finally:
        os.unlink(path)


def test_binary_load_file_not_found():
    """load_binary raises FileNotFoundError for a missing path."""
    sl = _make_layer()
    with pytest.raises(FileNotFoundError):
        sl.load_binary("/nonexistent/path/xyz.tpyb")


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
