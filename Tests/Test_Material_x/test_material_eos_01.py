"""
Tests for TidalPy.Material_x.eos.material_eos — material EOS model hierarchy.

Covers the constant-density, Birch-Murnaghan, Vinet, and interpolated models:
density evaluation, analytic-inversion cross-checks (density-from-pressure round
trips through the forward pressure law), the thermal terms (thermal pressure and
thermal expansion), the bulk modulus, the name factory, config dicts, binary
round-trips, and isinstance checks.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import math
import os
import tempfile

import pytest


def _import_eos():
    from TidalPy.Material_x.eos import material_eos as _mod
    return _mod


# Reference material parameters (silicate-ish, MKS)
_RHO0 = 3500.0      # [kg/m^3]
_K0   = 1.30e11     # [Pa]
_K0P  = 4.5         # [dimensionless]


# =====================================================================================================================
# ConstantDensityEOS
# =====================================================================================================================
def test_constant_density():
    eos = _import_eos().ConstantDensityEOS(reference_density=_RHO0)
    assert eos.reference_density == pytest.approx(_RHO0)
    # Density is independent of pressure / radius.
    assert eos.calc_density(0.0) == pytest.approx(_RHO0)
    assert eos.calc_density(5.0e10) == pytest.approx(_RHO0)
    assert eos.calc_density(0.0, 0.0, 1.0e6) == pytest.approx(_RHO0)


# =====================================================================================================================
# Analytic pressure laws and inversion
# =====================================================================================================================
def test_pressure_laws_zero_at_reference():
    mod = _import_eos()
    # At eta = 1 (rho = rho0) the pressure is exactly zero for both laws.
    assert mod.birch_murnaghan_pressure(1.0, _K0, _K0P) == pytest.approx(0.0, abs=1e-3)
    assert mod.vinet_pressure(1.0, _K0, _K0P) == pytest.approx(0.0, abs=1e-3)


def test_birch_murnaghan_density_at_zero_pressure():
    eos = _import_eos().BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    assert eos.calc_density(0.0) == pytest.approx(_RHO0, rel=1e-9)


@pytest.mark.parametrize("pressure", [1.0e9, 1.0e10, 5.0e10, 1.5e11])
def test_birch_murnaghan_inversion_roundtrip(pressure):
    """density(P) inverted then pushed back through the forward law recovers P."""
    mod = _import_eos()
    eos = mod.BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    rho = eos.calc_density(pressure)
    eta = rho / _RHO0
    assert rho > _RHO0  # compression
    assert mod.birch_murnaghan_pressure(eta, _K0, _K0P) == pytest.approx(pressure, rel=1e-6)


@pytest.mark.parametrize("pressure", [1.0e9, 1.0e10, 5.0e10, 1.5e11])
def test_vinet_inversion_roundtrip(pressure):
    mod = _import_eos()
    eos = mod.VinetEOS(_RHO0, _K0, _K0P)
    rho = eos.calc_density(pressure)
    eta = rho / _RHO0
    assert rho > _RHO0
    assert mod.vinet_pressure(eta, _K0, _K0P) == pytest.approx(pressure, rel=1e-6)


def test_density_monotonic_in_pressure():
    mod = _import_eos()
    for cls in (mod.BirchMurnaghanEOS, mod.VinetEOS):
        eos = cls(_RHO0, _K0, _K0P)
        rhos = [eos.calc_density(p) for p in (0.0, 1e10, 5e10, 1e11)]
        assert all(b > a for a, b in zip(rhos, rhos[1:]))


def test_bm_and_vinet_agree_at_small_compression():
    """BM and Vinet should give similar densities at modest pressure."""
    mod = _import_eos()
    bm = mod.BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    vi = mod.VinetEOS(_RHO0, _K0, _K0P)
    p = 5.0e9
    assert bm.calc_density(p) == pytest.approx(vi.calc_density(p), rel=0.02)


# =====================================================================================================================
# InterpolatedEOS
# =====================================================================================================================
def test_interpolated_density():
    eos = _import_eos().InterpolatedEOS([0.0, 1.0e6, 2.0e6], [5000.0, 4000.0, 3000.0])
    assert eos.num_points == 3
    # Exact node values.
    assert eos.calc_density(0.0, 0.0, 0.0)     == pytest.approx(5000.0)
    assert eos.calc_density(0.0, 0.0, 2.0e6)   == pytest.approx(3000.0)
    # Midpoint linear interpolation.
    assert eos.calc_density(0.0, 0.0, 0.5e6)   == pytest.approx(4500.0)
    # Boundary clamping.
    assert eos.calc_density(0.0, 0.0, -1.0e6)  == pytest.approx(5000.0)
    assert eos.calc_density(0.0, 0.0, 9.0e6)   == pytest.approx(3000.0)


def test_interpolated_length_mismatch_raises():
    with pytest.raises(ValueError):
        _import_eos().InterpolatedEOS([0.0, 1.0e6], [5000.0])


# =====================================================================================================================
# Factory
# =====================================================================================================================
@pytest.mark.parametrize("name,cls_attr", [
    ("constant", "ConstantDensityEOS"),
    ("uniform", "ConstantDensityEOS"),
    ("bm", "BirchMurnaghanEOS"),
    ("birch_murnaghan", "BirchMurnaghanEOS"),
    ("vinet", "VinetEOS"),
    ("interp", "InterpolatedEOS"),
])
def test_factory_aliases(name, cls_attr):
    mod = _import_eos()
    cfg = {"reference_density_kg_m3": _RHO0, "reference_bulk_modulus_pa": _K0,
           "bulk_modulus_derivative": _K0P, "radius_m": [0.0, 1.0e6],
           "density_kg_m3": [5000.0, 4000.0]}
    eos = mod.make_material_eos(name, cfg)
    assert isinstance(eos, getattr(mod, cls_attr))


def test_factory_unknown_name_raises():
    with pytest.raises(ValueError):
        _import_eos().make_material_eos("not_a_model")


def test_factory_returns_usable_model():
    mod = _import_eos()
    eos = mod.make_material_eos("bm", {"reference_density_kg_m3": _RHO0,
                                       "reference_bulk_modulus_pa": _K0,
                                       "bulk_modulus_derivative": _K0P})
    assert eos.calc_density(0.0) == pytest.approx(_RHO0, rel=1e-9)
    assert eos.reference_bulk_modulus == pytest.approx(_K0)


# =====================================================================================================================
# Config dict
# =====================================================================================================================
def test_config_dict_bm():
    eos = _import_eos().BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    cfg = eos.get_config_dict()
    assert cfg["model"] == "birch_murnaghan"
    assert cfg["reference_density_kg_m3"]   == pytest.approx(_RHO0)
    assert cfg["reference_bulk_modulus_pa"] == pytest.approx(_K0)
    assert cfg["bulk_modulus_derivative"]   == pytest.approx(_K0P)
    # The inversion settings are config-driven (carry their C++ defaults here).
    assert "invert_rtol" in cfg
    assert "invert_max_iters" in cfg
    assert cfg["invert_rtol"] > 0.0
    assert cfg["invert_max_iters"] > 0


# =====================================================================================================================
# Configurable inversion settings
# =====================================================================================================================
@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
def test_invert_settings_configurable(cls_name):
    mod = _import_eos()
    cls = getattr(mod, cls_name)
    # Defaults present and positive.
    default = cls(_RHO0, _K0, _K0P)
    assert default.invert_rtol > 0.0
    assert default.invert_max_iters > 0
    # Caller-supplied overrides are stored.
    tuned = cls(_RHO0, _K0, _K0P, invert_rtol=1.0e-8, invert_max_iters=80)
    assert tuned.invert_rtol == pytest.approx(1.0e-8)
    assert tuned.invert_max_iters == 80
    # A looser-but-still-tight tolerance gives effectively the same density.
    p = 5.0e10
    assert tuned.calc_density(p) == pytest.approx(default.calc_density(p), rel=1e-6)


def test_invert_settings_via_factory_and_binary():
    mod = _import_eos()
    eos = mod.make_material_eos("bm", {"reference_density_kg_m3": _RHO0,
                                       "reference_bulk_modulus_pa": _K0,
                                       "bulk_modulus_derivative": _K0P,
                                       "invert_rtol": 1.0e-9,
                                       "invert_max_iters": 75})
    assert eos.invert_rtol == pytest.approx(1.0e-9)
    assert eos.invert_max_iters == 75
    # Inversion settings survive a binary round-trip.
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        eos.save_binary(path)
        eos2 = mod.make_material_eos("bm")
        eos2.load_binary(path)
        assert eos2.invert_rtol == pytest.approx(1.0e-9)
        assert eos2.invert_max_iters == 75
    finally:
        os.unlink(path)


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
@pytest.mark.parametrize("factory", [
    lambda m: m.ConstantDensityEOS(_RHO0),
    lambda m: m.BirchMurnaghanEOS(_RHO0, _K0, _K0P),
    lambda m: m.VinetEOS(_RHO0, _K0, _K0P),
    lambda m: m.InterpolatedEOS([0.0, 1.0e6, 2.0e6], [5000.0, 4000.0, 3000.0]),
])
def test_binary_roundtrip(factory):
    mod = _import_eos()
    eos = factory(mod)
    p_test = 3.0e10
    rho_before = eos.calc_density(p_test, 0.0, 0.5e6)
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        eos.save_binary(path)
        eos2 = mod.make_material_eos(eos.model_name)
        eos2.load_binary(path)
        assert eos2.model_name == eos.model_name
        assert eos2.calc_density(p_test, 0.0, 0.5e6) == pytest.approx(rho_before, rel=1e-12)
    finally:
        os.unlink(path)


# =====================================================================================================================
# isinstance
# =====================================================================================================================
def test_isinstance_chain():
    from TidalPy.Utilities_x.classes_x.classes import PhysicsBase, TidalPyBaseClass
    mod = _import_eos()
    eos = mod.BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    assert isinstance(eos, mod.MaterialEOSBase)
    assert isinstance(eos, PhysicsBase)
    assert isinstance(eos, TidalPyBaseClass)


def test_config_dict_interpolated_roundtrip():
    """An interpolated EOS emits its tables (optional ones only when supplied) and rebuilds through the factory."""
    mod = _import_eos()
    radii = [0.0, 1.0e6, 2.0e6]
    densities = [5000.0, 4000.0, 3000.0]
    eos = mod.InterpolatedEOS(radii, densities)
    cfg = eos.get_config_dict()
    assert cfg["model"] == eos.model_name
    assert cfg["radius_m"] == pytest.approx(radii)
    assert cfg["density_kg_m3"] == pytest.approx(densities)
    for optional_key in ("shear_modulus_pa", "bulk_modulus_pa", "shear_viscosity_pas", "bulk_viscosity_pas"):
        assert optional_key not in cfg
    rebuilt = mod.make_material_eos(cfg["model"], {key: value for key, value in cfg.items() if key != "model"})
    assert rebuilt.num_points == 3
    assert rebuilt.get_config_dict() == cfg


# =====================================================================================================================
# Thermal terms
# =====================================================================================================================
_ALPHA = 3.0e-5     # [1/K]
_T_REF = 300.0      # [K]


def _thermal_models(mod):
    """One thermal instance of every model, keyed by class name."""
    return {
        "ConstantDensityEOS": mod.ConstantDensityEOS(_RHO0, thermal_expansion=_ALPHA),
        "BirchMurnaghanEOS": mod.BirchMurnaghanEOS(_RHO0, _K0, _K0P, thermal_expansion=_ALPHA),
        "VinetEOS": mod.VinetEOS(_RHO0, _K0, _K0P, thermal_expansion=_ALPHA),
        "InterpolatedEOS": mod.InterpolatedEOS(
            [0.0, 1.0e6, 2.0e6], [5000.0, 4000.0, 3000.0], thermal_expansion=_ALPHA),
    }


_MODEL_NAMES = ["ConstantDensityEOS", "BirchMurnaghanEOS", "VinetEOS", "InterpolatedEOS"]


def test_default_models_are_athermal():
    mod = _import_eos()
    eos = mod.BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    assert eos.thermal_expansion == 0.0
    assert eos.reference_temperature == pytest.approx(_T_REF)
    assert eos.calc_density(1.0e10, 2000.0) == eos.calc_density(1.0e10)


@pytest.mark.parametrize("cls_name", _MODEL_NAMES)
def test_density_at_the_reference_temperature_is_athermal(cls_name):
    eos = _thermal_models(_import_eos())[cls_name]
    assert eos.thermal_expansion == pytest.approx(_ALPHA)
    athermal = eos.calc_density(2.0e10, None, 0.5e6)
    assert eos.calc_density(2.0e10, _T_REF, 0.5e6) == pytest.approx(athermal, rel=1e-13)
    # A hotter material is less dense at the same pressure.
    assert eos.calc_density(2.0e10, 2000.0, 0.5e6) < athermal


@pytest.mark.parametrize("cls_name", ["ConstantDensityEOS", "InterpolatedEOS"])
@pytest.mark.parametrize("temperature", [100.0, 1500.0, 4000.0])
def test_thermal_expansion_factor(cls_name, temperature):
    """Models with no pressure law scale their density by exp(-alpha (T - T_ref))."""
    eos = _thermal_models(_import_eos())[cls_name]
    athermal = eos.calc_density(0.0, None, 0.5e6)
    expected = athermal * math.exp(-_ALPHA * (temperature - _T_REF))
    assert eos.calc_density(0.0, temperature, 0.5e6) == pytest.approx(expected, rel=1e-13)


@pytest.mark.parametrize("cls_name,law", [("BirchMurnaghanEOS", "birch_murnaghan_pressure"),
                                          ("VinetEOS", "vinet_pressure")])
@pytest.mark.parametrize("pressure", [0.0, 1.0e10, 1.0e11])
@pytest.mark.parametrize("temperature", [100.0, 1500.0, 4000.0])
def test_thermal_pressure_roundtrip(cls_name, law, pressure, temperature):
    """The cold law at the solved compression returns the pressure less alpha0 K0 (T - T_ref)."""
    mod = _import_eos()
    eos = _thermal_models(mod)[cls_name]
    eta = eos.calc_density(pressure, temperature) / _RHO0
    cold_pressure = pressure - _ALPHA * _K0 * (temperature - _T_REF)
    assert getattr(mod, law)(eta, _K0, _K0P) == pytest.approx(cold_pressure, rel=1e-6, abs=1.0)


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
def test_free_surface_expansion_is_alpha_delta_t(cls_name):
    """At zero pressure a small temperature rise lowers the density by alpha dT to first order."""
    eos = _thermal_models(_import_eos())[cls_name]
    delta_temp = 10.0
    relative_change = 1.0 - eos.calc_density(0.0, _T_REF + delta_temp) / _RHO0
    assert relative_change == pytest.approx(_ALPHA * delta_temp, rel=1e-2)


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
@pytest.mark.parametrize("pressure", [0.0, 1.0e10, 1.0e11])
@pytest.mark.parametrize("temperature", [None, 2500.0])
def test_bulk_modulus_matches_the_density_derivative(cls_name, pressure, temperature):
    """K = rho dP/drho, checked against a centered difference of calc_density."""
    eos = _thermal_models(_import_eos())[cls_name]
    step = 1.0e6
    density = eos.calc_density(pressure, temperature)
    slope = (eos.calc_density(pressure + step, temperature) - eos.calc_density(pressure - step, temperature))
    expected = density * 2.0 * step / slope
    assert eos.calc_bulk_modulus(pressure, temperature) == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
def test_bulk_modulus_at_the_reference_state(cls_name):
    eos = _thermal_models(_import_eos())[cls_name]
    assert eos.calc_bulk_modulus(0.0) == pytest.approx(_K0, rel=1e-12)
    # Compression stiffens the material and heating softens it.
    assert eos.calc_bulk_modulus(5.0e10) > _K0
    assert eos.calc_bulk_modulus(0.0, 2000.0) < _K0


def test_bulk_modulus_of_models_without_a_pressure_law():
    mod = _import_eos()
    assert math.isnan(mod.ConstantDensityEOS(_RHO0).calc_bulk_modulus(1.0e10, 1000.0))
    table = mod.InterpolatedEOS([0.0, 2.0e6], [5000.0, 3000.0], bulk_modulus=[3.0e11, 1.0e11])
    assert table.calc_bulk_modulus(0.0, None, 1.0e6) == pytest.approx(2.0e11)


@pytest.mark.parametrize("cls_name", _MODEL_NAMES)
def test_thermal_terms_survive_config_and_binary_roundtrips(cls_name):
    mod = _import_eos()
    eos = _thermal_models(mod)[cls_name]
    density_before = eos.calc_density(3.0e10, 1800.0, 0.5e6)

    cfg = eos.get_config_dict()
    assert cfg["thermal_expansion_1_k"] == pytest.approx(_ALPHA)
    assert cfg["reference_temperature_k"] == pytest.approx(_T_REF)
    rebuilt = mod.make_material_eos(cfg["model"], {key: value for key, value in cfg.items() if key != "model"})
    assert rebuilt.get_config_dict() == cfg
    assert rebuilt.calc_density(3.0e10, 1800.0, 0.5e6) == pytest.approx(density_before, rel=1e-13)

    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        eos.save_binary(path)
        loaded = mod.make_material_eos(eos.model_name)
        loaded.load_binary(path)
        assert loaded.thermal_expansion == pytest.approx(_ALPHA)
        assert loaded.calc_density(3.0e10, 1800.0, 0.5e6) == pytest.approx(density_before, rel=1e-13)
    finally:
        os.unlink(path)


def test_reference_temperature_is_configurable():
    mod = _import_eos()
    eos = mod.make_material_eos("constant", {
        "reference_density_kg_m3": _RHO0,
        "thermal_expansion_1_k": _ALPHA,
        "reference_temperature_k": 1600.0})
    assert eos.reference_temperature == pytest.approx(1600.0)
    assert eos.calc_density(0.0, 1600.0) == pytest.approx(_RHO0, rel=1e-13)
