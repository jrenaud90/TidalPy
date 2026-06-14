"""
Tests for TidalPy.Material_x.eos.material_eos — material EOS model hierarchy.

Covers the constant-density, Birch-Murnaghan, Vinet, and interpolated models:
density evaluation, analytic-inversion cross-checks (density-from-pressure round
trips through the forward pressure law), the name factory, config dicts, binary
round-trips, and isinstance checks.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

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
    eos = _import_eos().ConstantDensityEOS(reference_density_kg_m3=_RHO0)
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
