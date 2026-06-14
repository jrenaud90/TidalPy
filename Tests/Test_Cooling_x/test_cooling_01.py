"""
Tests for TidalPy.cooling_x.cooling.

Covers construction, model names and aliases, cooling physics (cross-checked
against the legacy ``TidalPy/cooling/cooling_models.py``), the ``make_cooling``
factory, the ``CoolingResult`` container, vectorization, config dicts, binary
round-trip, TOML config save, and the isinstance chain.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import os
import tempfile

import numpy as np
import pytest


# =====================================================================================================================
# Imports
# =====================================================================================================================
def _import_cooling():
    try:
        from TidalPy import cooling_x as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.cooling_x not compiled — run uv pip install first."
        )


def _legacy():
    from TidalPy.cooling import cooling_models as legacy
    return legacy


# Representative MKS inputs for a convecting silicate mantle (high Rayleigh number).
_DT    = 1000.0    # temperature drop [K]
_THICK = 1.0e6     # layer thickness [m]
_G     = 9.8       # gravity [m/s^2]
_RHO   = 3300.0    # density [kg/m^3]
_VISC  = 1.0e21    # viscosity [Pa·s]
_K     = 4.0       # thermal conductivity [W/m/K]
_DIFF  = 1.0e-6    # thermal diffusivity [m^2/s]
_EXP   = 3.0e-5    # thermal expansivity [1/K]

# Default convection scaling parameters (match the C++ / legacy defaults).
_ALPHA = 1.0
_BETA  = 1.0 / 3.0
_RACR  = 1100.0


def _inputs():
    """Return the 8 calc_cooling inputs in order."""
    return (_DT, _THICK, _G, _RHO, _VISC, _K, _DIFF, _EXP)


# =====================================================================================================================
# Construction & names
# =====================================================================================================================
def test_model_names():
    """Each model reports its canonical lower-case model name."""
    mod = _import_cooling()
    assert mod.OffCooling().model_name == "off"
    assert mod.ConvectiveCooling().model_name == "convection"
    assert mod.ConductiveCooling().model_name == "conduction"


def test_abstract_base_not_instantiable():
    """CoolingBase is abstract and cannot be constructed directly."""
    mod = _import_cooling()
    with pytest.raises(TypeError):
        mod.CoolingBase()


# =====================================================================================================================
# Cooling physics (cross-check vs legacy cooling_models)
# =====================================================================================================================
def test_off_cooling():
    """Off model: zero flux, boundary layer half the thickness, Ra=0, Nu=1."""
    mod = _import_cooling()
    res = mod.OffCooling().calc_cooling(*_inputs())
    assert res.cooling_flux == 0.0
    assert res.boundary_layer_thickness == pytest.approx(0.5 * _THICK)
    assert res.rayleigh == 0.0
    assert res.nusselt == 1.0


def test_conduction_matches_legacy():
    """Conduction matches the legacy conduction() model."""
    mod = _import_cooling()
    legacy = _legacy()
    res = mod.ConductiveCooling().calc_cooling(*_inputs())
    flux, blt, ray, nu = legacy.conduction(_DT, _K, _THICK)
    assert res.cooling_flux == pytest.approx(flux)
    assert res.boundary_layer_thickness == pytest.approx(blt)
    assert res.rayleigh == pytest.approx(ray)
    assert res.nusselt == pytest.approx(nu)
    # Closed form: q = k * dT / thickness.
    assert res.cooling_flux == pytest.approx(_K * _DT / _THICK)


def test_convection_matches_legacy_high_rayleigh():
    """Convection matches the legacy model in the active (Nu > 2) regime."""
    mod = _import_cooling()
    legacy = _legacy()
    res = mod.ConvectiveCooling().calc_cooling(*_inputs())
    flux, blt, ray, nu = legacy.convection(
        _DT, _VISC, _K, _DIFF, _EXP, _THICK, _G, _RHO, _ALPHA, _BETA, _RACR)
    assert res.rayleigh == pytest.approx(ray, rel=1e-12)
    assert res.nusselt == pytest.approx(nu, rel=1e-12)
    assert res.boundary_layer_thickness == pytest.approx(blt, rel=1e-12)
    assert res.cooling_flux == pytest.approx(flux, rel=1e-12)
    # This regime really is convecting.
    assert res.nusselt > 2.0
    assert res.rayleigh > _RACR


def test_convection_matches_legacy_low_rayleigh():
    """Convection matches the legacy model in the floored (Nu == 2) regime."""
    mod = _import_cooling()
    legacy = _legacy()
    # High viscosity, thin layer -> sub-critical Rayleigh -> Nu floored at 2.
    dt, visc, thick = 100.0, 1.0e24, 1.0e4
    res = mod.ConvectiveCooling().calc_cooling(dt, thick, _G, _RHO, visc, _K, _DIFF, _EXP)
    flux, blt, ray, nu = legacy.convection(
        dt, visc, _K, _DIFF, _EXP, thick, _G, _RHO, _ALPHA, _BETA, _RACR)
    assert res.nusselt == pytest.approx(nu)
    assert res.nusselt == pytest.approx(2.0)
    assert res.cooling_flux == pytest.approx(flux, rel=1e-12)


def test_convection_parameters_affect_result():
    """A larger critical Rayleigh number reduces the Nusselt number."""
    mod = _import_cooling()
    base = mod.ConvectiveCooling().calc_cooling(*_inputs())
    stiff = mod.ConvectiveCooling(critical_rayleigh=1.0e5).calc_cooling(*_inputs())
    assert stiff.nusselt < base.nusselt


# =====================================================================================================================
# Parameter getters
# =====================================================================================================================
def test_convective_parameters():
    mod = _import_cooling()
    cv = mod.ConvectiveCooling(0.5, 0.25, 2000.0)
    assert cv.convection_alpha == pytest.approx(0.5)
    assert cv.convection_beta == pytest.approx(0.25)
    assert cv.critical_rayleigh == pytest.approx(2000.0)


# =====================================================================================================================
# CoolingResult container
# =====================================================================================================================
def test_cooling_result_container():
    mod = _import_cooling()
    res = mod.ConductiveCooling().calc_cooling(*_inputs())
    d = res.to_dict()
    assert set(d) == {"cooling_flux", "boundary_layer_thickness", "rayleigh", "nusselt"}
    # Iteration yields the four outputs in order.
    flux, blt, ray, nu = res
    assert flux == res.cooling_flux
    assert blt == res.boundary_layer_thickness
    assert ray == res.rayleigh
    assert nu == res.nusselt
    assert "CoolingResult" in repr(res)


# =====================================================================================================================
# Factory
# =====================================================================================================================
@pytest.mark.parametrize("alias,canonical", [
    ("off",         "off"),
    ("none",        "off"),
    ("OFF",         "off"),
    ("convection",  "convection"),
    ("convective",  "convection"),
    ("Convection",  "convection"),
    ("conduction",  "conduction"),
    ("conductive",  "conduction"),
    ("CONDUCTION",  "conduction"),
])
def test_make_cooling_aliases(alias, canonical):
    """The factory resolves aliases (case-insensitive) to canonical models."""
    mod = _import_cooling()
    assert mod.make_cooling(alias).model_name == canonical


def test_make_cooling_unknown_raises():
    """An unknown model name raises ValueError (from the C++ enum mapping)."""
    mod = _import_cooling()
    with pytest.raises(ValueError):
        mod.make_cooling("not_a_real_model")


def test_make_cooling_with_config():
    """The factory forwards convection config parameters."""
    mod = _import_cooling()
    cv = mod.make_cooling("convection", {
        "convection_alpha": 0.6, "convection_beta": 0.28, "critical_rayleigh": 1600.0})
    assert cv.convection_alpha == pytest.approx(0.6)
    assert cv.convection_beta == pytest.approx(0.28)
    assert cv.critical_rayleigh == pytest.approx(1600.0)


def test_make_cooling_returns_rich_subclass():
    """The factory returns the correct concrete Python subclass (not the base)."""
    mod = _import_cooling()
    assert type(mod.make_cooling("off")).__name__ == "OffCooling"
    assert type(mod.make_cooling("convection")).__name__ == "ConvectiveCooling"
    assert type(mod.make_cooling("conduction")).__name__ == "ConductiveCooling"


def test_make_cooling_adopted_object_is_usable():
    """An object built by the C++ enum factory is fully usable and round-trips."""
    mod = _import_cooling()
    cv = mod.make_cooling("convection", {"critical_rayleigh": 1600.0})
    assert cv.critical_rayleigh == pytest.approx(1600.0)
    res = cv.calc_cooling(*_inputs())
    assert res.cooling_flux > 0.0
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "adopted.tpyb")
        cv.save_binary(path)
        restored = mod.ConvectiveCooling()
        restored.load_binary(path)
    assert restored.get_config_dict() == cv.get_config_dict()


# =====================================================================================================================
# Config dict
# =====================================================================================================================
def test_config_dict_keys():
    """Each model's config dict carries the expected keys."""
    mod = _import_cooling()
    assert set(mod.OffCooling().get_config_dict()) == {"model_name"}
    assert set(mod.ConductiveCooling().get_config_dict()) == {"model_name"}
    assert set(mod.ConvectiveCooling().get_config_dict()) == {
        "model_name", "convection_alpha", "convection_beta", "critical_rayleigh"}


def test_save_config_writes_toml():
    """save_config writes a TOML file that round-trips through toml.load."""
    import toml
    mod = _import_cooling()
    cv = mod.ConvectiveCooling(0.6, 0.28, 1600.0)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "cool.toml")
        cv.save_config(path)
        loaded = toml.load(path)
    assert loaded["model_name"] == "convection"
    assert loaded["convection_alpha"] == pytest.approx(0.6)
    assert loaded["critical_rayleigh"] == pytest.approx(1600.0)


# =====================================================================================================================
# Vectorized methods and direct convenience functions
# =====================================================================================================================
def test_convenience_scalar_matches_class():
    """The convenience function matches the class for scalar input."""
    mod = _import_cooling()
    got = mod.convective(_DT, _VISC, _THICK, _G, _RHO, _K, _DIFF, _EXP)
    expected = mod.ConvectiveCooling().calc_cooling(*_inputs())
    assert isinstance(got.cooling_flux, float)
    assert got.cooling_flux == pytest.approx(expected.cooling_flux)
    assert got.nusselt == pytest.approx(expected.nusselt)


def test_convenience_vectorize_temperature():
    """Array delta_temp at constant state returns matching ndarrays."""
    mod = _import_cooling()
    inst = mod.ConvectiveCooling()
    temps = np.array([100.0, 500.0, 1000.0, 2000.0])
    got = mod.convective(temps, _VISC, _THICK, _G, _RHO, _K, _DIFF, _EXP)
    assert isinstance(got.cooling_flux, np.ndarray)
    assert got.cooling_flux.shape == (4,)
    assert got.cooling_flux.dtype == np.float64
    for i in range(4):
        single = inst.calc_cooling(temps[i], _THICK, _G, _RHO, _VISC, _K, _DIFF, _EXP)
        assert got.cooling_flux[i] == pytest.approx(single.cooling_flux)
        assert got.nusselt[i] == pytest.approx(single.nusselt)


def test_convenience_vectorize_viscosity():
    """Array viscosity at constant state returns matching ndarrays."""
    mod = _import_cooling()
    inst = mod.ConvectiveCooling()
    viscs = np.array([1.0e19, 1.0e20, 1.0e21])
    got = mod.convective(_DT, viscs, _THICK, _G, _RHO, _K, _DIFF, _EXP)
    assert got.cooling_flux.shape == (3,)
    for i in range(3):
        single = inst.calc_cooling(_DT, _THICK, _G, _RHO, viscs[i], _K, _DIFF, _EXP)
        assert got.cooling_flux[i] == pytest.approx(single.cooling_flux)


def test_convenience_vectorize_all_and_broadcast():
    """All-array and mixed (broadcast) inputs return the broadcast shape."""
    mod = _import_cooling()
    temps = np.array([500.0, 1000.0, 2000.0])
    viscs = np.array([1.0e20, 1.0e21, 1.0e22])
    got_all = mod.convective(temps, viscs, _THICK, _G, _RHO, _K, _DIFF, _EXP)
    assert got_all.cooling_flux.shape == (3,)
    inst = mod.ConvectiveCooling()
    for i in range(3):
        single = inst.calc_cooling(temps[i], _THICK, _G, _RHO, viscs[i], _K, _DIFF, _EXP)
        assert got_all.cooling_flux[i] == pytest.approx(single.cooling_flux)


def test_convenience_preserves_2d_shape():
    """A 2-D input shape is preserved in the output arrays."""
    mod = _import_cooling()
    temps = np.array([[100.0, 500.0], [1000.0, 2000.0]])
    got = mod.convective(temps, _VISC, _THICK, _G, _RHO, _K, _DIFF, _EXP)
    assert got.cooling_flux.shape == (2, 2)
    assert got.nusselt.shape == (2, 2)


def test_conductive_off_convenience():
    """The conduction and off convenience functions match their classes."""
    mod = _import_cooling()
    cd = mod.conductive(_DT, _THICK, _K)
    assert cd.cooling_flux == pytest.approx(_K * _DT / _THICK)
    off = mod.cooling_off(_DT, _THICK)
    assert off.cooling_flux == 0.0
    assert off.boundary_layer_thickness == pytest.approx(0.5 * _THICK)
    # Conduction convenience accepts arrays for delta_temp.
    arr = mod.conductive(np.array([100.0, 200.0]), _THICK, _K)
    assert arr.cooling_flux.shape == (2,)


def test_class_vectorize_methods_match_scalar():
    """The class-level vectorized methods match element-wise scalar calls."""
    mod = _import_cooling()
    cv = mod.ConvectiveCooling()
    temps = np.array([100.0, 1000.0, 2000.0])
    viscs = np.array([1.0e20, 1.0e21, 1.0e22])

    by_temp = cv.calc_cooling_vectorize_temperature(temps, _THICK, _G, _RHO, _VISC, _K, _DIFF, _EXP)
    by_visc = cv.calc_cooling_vectorize_viscosity(_DT, _THICK, _G, _RHO, viscs, _K, _DIFF, _EXP)
    by_all  = cv.calc_cooling_vectorize_all(temps, _THICK, _G, _RHO, viscs, _K, _DIFF, _EXP)
    for i in range(3):
        s_t = cv.calc_cooling(temps[i], _THICK, _G, _RHO, _VISC, _K, _DIFF, _EXP)
        s_v = cv.calc_cooling(_DT, _THICK, _G, _RHO, viscs[i], _K, _DIFF, _EXP)
        s_a = cv.calc_cooling(temps[i], _THICK, _G, _RHO, viscs[i], _K, _DIFF, _EXP)
        assert by_temp.cooling_flux[i] == pytest.approx(s_t.cooling_flux)
        assert by_visc.cooling_flux[i] == pytest.approx(s_v.cooling_flux)
        assert by_all.cooling_flux[i] == pytest.approx(s_a.cooling_flux)


def test_class_vectorize_size_mismatch_raises():
    """Mismatched input lengths raise ValueError from the C++ layer."""
    mod = _import_cooling()
    cv = mod.ConvectiveCooling()
    with pytest.raises(ValueError):
        cv.calc_cooling_vectorize_all(
            np.ones(3), _THICK, _G, _RHO, np.ones(2), _K, _DIFF, _EXP)


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
@pytest.mark.parametrize("name,cls,args", [
    ("off", "OffCooling", ()),
    ("conduction", "ConductiveCooling", ()),
    ("convection", "ConvectiveCooling", (0.6, 0.28, 1600.0)),
])
def test_binary_round_trip(name, cls, args):
    """save_binary / load_binary preserves model name and parameters."""
    mod = _import_cooling()
    cls_obj = getattr(mod, cls)
    original = cls_obj(*args)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, f"{name}.tpyb")
        original.save_binary(path)
        restored = cls_obj()
        restored.load_binary(path)
    assert restored.model_name == name
    assert restored.get_config_dict() == original.get_config_dict()


def test_load_binary_missing_file_raises():
    """Loading a non-existent binary raises FileNotFoundError."""
    mod = _import_cooling()
    with pytest.raises(FileNotFoundError):
        mod.ConvectiveCooling().load_binary("does_not_exist_12345.tpyb")


# =====================================================================================================================
# isinstance chain
# =====================================================================================================================
@pytest.mark.parametrize("cls", ["OffCooling", "ConvectiveCooling", "ConductiveCooling"])
def test_isinstance_chain(cls):
    """Every model is a CoolingBase, PhysicsBase and TidalPyBaseClass."""
    mod = _import_cooling()
    from TidalPy.Utilities_x.classes_x.classes import PhysicsBase, TidalPyBaseClass
    inst = getattr(mod, cls)()
    assert isinstance(inst, mod.CoolingBase)
    assert isinstance(inst, PhysicsBase)
    assert isinstance(inst, TidalPyBaseClass)
