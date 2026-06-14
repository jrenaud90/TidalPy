"""
Tests for TidalPy.structures_x.worlds — BaseWorld and StarWorld (Phase 8b).

Covers construction, bulk geometry, equilibrium temperature, Stefan-Boltzmann
luminosity (stars), config dict, binary round-trip, and isinstance checks.

Requires the Cython extensions to be compiled first::

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
    from TidalPy.structures_x.worlds import base as _mod
    return _mod


def _import_stellar():
    from TidalPy.structures_x.worlds import stellar as _mod
    return _mod


# Reference values (Earth, MKS)
_R_EARTH = 6.371e6      # [m]
_M_EARTH = 5.972e24     # [kg]
_G       = 6.674e-11    # [m^3 kg^-1 s^-2]
_SIGMA   = 5.670374419e-8  # [W/m^2/K^4]
_R_SUN   = 6.957e8      # [m]
_M_SUN   = 1.989e30     # [kg]
_T_SUN   = 5772.0       # [K]


def _make_earth(**kw):
    mod = _import_base()
    defaults = dict(name="Earth", radius_m=_R_EARTH, mass_kg=_M_EARTH,
                    world_type="terrestrial", albedo=0.3, emissivity=1.0)
    defaults.update(kw)
    return mod.BaseWorld(**defaults)


# =====================================================================================================================
# Construction and properties
# =====================================================================================================================
def test_base_world_construction():
    w = _make_earth()
    assert w.name           == "Earth"
    assert w.world_type     == "terrestrial"
    assert w.radius         == pytest.approx(_R_EARTH)
    assert w.mass           == pytest.approx(_M_EARTH)
    assert w.albedo         == pytest.approx(0.3)
    assert w.emissivity     == pytest.approx(1.0)
    assert w.obliquity      == pytest.approx(0.0)
    assert w.spin_frequency == pytest.approx(0.0)


def test_base_world_surface_gravity():
    w = _make_earth()
    expected = _G * _M_EARTH / (_R_EARTH ** 2)
    assert w.calc_surface_gravity() == pytest.approx(expected, rel=1e-3)


def test_base_world_escape_velocity():
    w = _make_earth()
    expected = math.sqrt(2.0 * _G * _M_EARTH / _R_EARTH)
    assert w.calc_escape_velocity() == pytest.approx(expected, rel=1e-3)


def test_base_world_mean_density():
    w = _make_earth()
    vol = (4.0 / 3.0) * math.pi * _R_EARTH ** 3
    assert w.calc_mean_density() == pytest.approx(_M_EARTH / vol, rel=1e-9)


def test_base_world_equilibrium_temperature():
    # Earth at ~1361 W/m^2 with A=0.3, eps=1 => ~255 K
    w = _make_earth()
    flux = 1361.0
    expected = ((1.0 - 0.3) * flux / (4.0 * 1.0 * _SIGMA)) ** 0.25
    assert w.calc_equilibrium_temperature(flux) == pytest.approx(expected, rel=1e-3)


def test_base_world_equilibrium_temperature_zero_flux():
    w = _make_earth()
    assert w.calc_equilibrium_temperature(0.0) == pytest.approx(0.0)
    assert w.calc_equilibrium_temperature(-5.0) == pytest.approx(0.0)


def test_base_world_setters():
    w = _make_earth()
    w.set_spin_frequency(7.29e-5)
    w.set_obliquity(0.41)
    assert w.spin_frequency == pytest.approx(7.29e-5)
    assert w.obliquity      == pytest.approx(0.41)


# =====================================================================================================================
# Config + binary
# =====================================================================================================================
def test_base_world_config_dict():
    w   = _make_earth(albedo=0.25, emissivity=0.9)
    cfg = w.get_config_dict()
    for key in ("name", "world_type", "radius_m", "mass_kg", "albedo",
                "emissivity", "obliquity_rad", "spin_frequency_rad_s"):
        assert key in cfg
    assert cfg["name"]       == "Earth"
    assert cfg["albedo"]     == pytest.approx(0.25)
    assert cfg["emissivity"] == pytest.approx(0.9)


def test_base_world_binary_roundtrip():
    mod = _import_base()
    w1  = _make_earth(albedo=0.31, emissivity=0.95, obliquity_rad=0.41,
                      spin_frequency_rad_s=7.29e-5)
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        w1.save_binary(path)
        w2 = mod.BaseWorld("placeholder", 1.0, 1.0)
        w2.load_binary(path)
        assert w2.name           == "Earth"
        assert w2.world_type     == "terrestrial"
        assert w2.radius         == pytest.approx(_R_EARTH)
        assert w2.mass           == pytest.approx(_M_EARTH)
        assert w2.albedo         == pytest.approx(0.31)
        assert w2.emissivity     == pytest.approx(0.95)
        assert w2.obliquity      == pytest.approx(0.41)
        assert w2.spin_frequency == pytest.approx(7.29e-5)
    finally:
        os.unlink(path)


def test_base_world_load_file_not_found():
    w = _make_earth()
    with pytest.raises(FileNotFoundError):
        w.load_binary("/nonexistent/path/world.tpyb")


# =====================================================================================================================
# StarWorld
# =====================================================================================================================
def test_star_construction_and_luminosity():
    mod = _import_stellar()
    star = mod.StarWorld("Sun", _R_SUN, _M_SUN, effective_temperature_k=_T_SUN)
    assert star.world_type == "star"
    assert star.effective_temperature == pytest.approx(_T_SUN)
    # Luminosity derived from T via Stefan-Boltzmann should be ~3.8e26 W.
    expected_L = 4.0 * math.pi * _R_SUN ** 2 * _SIGMA * _T_SUN ** 4
    assert star.luminosity == pytest.approx(expected_L, rel=1e-3)
    assert star.luminosity == pytest.approx(3.83e26, rel=0.05)


def test_star_temperature_luminosity_roundtrip():
    mod = _import_stellar()
    star = mod.StarWorld("Sun", _R_SUN, _M_SUN, effective_temperature_k=_T_SUN)
    L = star.luminosity
    assert star.calc_temperature_from_luminosity(L) == pytest.approx(_T_SUN, rel=1e-6)


def test_star_setters_keep_consistent():
    mod = _import_stellar()
    star = mod.StarWorld("Sun", _R_SUN, _M_SUN, effective_temperature_k=_T_SUN)
    star.set_effective_temperature(6000.0)
    assert star.effective_temperature == pytest.approx(6000.0)
    expected_L = 4.0 * math.pi * _R_SUN ** 2 * _SIGMA * 6000.0 ** 4
    assert star.luminosity == pytest.approx(expected_L, rel=1e-6)


def test_star_binary_roundtrip():
    mod = _import_stellar()
    s1  = mod.StarWorld("Sun", _R_SUN, _M_SUN, effective_temperature_k=_T_SUN)
    L_before = s1.luminosity
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        s1.save_binary(path)
        s2 = mod.StarWorld("placeholder", 1.0, 1.0)
        s2.load_binary(path)
        assert s2.name == "Sun"
        assert s2.effective_temperature == pytest.approx(_T_SUN)
        assert s2.luminosity == pytest.approx(L_before, rel=1e-12)
    finally:
        os.unlink(path)


# =====================================================================================================================
# isinstance
# =====================================================================================================================
def test_world_isinstance_chain():
    from TidalPy.Utilities_x.classes_x.classes import StructureBase, TidalPyBaseClass
    w = _make_earth()
    assert isinstance(w, StructureBase)
    assert isinstance(w, TidalPyBaseClass)


def test_star_is_base_world():
    from TidalPy.structures_x.worlds.base import BaseWorld
    mod = _import_stellar()
    star = mod.StarWorld("Sun", _R_SUN, _M_SUN)
    assert isinstance(star, BaseWorld)
