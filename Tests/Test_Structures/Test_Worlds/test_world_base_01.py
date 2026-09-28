"""BaseWorld and StarWorld: construction, bulk geometry, equilibrium temperature, luminosity, config, and binary."""
import math

import pytest

from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Utilities.classes.classes import StructureBase, TidalPyBaseClass


_R_EARTH = 6.371e6         # [m]
_M_EARTH = 5.972e24        # [kg]
_G       = 6.674e-11       # [m^3 kg^-1 s^-2]
_SIGMA   = 5.670374419e-8  # [W/m^2/K^4]
_R_SUN   = 6.957e8         # [m]
_M_SUN   = 1.989e30        # [kg]
_T_SUN   = 5772.0          # [K]


def _make_earth(**kw):
    defaults = dict(name="Earth", radius=_R_EARTH, mass=_M_EARTH,
                    world_type="terrestrial", albedo=0.3, emissivity=1.0)
    defaults.update(kw)
    return BaseWorld(**defaults)


def _make_sun():
    return StarWorld("Sun", _R_SUN, _M_SUN, effective_temperature=_T_SUN)


# =====================================================================================================================
# BaseWorld
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


@pytest.mark.parametrize("method, expected, rel", [
    ("calc_surface_gravity", _G * _M_EARTH / _R_EARTH ** 2, 1e-3),
    ("calc_escape_velocity", math.sqrt(2.0 * _G * _M_EARTH / _R_EARTH), 1e-3),
    ("calc_mean_density", _M_EARTH / ((4.0 / 3.0) * math.pi * _R_EARTH ** 3), 1e-9),
])
def test_base_world_bulk_quantities(method, expected, rel):
    assert getattr(_make_earth(), method)() == pytest.approx(expected, rel=rel)


@pytest.mark.parametrize("flux, expected", [
    (1361.0, ((1.0 - 0.3) * 1361.0 / (4.0 * 1.0 * _SIGMA)) ** 0.25),
    (0.0, 0.0),
    (-5.0, 0.0),
], ids=["earth", "zero_flux", "negative_flux"])
def test_base_world_equilibrium_temperature(flux, expected):
    assert _make_earth().calc_equilibrium_temperature(flux) == pytest.approx(expected, rel=1e-3)


def test_base_world_setters():
    w = _make_earth()
    w.set_spin_frequency(7.29e-5)
    w.set_obliquity(0.41)
    assert w.spin_frequency == pytest.approx(7.29e-5)
    assert w.obliquity      == pytest.approx(0.41)


def test_base_world_config_dict():
    w   = _make_earth(albedo=0.25, emissivity=0.9)
    cfg = w.get_config_dict()
    for key in ("schema_version", "name", "type", "radius_m", "mass_kg", "albedo",
                "emissivity", "obliquity_rad", "spin_frequency_rad_s"):
        assert key in cfg
    assert cfg["name"]       == "Earth"
    assert cfg["albedo"]     == pytest.approx(0.25)
    assert cfg["emissivity"] == pytest.approx(0.9)


def test_base_world_binary_roundtrip(tmp_path):
    path = str(tmp_path / "world.tpyb")
    _make_earth(albedo=0.31, emissivity=0.95, obliquity=0.41, spin_frequency=7.29e-5).save_binary(path)
    w2 = BaseWorld("placeholder", 1.0, 1.0)
    w2.load_binary(path)
    assert w2.name           == "Earth"
    assert w2.world_type     == "terrestrial"
    assert w2.radius         == pytest.approx(_R_EARTH)
    assert w2.mass           == pytest.approx(_M_EARTH)
    assert w2.albedo         == pytest.approx(0.31)
    assert w2.emissivity     == pytest.approx(0.95)
    assert w2.obliquity      == pytest.approx(0.41)
    assert w2.spin_frequency == pytest.approx(7.29e-5)


def test_base_world_load_file_not_found():
    with pytest.raises(FileNotFoundError):
        _make_earth().load_binary("/nonexistent/path/world.tpyb")


# =====================================================================================================================
# StarWorld
# =====================================================================================================================
def test_star_construction_and_luminosity():
    star = _make_sun()
    assert star.world_type == "star"
    assert star.effective_temperature == pytest.approx(_T_SUN)
    expected_L = 4.0 * math.pi * _R_SUN ** 2 * _SIGMA * _T_SUN ** 4
    assert star.luminosity == pytest.approx(expected_L, rel=1e-3)
    assert star.luminosity == pytest.approx(3.83e26, rel=0.05)


def test_star_temperature_luminosity_roundtrip():
    star = _make_sun()
    assert star.calc_temperature_from_luminosity(star.luminosity) == pytest.approx(_T_SUN, rel=1e-6)


def test_star_setters_keep_consistent():
    star = _make_sun()
    star.set_effective_temperature(6000.0)
    assert star.effective_temperature == pytest.approx(6000.0)
    expected_L = 4.0 * math.pi * _R_SUN ** 2 * _SIGMA * 6000.0 ** 4
    assert star.luminosity == pytest.approx(expected_L, rel=1e-6)


def test_star_binary_roundtrip(tmp_path):
    path = str(tmp_path / "star.tpyb")
    s1 = _make_sun()
    L_before = s1.luminosity
    s1.save_binary(path)
    s2 = StarWorld("placeholder", 1.0, 1.0)
    s2.load_binary(path)
    assert s2.name == "Sun"
    assert s2.effective_temperature == pytest.approx(_T_SUN)
    assert s2.luminosity == pytest.approx(L_before, rel=1e-12)


def test_world_isinstance_chain():
    w = _make_earth()
    assert isinstance(w, StructureBase)
    assert isinstance(w, TidalPyBaseClass)
    assert isinstance(StarWorld("Sun", _R_SUN, _M_SUN), BaseWorld)
