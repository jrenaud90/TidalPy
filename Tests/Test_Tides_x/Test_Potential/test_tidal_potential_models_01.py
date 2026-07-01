"""Tests for the tidal-potential model hierarchy (TidalPy.Tides_x.potential.tidal_potential).

Covers the name factory (aliases, unknown-name ValueError, rich-subclass return), num_modes, the
config dict, a binary round-trip, the isinstance chain, and that calc_modes reproduces the legacy
synchronous-low-e and NSR-modes potentials.
"""
import os
import tempfile

import numpy as np
import pytest


def _import():
    import TidalPy.Tides_x.potential.tidal_potential as mod
    return mod


# =====================================================================================================================
# Factory: names, aliases, rich-subclass return
# =====================================================================================================================
@pytest.mark.parametrize("name,cls_attr", [
    ("sync_low_e", "SyncLowEPotential"),
    ("synchronous_low_e", "SyncLowEPotential"),
    ("simple", "SyncLowEPotential"),
    ("nsr_modes", "NSRModesPotential"),
    ("nsr", "NSRModesPotential"),
    ("NSR_MODES", "NSRModesPotential"),   # case-insensitive
])
def test_make_tidal_potential_returns_subclass(name, cls_attr):
    mod = _import()
    model = mod.make_tidal_potential(name)
    assert isinstance(model, getattr(mod, cls_attr))
    assert isinstance(model, mod.TidalPotentialBase)


def test_make_tidal_potential_unknown_name_raises():
    with pytest.raises(ValueError):
        _import().make_tidal_potential("not_a_model")


def test_num_modes():
    mod = _import()
    assert mod.make_tidal_potential("sync_low_e").num_modes == 1
    assert mod.make_tidal_potential("nsr_modes").num_modes == 9


# =====================================================================================================================
# Config dict + use_static
# =====================================================================================================================
def test_config_dict():
    mod = _import()
    sync = mod.make_tidal_potential("sync_low_e")
    assert sync.get_config_dict()["model"] == "sync_low_e"
    nsr = mod.make_tidal_potential("nsr_modes", {"use_static": True})
    cfg = nsr.get_config_dict()
    assert cfg["model"] == "nsr_modes"
    assert cfg["use_static"] is True
    assert nsr.use_static is True


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
@pytest.mark.parametrize("name,config", [
    ("sync_low_e", None),
    ("nsr_modes", {"use_static": True}),
    ("nsr_modes", {"use_static": False}),
])
def test_binary_round_trip(name, config):
    mod = _import()
    model = mod.make_tidal_potential(name, config)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, f"{name}.tpyb")
        model.save_binary(path)
        reloaded = mod.make_tidal_potential(name)
        reloaded.load_binary(path)
    assert reloaded.get_config_dict() == model.get_config_dict()


# =====================================================================================================================
# isinstance chain
# =====================================================================================================================
def test_isinstance_chain():
    from TidalPy.Utilities_x.classes_x.classes import PhysicsBase, TidalPyBaseClass
    mod = _import()
    model = mod.make_tidal_potential("nsr_modes")
    assert isinstance(model, mod.NSRModesPotential)
    assert isinstance(model, mod.TidalPotentialBase)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)


# =====================================================================================================================
# calc_modes physics vs the legacy providers
# =====================================================================================================================
def test_sync_low_e_matches_legacy():
    from TidalPy.tides.potential.synchronous_low_e import tidal_potential as legacy
    mod = _import()
    sync = mod.make_tidal_potential("sync_low_e")
    n = 2.0 * np.pi / 86400.0
    radius, lon, colat, t, ecc, host, a = 6.0e6, 0.7, 1.3, 1234.0, 0.05, 1.0e26, 1.0e9
    _, pots = sync.calc_modes(n, 0.0, ecc, host, a, radius, colat, lon, t)
    _, _, ptup = legacy(radius, lon, colat, t, n, ecc, host, a)
    assert np.allclose(np.asarray(pots[0]), np.asarray(ptup['n']), rtol=1e-12)


def test_nsr_modes_matches_legacy():
    from TidalPy.tides.potential.nsr_modes_med_eccen_no_obliquity import tidal_potential as legacy
    mod = _import()
    nsr = mod.make_tidal_potential("nsr_modes")
    n = 2.0 * np.pi / 86400.0
    o = 1.5 * n
    radius, lon, colat, t, ecc, host, a = 6.0e6, 0.7, 1.3, 1234.0, 0.05, 1.0e26, 1.0e9
    modes, pots = nsr.calc_modes(n, o, ecc, host, a, radius, colat, lon, t)
    freqs_by_name, modes_by_name, ptup_by_mode = legacy(radius, lon, colat, t, n, o, ecc, host, a)
    mode_names = ('n', '2n', '3n', '2o+n', '2o-n', '2o-2n', '2o-3n', '2o-4n', '2o-5n')
    for i, mode_name in enumerate(mode_names):
        assert np.isclose(modes[i], modes_by_name[mode_name], rtol=1e-12)
        assert np.allclose(np.asarray(pots[i]), np.asarray(ptup_by_mode[mode_name]), rtol=1e-10, atol=1e-30)
