"""Tests for the global (1D) tidal dissipation model hierarchy (TidalPy.Tides_x.classes.tide).

Covers the name factory (incl. aliases, unknown-name ValueError, and rich-subclass return),
the complex Love number / -Im[k] laws for each model, per-degree parameters, config dicts,
binary round-trips, and the isinstance chain.
"""
import math
import os
import tempfile
from math import isclose

import pytest

from TidalPy.Tides_x.love.love import LoveNumbers


def _import():
    import TidalPy.Tides_x.classes as mod
    return mod


# =====================================================================================================================
# Factory: names, aliases, rich-subclass return
# =====================================================================================================================
@pytest.mark.parametrize("name,cls_attr", [
    ("rheology", "RheologyTide"),
    ("cpl", "FixedQTide"),
    ("fixed_q", "FixedQTide"),
    ("constant_phase_lag", "FixedQTide"),
    ("ctl", "FixedLagTide"),
    ("fixed_dt", "FixedLagTide"),
    ("constant_time_lag", "FixedLagTide"),
    ("ctl_q", "CTLQTide"),
    ("fixed_dt_q", "CTLQTide"),
    ("CTL_Q", "CTLQTide"),  # case-insensitive
])
def test_make_tide_returns_subclass(name, cls_attr):
    mod = _import()
    model = mod.make_tide(name)
    assert isinstance(model, getattr(mod, cls_attr))
    assert isinstance(model, mod.TideBase)


def test_make_tide_unknown_name_raises():
    with pytest.raises(ValueError):
        _import().make_tide("not_a_model")


# =====================================================================================================================
# needs_radial_solve flag
# =====================================================================================================================
@pytest.mark.parametrize("name,expected", [
    ("rheology", True),
    ("cpl", False),
    ("ctl", False),
    ("ctl_q", False),
])
def test_needs_radial_solve(name, expected):
    model = _import().make_tide(name)
    assert model.needs_radial_solve is expected


# =====================================================================================================================
# FixedQ: k_l(omega) = k_l * (1 - i / Q_l); -Im[k] = k_l / Q_l (frequency independent)
# =====================================================================================================================
@pytest.mark.parametrize("frequency", [1.0e-6, 1.0e-4, 3.0e-3])
def test_fixed_q_love(frequency):
    mod = _import()
    model = mod.make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})
    love = model.calc_love_numbers(2, frequency)
    assert isinstance(love, LoveNumbers)
    assert isclose(love.k.real, 0.3, rel_tol=1e-12)
    assert isclose(love.k.imag, -0.3 / 50.0, rel_tol=1e-12)
    # Analytic models do not produce displacement Love numbers.
    assert math.isnan(love.h.real) and math.isnan(love.l.real)
    assert isclose(model.calc_neg_imk(2, frequency), 0.3 / 50.0, rel_tol=1e-12)


def test_fixed_q_zero_q_is_elastic():
    """A zero/unset Q_l must not divide by zero; it should give no dissipation."""
    model = _import().make_tide("cpl", {"fixed_k": [0.3]})  # no fixed_q -> 0
    love = model.calc_love_numbers(2, 1.0e-5)
    assert isclose(love.k.real, 0.3, rel_tol=1e-12)
    assert love.k.imag == 0.0
    assert model.calc_neg_imk(2, 1.0e-5) == 0.0


# =====================================================================================================================
# FixedLag (CTL): -Im[k] = k_l * omega * dt_l
# =====================================================================================================================
@pytest.mark.parametrize("frequency", [1.0e-6, 1.0e-4, 3.0e-3])
def test_fixed_lag_love(frequency):
    model = _import().make_tide("ctl", {"fixed_k": [0.3], "fixed_dt": [100.0]})
    assert isclose(model.calc_neg_imk(2, frequency), 0.3 * frequency * 100.0, rel_tol=1e-12)


# =====================================================================================================================
# CTL_Q: -Im[k] = k_l * omega * dt_l / Q_l
# =====================================================================================================================
@pytest.mark.parametrize("frequency", [1.0e-6, 1.0e-4])
def test_ctlq_love(frequency):
    model = _import().make_tide("ctl_q", {"fixed_k": [0.3], "fixed_dt": [100.0], "fixed_q": [20.0]})
    expected = 0.3 * frequency * 100.0 / 20.0
    assert isclose(model.calc_neg_imk(2, frequency), expected, rel_tol=1e-12)


# =====================================================================================================================
# Rheology: passthrough of the supplied solver Love number
# =====================================================================================================================
def test_rheology_passthrough():
    """The rheology model returns the full radial-solver Love suite (k, h, l) unchanged."""
    model = _import().make_tide("rheology")
    supplied = LoveNumbers(k=0.25 - 0.013j, h=0.90 - 0.05j, l=0.30 - 0.02j)
    result = model.calc_love_numbers(2, 1.0e-5, supplied)
    assert result.k == supplied.k
    assert result.h == supplied.h
    assert result.l == supplied.l
    assert isclose(model.calc_neg_imk(2, 1.0e-5, supplied), 0.013, rel_tol=1e-12)


# =====================================================================================================================
# Per-degree parameters (l = 2 and l = 3 distinct)
# =====================================================================================================================
def test_per_degree_parameters():
    model = _import().make_tide("cpl", {"fixed_k": [0.30, 0.10], "fixed_q": [50.0, 80.0]})
    assert isclose(model.get_fixed_k(2), 0.30)
    assert isclose(model.get_fixed_k(3), 0.10)
    assert isclose(model.calc_neg_imk(2, 1.0e-5), 0.30 / 50.0, rel_tol=1e-12)
    assert isclose(model.calc_neg_imk(3, 1.0e-5), 0.10 / 80.0, rel_tol=1e-12)
    # A degree with no provided parameters contributes nothing.
    assert model.get_fixed_k(4) == 0.0
    assert model.calc_neg_imk(4, 1.0e-5) == 0.0


# =====================================================================================================================
# Config dict
# =====================================================================================================================
def test_config_dict_fixed_q():
    model = _import().make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})
    d = model.get_config_dict()
    assert d["model"] == "fixed_q"
    # Per-degree lists run l = 2..10 (length 9); index 0 == l=2.
    assert len(d["fixed_k"]) == 9
    assert isclose(d["fixed_k"][0], 0.3)
    assert isclose(d["fixed_q"][0], 50.0)


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
@pytest.mark.parametrize("name,config", [
    ("rheology", None),
    ("cpl", {"fixed_k": [0.3, 0.1], "fixed_q": [50.0, 80.0]}),
    ("ctl", {"fixed_k": [0.3], "fixed_dt": [120.0]}),
    ("ctl_q", {"fixed_k": [0.3, 0.05], "fixed_dt": [120.0, 90.0], "fixed_q": [40.0, 60.0]}),
])
def test_binary_round_trip(name, config):
    mod = _import()
    model = mod.make_tide(name, config)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, f"{name}.tpyb")
        model.save_binary(path)
        reloaded = mod.make_tide(name)
        reloaded.load_binary(path)
    solver_love = LoveNumbers(k=0.2 - 0.01j)
    assert reloaded.get_config_dict() == model.get_config_dict()
    assert reloaded.calc_neg_imk(2, 1.0e-4, solver_love) == model.calc_neg_imk(2, 1.0e-4, solver_love)


# =====================================================================================================================
# isinstance chain
# =====================================================================================================================
def test_isinstance_chain():
    from TidalPy.Utilities_x.classes_x.classes import PhysicsBase, TidalPyBaseClass
    mod = _import()
    model = mod.make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})
    assert isinstance(model, mod.FixedQTide)
    assert isinstance(model, mod.TideBase)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)
