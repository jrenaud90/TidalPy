"""The global (1D) tide model hierarchy: factory, Love-number laws, per-degree parameters, config, and binary I/O."""
import math
import os
import tempfile
from math import isclose

import pytest

import TidalPy
import TidalPy.Tides.classes as tide_classes
from TidalPy.Tides.love.love import LoveNumbers
from TidalPy.Utilities.classes.classes import PhysicsBase, TidalPyBaseClass


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
    ("constant_time_lag_and_q", "CTLQTide"),
    ("CTL_Q", "CTLQTide"),  # case-insensitive
])
def test_make_tide_returns_subclass(name, cls_attr):
    """Each name and alias builds its subclass."""
    model = tide_classes.make_tide(name)
    assert isinstance(model, getattr(tide_classes, cls_attr))
    assert isinstance(model, tide_classes.TideBase)


def test_make_tide_unknown_name_raises():
    with pytest.raises(ValueError):
        tide_classes.make_tide("not_a_model")


@pytest.mark.parametrize("name,expected", [
    ("rheology", True),
    ("cpl", False),
    ("ctl", False),
    ("ctl_q", False),
])
def test_needs_radial_solve(name, expected):
    """Only the rheology model needs a radial solve."""
    assert tide_classes.make_tide(name).needs_radial_solve is expected


@pytest.mark.parametrize("frequency", [1.0e-6, 1.0e-4, 3.0e-3])
def test_fixed_q_love(frequency):
    """Fixed Q: k = k_l (1 - i / Q_l) at every frequency, with no displacement Love numbers."""
    model = tide_classes.make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})
    love = model.calc_love_numbers(2, frequency)
    assert isinstance(love, LoveNumbers)
    assert isclose(love.k.real, 0.3, rel_tol=1e-12)
    assert isclose(love.k.imag, -0.3 / 50.0, rel_tol=1e-12)
    assert math.isnan(love.h.real) and math.isnan(love.l.real)
    assert isclose(model.calc_neg_imk(2, frequency), 0.3 / 50.0, rel_tol=1e-12)


def test_fixed_q_zero_q_is_elastic():
    """A zero Q_l gives no dissipation instead of dividing by zero."""
    model = tide_classes.make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [0.0]})
    love = model.calc_love_numbers(2, 1.0e-5)
    assert isclose(love.k.real, 0.3, rel_tol=1e-12)
    assert love.k.imag == 0.0
    assert model.calc_neg_imk(2, 1.0e-5) == 0.0


def test_partial_config_takes_the_configured_defaults():
    """A key left out of the config takes the [tides] value, as the world builder does."""
    configured_q = TidalPy.config["tides"]["fixed_q"][0]
    model = tide_classes.make_tide("cpl", {"fixed_k": [0.3]})
    assert isclose(model.calc_neg_imk(2, 1.0e-5), 0.3 / configured_q, rel_tol=1e-12)
    configured_k = TidalPy.config["tides"]["fixed_k"][0]
    model = tide_classes.make_tide("cpl", {"fixed_q": [50.0]})
    assert isclose(model.calc_neg_imk(2, 1.0e-5), configured_k / 50.0, rel_tol=1e-12)


@pytest.mark.parametrize("config", (
    {"fixed_q": [-50.0]},
    {"fixed_k": [float("nan")]},
    {"fixed_k": [0.3] * 10},
))
def test_bad_per_degree_parameters_are_rejected(config):
    with pytest.raises(ValueError):
        tide_classes.make_tide("cpl", config)


@pytest.mark.parametrize("frequency", [1.0e-6, 1.0e-4, 3.0e-3])
@pytest.mark.parametrize("name, config, fixed_q", [
    ("ctl", {"fixed_k": [0.3], "fixed_dt_s": [100.0]}, 1.0),
    ("ctl_q", {"fixed_k": [0.3], "fixed_dt_s": [100.0], "fixed_q": [20.0]}, 20.0),
])
def test_time_lag_love(name, config, fixed_q, frequency):
    """Time-lag models: -Im[k] = k_l omega dt_l, divided by Q_l for CTL_Q."""
    model = tide_classes.make_tide(name, config)
    assert isclose(model.calc_neg_imk(2, frequency), 0.3 * frequency * 100.0 / fixed_q, rel_tol=1e-12)


def test_rheology_passthrough():
    """The rheology model returns the supplied radial-solver Love numbers unchanged."""
    model = tide_classes.make_tide("rheology")
    supplied = LoveNumbers(k=0.25 - 0.013j, h=0.90 - 0.05j, l=0.30 - 0.02j)
    result = model.calc_love_numbers(2, 1.0e-5, supplied)
    assert result.k == supplied.k
    assert result.h == supplied.h
    assert result.l == supplied.l
    assert isclose(model.calc_neg_imk(2, 1.0e-5, supplied), 0.013, rel_tol=1e-12)


def test_per_degree_parameters():
    """Degrees 2 and 3 use their own parameters; an unset degree contributes nothing."""
    model = tide_classes.make_tide("cpl", {"fixed_k": [0.30, 0.10], "fixed_q": [50.0, 80.0]})
    assert isclose(model.get_fixed_k(2), 0.30)
    assert isclose(model.get_fixed_k(3), 0.10)
    assert isclose(model.calc_neg_imk(2, 1.0e-5), 0.30 / 50.0, rel_tol=1e-12)
    assert isclose(model.calc_neg_imk(3, 1.0e-5), 0.10 / 80.0, rel_tol=1e-12)
    assert model.get_fixed_k(4) == 0.0
    assert model.calc_neg_imk(4, 1.0e-5) == 0.0


@pytest.mark.parametrize("name, config, fixed_q, fixed_dt", [
    ("rheology", None, math.nan, math.nan),
    ("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}, 50.0, math.nan),
    ("ctl", {"fixed_k": [0.3], "fixed_dt_s": [100.0]}, math.nan, 100.0),
    ("ctl_q", {"fixed_k": [0.3], "fixed_dt_s": [100.0], "fixed_q": [20.0]}, 20.0, 100.0),
])
def test_fixed_q_and_dt_on_every_model(name, config, fixed_q, fixed_dt):
    """Every tide model answers get_fixed_q and get_fixed_dt, NaN for a parameter it does not carry."""
    model = tide_classes.make_tide(name, config)
    for found, expected in ((model.get_fixed_q(2), fixed_q), (model.get_fixed_dt(2), fixed_dt)):
        if math.isnan(expected):
            assert math.isnan(found)
        else:
            assert isclose(found, expected)


def test_config_dict_fixed_q():
    """The config dict names the model and holds its per-degree lists as given, from l = 2."""
    config = tide_classes.make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}).get_config_dict()
    assert config["model"] == "fixed_q"
    assert config["fixed_k"] == [0.3]
    assert isclose(config["fixed_k"][0], 0.3)
    assert isclose(config["fixed_q"][0], 50.0)


@pytest.mark.parametrize("name,config", [
    ("rheology", None),
    ("cpl", {"fixed_k": [0.3, 0.1], "fixed_q": [50.0, 80.0]}),
    ("ctl", {"fixed_k": [0.3], "fixed_dt_s": [120.0]}),
    ("ctl_q", {"fixed_k": [0.3, 0.05], "fixed_dt_s": [120.0, 90.0], "fixed_q": [40.0, 60.0]}),
])
def test_binary_round_trip(name, config):
    """A saved and reloaded model has the same config and -Im[k]."""
    model = tide_classes.make_tide(name, config)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, f"{name}.tpyb")
        model.save_binary(path)
        reloaded = tide_classes.make_tide(name)
        reloaded.load_binary(path)
    solver_love = LoveNumbers(k=0.2 - 0.01j)
    assert reloaded.get_config_dict() == model.get_config_dict()
    assert reloaded.calc_neg_imk(2, 1.0e-4, solver_love) == model.calc_neg_imk(2, 1.0e-4, solver_love)


def test_isinstance_chain():
    """A tide model is a TideBase, PhysicsBase, and TidalPyBaseClass."""
    model = tide_classes.make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})
    assert isinstance(model, tide_classes.FixedQTide)
    assert isinstance(model, tide_classes.TideBase)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)


@pytest.mark.parametrize("model_name", ["cpl", "ctl", "ctl_q"])
def test_the_lag_depends_on_the_frequency_magnitude(model_name):
    """The analytic models give a negative frequency the lag of its magnitude, as the collapse (which passes |omega|)
    and the cpl and ctl Love methods do."""
    from TidalPy.Tides.classes.tide import make_tide, tide_config_keys
    lists = {"fixed_k": [0.3], "fixed_q": [50.0], "fixed_dt_s": [600.0]}
    model = make_tide(model_name, {key: lists[key] for key in tide_config_keys(model_name)})
    frequency = 4.1e-5
    assert model.calc_neg_imk(2, -frequency) == model.calc_neg_imk(2, frequency) > 0.0
