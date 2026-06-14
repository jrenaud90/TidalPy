"""
Tests for TidalPy.partial_melt_x.partial_melt — partial-melt model hierarchy.

Covers the Off / Spohn / Henning models: melt-fraction formula, post-melt
viscosity/shear (cross-checked against the legacy melting_models formulas), the
name factory (incl. unknown-name ValueError and rich-subclass return), config
dicts, binary round-trips, and the isinstance chain.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import math
import os
import tempfile

import pytest


def _import():
    from TidalPy.partial_melt_x import partial_melt as _mod
    return _mod


# Reference material (rock-like, MKS)
_SOLIDUS   = 1600.0
_LIQUIDUS  = 2000.0
_LIQ_SHEAR = 1.0e-5
_PREMELT_VISC  = 1.0e22
_PREMELT_SHEAR = 6.0e10
_LIQ_VISC      = 0.2


def _melt_fraction(T, solidus=_SOLIDUS, liquidus=_LIQUIDUS):
    return min(1.0, max(0.0, (T - solidus) / (liquidus - solidus)))


# =====================================================================================================================
# Melt fraction
# =====================================================================================================================
@pytest.mark.parametrize("T", [1400.0, 1600.0, 1700.0, 1800.0, 2000.0, 2200.0])
def test_melt_fraction_formula(T):
    m = _import().OffPartialMelt(solidus_k=_SOLIDUS, liquidus_k=_LIQUIDUS)
    assert m.calc_melt_fraction(T) == pytest.approx(_melt_fraction(T))


def test_melt_fraction_degenerate_envelope():
    # solidus >= liquidus -> fully solid (phi = 0).
    m = _import().OffPartialMelt(solidus_k=2000.0, liquidus_k=2000.0)
    assert m.calc_melt_fraction(1900.0) == 0.0


# =====================================================================================================================
# Off model
# =====================================================================================================================
def test_off_returns_premelt():
    m = _import().OffPartialMelt(solidus_k=_SOLIDUS, liquidus_k=_LIQUIDUS, liquid_shear_pa=_LIQ_SHEAR)
    phi, visc, shear = m.calc_partial_melt(1800.0, _PREMELT_VISC, _PREMELT_SHEAR, _LIQ_VISC)
    assert phi == pytest.approx(_melt_fraction(1800.0))
    assert visc == pytest.approx(_PREMELT_VISC)
    assert shear == pytest.approx(_PREMELT_SHEAR)


# =====================================================================================================================
# Spohn model
# =====================================================================================================================
@pytest.mark.parametrize("T", [1700.0, 1900.0, 2100.0])
def test_spohn_formula(T):
    m = _import().SpohnPartialMelt(solidus_k=_SOLIDUS, liquidus_k=_LIQUIDUS, liquid_shear_pa=_LIQ_SHEAR)
    phi, visc, shear = m.calc_partial_melt(T, _PREMELT_VISC, _PREMELT_SHEAR, _LIQ_VISC)
    exp_visc  = max(_LIQ_VISC,  10.0 ** ((27000.0 / T) - 1.0))
    exp_shear = max(_LIQ_SHEAR, 10.0 ** ((82000.0 / T) - 40.6))
    assert visc == pytest.approx(exp_visc, rel=1e-9)
    assert shear == pytest.approx(exp_shear, rel=1e-9)


# =====================================================================================================================
# Henning model (three regimes) — cross-check the legacy formula
# =====================================================================================================================
def _henning_expected(T):
    crit, width = 0.5, 0.05
    vslope1, vfall = 13.5, 370.0
    sp1, sp2, sfall = 40000.0, 25.0, 700.0
    phi = _melt_fraction(T)
    crit_plus = crit + width
    break_temp = _SOLIDUS + crit * (_LIQUIDUS - _SOLIDUS)
    if phi <= 0.0:
        visc, shear = _PREMELT_VISC, _PREMELT_SHEAR
    elif phi < crit:
        visc = _PREMELT_VISC * math.exp(-vslope1 * phi)
        shear = _PREMELT_SHEAR * math.exp((sp1 / T) - sp2)
    elif phi <= crit_plus:
        visc = _PREMELT_VISC * math.exp(-vslope1 * crit) * math.exp(-vfall * (phi - crit))
        shear = _PREMELT_SHEAR * math.exp((sp1 / break_temp) - sp2) * math.exp(-sfall * (phi - crit))
    else:
        visc, shear = _LIQ_VISC, _LIQ_SHEAR
    return max(_LIQ_VISC, visc), max(_LIQ_SHEAR, shear)


@pytest.mark.parametrize("T", [1500.0, 1700.0, 1810.0, 1900.0])
def test_henning_regimes(T):
    m = _import().HenningPartialMelt(solidus_k=_SOLIDUS, liquidus_k=_LIQUIDUS, liquid_shear_pa=_LIQ_SHEAR)
    phi, visc, shear = m.calc_partial_melt(T, _PREMELT_VISC, _PREMELT_SHEAR, _LIQ_VISC)
    exp_visc, exp_shear = _henning_expected(T)
    assert visc == pytest.approx(exp_visc, rel=1e-9)
    assert shear == pytest.approx(exp_shear, rel=1e-9)


def test_henning_liquid_regime_floors():
    m = _import().HenningPartialMelt(solidus_k=_SOLIDUS, liquidus_k=_LIQUIDUS, liquid_shear_pa=_LIQ_SHEAR)
    _, visc, shear = m.calc_partial_melt(1950.0, _PREMELT_VISC, _PREMELT_SHEAR, _LIQ_VISC)
    assert visc == pytest.approx(_LIQ_VISC)
    assert shear == pytest.approx(_LIQ_SHEAR)


def test_henning_weakens_with_temperature():
    m = _import().HenningPartialMelt(solidus_k=_SOLIDUS, liquidus_k=_LIQUIDUS, liquid_shear_pa=_LIQ_SHEAR)
    _, v1, s1 = m.calc_partial_melt(1650.0, _PREMELT_VISC, _PREMELT_SHEAR, _LIQ_VISC)
    _, v2, s2 = m.calc_partial_melt(1750.0, _PREMELT_VISC, _PREMELT_SHEAR, _LIQ_VISC)
    assert v2 < v1
    assert s2 < s1


# =====================================================================================================================
# Factory
# =====================================================================================================================
@pytest.mark.parametrize("name,cls_attr", [
    ("off", "OffPartialMelt"),
    ("none", "OffPartialMelt"),
    ("spohn", "SpohnPartialMelt"),
    ("fischer", "SpohnPartialMelt"),
    ("FISCHER_SPOHN", "SpohnPartialMelt"),
    ("henning", "HenningPartialMelt"),
])
def test_factory_names(name, cls_attr):
    mod = _import()
    m = mod.make_partial_melt(name)
    assert isinstance(m, getattr(mod, cls_attr))


def test_factory_unknown_raises():
    with pytest.raises(ValueError):
        _import().make_partial_melt("not_a_model")


def test_factory_config_override():
    m = _import().make_partial_melt("henning", {"solidus_k": 1500.0, "crit_melt_frac": 0.4})
    d = m.get_config_dict()
    assert d["solidus_k"] == pytest.approx(1500.0)
    assert d["crit_melt_frac"] == pytest.approx(0.4)


# =====================================================================================================================
# Config dict
# =====================================================================================================================
def test_config_dict_keys():
    m = _import().HenningPartialMelt()
    d = m.get_config_dict()
    for key in ("model", "solidus_k", "liquidus_k", "liquid_shear_pa",
                "crit_melt_frac", "hn_visc_slope_1", "hn_shear_param_1"):
        assert key in d
    assert d["model"] == "henning"


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
@pytest.mark.parametrize("name", ["off", "spohn", "henning"])
def test_binary_round_trip(name):
    mod = _import()
    m = mod.make_partial_melt(name, {"solidus_k": 1550.0, "liquidus_k": 1950.0})
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, f"{name}.tpyb")
        m.save_binary(path)
        reloaded = mod.make_partial_melt(name)
        reloaded.load_binary(path)
    assert reloaded.get_config_dict() == m.get_config_dict()
    # A representative evaluation survives the round-trip.
    assert reloaded.calc_partial_melt(1700.0, _PREMELT_VISC, _PREMELT_SHEAR, _LIQ_VISC) == \
        m.calc_partial_melt(1700.0, _PREMELT_VISC, _PREMELT_SHEAR, _LIQ_VISC)


# =====================================================================================================================
# isinstance chain
# =====================================================================================================================
def test_isinstance_chain():
    from TidalPy.Utilities_x.classes_x.classes import PhysicsBase, TidalPyBaseClass
    mod = _import()
    m = mod.HenningPartialMelt()
    assert isinstance(m, mod.PartialMeltBase)
    assert isinstance(m, PhysicsBase)
    assert isinstance(m, TidalPyBaseClass)
