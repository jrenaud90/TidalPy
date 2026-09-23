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
    m = _import().OffPartialMelt(solidus=_SOLIDUS, liquidus=_LIQUIDUS)
    assert m.calc_melt_fraction(T) == pytest.approx(_melt_fraction(T))


def test_melt_fraction_degenerate_envelope():
    # solidus >= liquidus -> fully solid (phi = 0).
    m = _import().OffPartialMelt(solidus=2000.0, liquidus=2000.0)
    assert m.calc_melt_fraction(1900.0) == 0.0


@pytest.mark.parametrize("cls_name", ["OffPartialMelt", "SpohnPartialMelt", "HenningPartialMelt"])
def test_non_finite_temperature_gives_nan(cls_name):
    """A non-finite temperature has no melt state: NaN out (Off keeps its unweakened strengths)."""
    m = getattr(_import(), cls_name)(solidus=_SOLIDUS, liquidus=_LIQUIDUS)
    phi, visc, shear = m.calc_partial_melt(math.nan, _PREMELT_VISC, _PREMELT_SHEAR)
    assert math.isnan(phi)
    if cls_name == "OffPartialMelt":
        assert (visc, shear) == (_PREMELT_VISC, _PREMELT_SHEAR)
    else:
        assert math.isnan(visc) and math.isnan(shear)


# =====================================================================================================================
# Off model
# =====================================================================================================================
def test_off_returns_premelt():
    m = _import().OffPartialMelt(solidus=_SOLIDUS, liquidus=_LIQUIDUS, liquid_shear=_LIQ_SHEAR)
    phi, visc, shear = m.calc_partial_melt(1800.0, _PREMELT_VISC, _PREMELT_SHEAR)
    assert phi == pytest.approx(_melt_fraction(1800.0))
    assert visc == pytest.approx(_PREMELT_VISC)
    assert shear == pytest.approx(_PREMELT_SHEAR)


# =====================================================================================================================
# Spohn model
# =====================================================================================================================
@pytest.mark.parametrize("T", [1700.0, 1900.0, 2100.0])
def test_spohn_formula(T):
    m = _import().SpohnPartialMelt(solidus=_SOLIDUS, liquidus=_LIQUIDUS, liquid_shear=_LIQ_SHEAR)
    phi, visc, shear = m.calc_partial_melt(T, _PREMELT_VISC, _PREMELT_SHEAR)
    exp_visc  = max(_LIQ_VISC,  10.0 ** ((27000.0 / T) - 1.0))
    exp_shear = max(_LIQ_SHEAR, 10.0 ** ((82000.0 / T) - 40.6))
    assert visc == pytest.approx(exp_visc, rel=1e-9)
    assert shear == pytest.approx(exp_shear, rel=1e-9)


@pytest.mark.parametrize("T", [200.0, 1000.0, _SOLIDUS])
def test_spohn_below_solidus_returns_premelt(T):
    """The Fischer-Spohn law applies only above the solidus; below it (where the law would overflow) nothing melts."""
    m = _import().SpohnPartialMelt(solidus=_SOLIDUS, liquidus=_LIQUIDUS, liquid_shear=_LIQ_SHEAR)
    phi, visc, shear = m.calc_partial_melt(T, _PREMELT_VISC, _PREMELT_SHEAR)
    assert phi == 0.0
    assert (visc, shear) == (_PREMELT_VISC, _PREMELT_SHEAR)


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
    m = _import().HenningPartialMelt(solidus=_SOLIDUS, liquidus=_LIQUIDUS, liquid_shear=_LIQ_SHEAR)
    phi, visc, shear = m.calc_partial_melt(T, _PREMELT_VISC, _PREMELT_SHEAR)
    exp_visc, exp_shear = _henning_expected(T)
    assert visc == pytest.approx(exp_visc, rel=1e-9)
    assert shear == pytest.approx(exp_shear, rel=1e-9)


def test_henning_liquid_regime_floors():
    m = _import().HenningPartialMelt(solidus=_SOLIDUS, liquidus=_LIQUIDUS, liquid_shear=_LIQ_SHEAR)
    _, visc, shear = m.calc_partial_melt(1950.0, _PREMELT_VISC, _PREMELT_SHEAR)
    assert visc == pytest.approx(_LIQ_VISC)
    assert shear == pytest.approx(_LIQ_SHEAR)


def test_liquid_viscosity_is_a_model_parameter():
    """The liquid viscosity is the model's own, reported and applied past breakdown and as the floor."""
    m = _import().HenningPartialMelt(solidus=_SOLIDUS, liquidus=_LIQUIDUS, liquid_viscosity=3.0)
    assert m.liquid_viscosity == 3.0
    assert m.get_config_dict()["liquid_viscosity_pas"] == 3.0
    _, visc, _ = m.calc_partial_melt(1950.0, _PREMELT_VISC, _PREMELT_SHEAR)
    assert visc == 3.0
    # Sub-critical weakening lowers the viscosity below the pre-melt value.
    _, visc_partial, _ = m.calc_partial_melt(1700.0, _PREMELT_VISC, _PREMELT_SHEAR)
    assert visc_partial == pytest.approx(_PREMELT_VISC * math.exp(-13.5 * 0.25), rel=1e-12)


# =====================================================================================================================
# Bulk-modulus weakening
# =====================================================================================================================
_PREMELT_BULK = 1.3e11
_LIQ_BULK     = 2.0e10


def _hashin_shtrikman(phi, framework_shear):
    return _PREMELT_BULK + phi / (1.0 / (_LIQ_BULK - _PREMELT_BULK)
                                  + (1.0 - phi) / (_PREMELT_BULK + (4.0 / 3.0) * framework_shear))


def test_bulk_weakening_is_off_by_default():
    m = _import().make_partial_melt("henning")
    assert m.bulk_melt_weakening is False
    assert m.calc_bulk_modulus_melt(1900.0, _PREMELT_BULK, 0.0) == _PREMELT_BULK


@pytest.mark.parametrize("T", [1500.0, 1640.0, 1800.0, 2000.0, 2100.0])
def test_bulk_weakening_follows_hashin_shtrikman(T):
    m = _import().HenningPartialMelt(
        solidus=_SOLIDUS, liquidus=_LIQUIDUS, bulk_melt_weakening=True, liquid_bulk_modulus=_LIQ_BULK)
    phi, _, shear = m.calc_partial_melt(T, _PREMELT_VISC, _PREMELT_SHEAR)
    bulk = m.calc_bulk_modulus_melt(T, _PREMELT_BULK, shear)
    if phi == 0.0:
        assert bulk == _PREMELT_BULK
    else:
        assert bulk == pytest.approx(_hashin_shtrikman(phi, shear), rel=1e-12)
    assert _LIQ_BULK <= bulk <= _PREMELT_BULK


def test_bulk_weakening_limits():
    """Weak while the framework holds, the Reuss average once it has collapsed, the melt's own value at phi = 1."""
    m = _import().OffPartialMelt(
        solidus=_SOLIDUS, liquidus=_LIQUIDUS, bulk_melt_weakening=True, liquid_bulk_modulus=_LIQ_BULK)
    phi = 0.1
    temperature = _SOLIDUS + phi * (_LIQUIDUS - _SOLIDUS)
    framework = m.calc_bulk_modulus_melt(temperature, _PREMELT_BULK, _PREMELT_SHEAR)
    reuss = 1.0 / ((1.0 - phi) / _PREMELT_BULK + phi / _LIQ_BULK)
    assert m.calc_bulk_modulus_melt(temperature, _PREMELT_BULK, 0.0) == pytest.approx(reuss, rel=1e-12)
    assert reuss < framework < _PREMELT_BULK
    # Melt affects the bulk modulus far less than Henning affects the shear modulus.
    assert framework / _PREMELT_BULK > 0.8
    assert m.calc_bulk_modulus_melt(_LIQUIDUS, _PREMELT_BULK, 0.0) == pytest.approx(_LIQ_BULK, rel=1e-12)


def test_henning_weakens_with_temperature():
    m = _import().HenningPartialMelt(solidus=_SOLIDUS, liquidus=_LIQUIDUS, liquid_shear=_LIQ_SHEAR)
    _, v1, s1 = m.calc_partial_melt(1650.0, _PREMELT_VISC, _PREMELT_SHEAR)
    _, v2, s2 = m.calc_partial_melt(1750.0, _PREMELT_VISC, _PREMELT_SHEAR)
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
    for key in ("model", "solidus_k", "liquidus_k", "liquid_shear_pa", "liquid_viscosity_pas",
                "bulk_melt_weakening", "liquid_bulk_modulus_pa",
                "crit_melt_frac", "hn_visc_slope_1", "hn_shear_param_1_k"):
        assert key in d
    assert d["model"] == "henning"


# Constructor keyword (and property) -> config key, for the parameters whose config key carries a unit.
_CONFIG_KEYS = {
    "fs_visc_power_slope":  "fs_visc_power_slope_k",
    "fs_shear_power_slope": "fs_shear_power_slope_k",
    "hn_shear_param_1":     "hn_shear_param_1_k",
}


@pytest.mark.parametrize("cls_name,params", [
    ("SpohnPartialMelt", dict(fs_visc_power_slope=25000.0, fs_visc_power_phase=1.5,
                              fs_shear_power_slope=80000.0, fs_shear_power_phase=39.0)),
    ("HenningPartialMelt", dict(crit_melt_frac=0.4, crit_melt_frac_width=0.08, hn_visc_slope_1=12.0,
                                hn_visc_falloff_slope=350.0, hn_shear_param_1=41000.0, hn_shear_param_2=24.0,
                                hn_shear_falloff_slope=650.0)),
])
def test_model_parameter_properties(cls_name, params):
    """Each model-specific parameter is a property named like its constructor keyword, emitted under its config key."""
    m = getattr(_import(), cls_name)(solidus=_SOLIDUS, liquidus=_LIQUIDUS, liquid_shear=_LIQ_SHEAR, **params)
    config = m.get_config_dict()
    for name, value in params.items():
        assert getattr(m, name) == value, name
        assert config[_CONFIG_KEYS.get(name, name)] == value, name


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
@pytest.mark.parametrize("name", ["off", "spohn", "henning"])
def test_binary_round_trip(name):
    mod = _import()
    m = mod.make_partial_melt(name, {"solidus_k": 1550.0, "liquidus_k": 1950.0, "liquid_viscosity_pas": 0.7,
                                     "bulk_melt_weakening": True, "liquid_bulk_modulus_pa": 1.5e10})
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, f"{name}.tpyb")
        m.save_binary(path)
        reloaded = mod.make_partial_melt(name)
        reloaded.load_binary(path)
    assert reloaded.get_config_dict() == m.get_config_dict()
    assert reloaded.bulk_melt_weakening is True
    assert reloaded.calc_bulk_modulus_melt(1700.0, 1.3e11, 1.0e9) == m.calc_bulk_modulus_melt(1700.0, 1.3e11, 1.0e9)
    # A representative evaluation survives the round-trip.
    assert reloaded.calc_partial_melt(1700.0, _PREMELT_VISC, _PREMELT_SHEAR) == \
        m.calc_partial_melt(1700.0, _PREMELT_VISC, _PREMELT_SHEAR)


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
