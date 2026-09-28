"""Partial-melt models (off, Spohn, Henning): melt fraction, weakening laws, bulk weakening, factory, and I/O."""

import math

import pytest

from TidalPy.PartialMelt import partial_melt
from TidalPy.Utilities.classes.classes import PhysicsBase, TidalPyBaseClass

_SOLIDUS = 1600.0
_LIQUIDUS = 2000.0
_LIQ_SHEAR = 1.0e-5
_PREMELT_VISC = 1.0e22
_PREMELT_SHEAR = 6.0e10
_LIQ_VISC = 0.2
_PREMELT_BULK = 1.3e11
_LIQ_BULK = 2.0e10

_CLASS_NAMES = ["OffPartialMelt", "SpohnPartialMelt", "HenningPartialMelt"]


def _melt_fraction(T, solidus=_SOLIDUS, liquidus=_LIQUIDUS):
    return min(1.0, max(0.0, (T - solidus) / (liquidus - solidus)))


def _model(cls_name, **kwargs):
    kwargs.setdefault("solidus", _SOLIDUS)
    kwargs.setdefault("liquidus", _LIQUIDUS)
    return getattr(partial_melt, cls_name)(**kwargs)


def _melt(model, T):
    """calc_partial_melt at the reference pre-melt viscosity and shear modulus."""
    return model.calc_partial_melt(T, _PREMELT_VISC, _PREMELT_SHEAR)


# =====================================================================================================================
# Melt fraction and Off model
# =====================================================================================================================
@pytest.mark.parametrize("T", [1400.0, 1600.0, 1700.0, 1800.0, 2000.0, 2200.0])
def test_melt_fraction_formula(T):
    assert _model("OffPartialMelt").calc_melt_fraction(T) == pytest.approx(_melt_fraction(T))


def test_melt_fraction_degenerate_envelope():
    """A solidus at or above the liquidus is fully solid."""
    assert _model("OffPartialMelt", solidus=2000.0, liquidus=2000.0).calc_melt_fraction(1900.0) == 0.0


@pytest.mark.parametrize("cls_name", _CLASS_NAMES)
def test_non_finite_temperature_gives_nan(cls_name):
    """A non-finite temperature gives NaN (Off keeps its unweakened strengths)."""
    phi, visc, shear = _melt(_model(cls_name), math.nan)
    assert math.isnan(phi)
    if cls_name == "OffPartialMelt":
        assert (visc, shear) == (_PREMELT_VISC, _PREMELT_SHEAR)
    else:
        assert math.isnan(visc) and math.isnan(shear)


def test_off_returns_premelt():
    phi, visc, shear = _melt(_model("OffPartialMelt", liquid_shear=_LIQ_SHEAR), 1800.0)
    assert phi == pytest.approx(_melt_fraction(1800.0))
    assert visc == pytest.approx(_PREMELT_VISC)
    assert shear == pytest.approx(_PREMELT_SHEAR)


# =====================================================================================================================
# Spohn model
# =====================================================================================================================
@pytest.mark.parametrize("T", [1700.0, 1900.0, 2100.0])
def test_spohn_formula(T):
    _, visc, shear = _melt(_model("SpohnPartialMelt", liquid_shear=_LIQ_SHEAR), T)
    assert visc == pytest.approx(max(_LIQ_VISC, 10.0 ** ((27000.0 / T) - 1.0)), rel=1e-9)
    assert shear == pytest.approx(max(_LIQ_SHEAR, 10.0 ** ((82000.0 / T) - 40.6)), rel=1e-9)


def test_spohn_anchored_at_the_solidus():
    """At another solidus the law keeps its solidus value and slope instead of the absolute 1600 K fit."""
    solidus, T = 250.0, 252.0
    model = _model("SpohnPartialMelt", solidus=solidus, liquidus=300.0, liquid_shear=_LIQ_SHEAR)
    phi, visc, shear = _melt(model, T)
    assert phi == pytest.approx(0.04)
    assert visc == pytest.approx(10.0 ** (15.875 + 27000.0 * (1.0 / T - 1.0 / solidus)), rel=1e-9)
    assert shear == pytest.approx(10.0 ** (10.65 + 82000.0 * (1.0 / T - 1.0 / solidus)), rel=1e-9)
    # The published fit would give 10^263 Pa here.
    assert shear < 1.0e11


@pytest.mark.parametrize("T", [200.0, 1000.0, _SOLIDUS])
def test_spohn_below_solidus_returns_premelt(T):
    """Below the solidus (where the law would overflow) nothing melts."""
    phi, visc, shear = _melt(_model("SpohnPartialMelt", liquid_shear=_LIQ_SHEAR), T)
    assert phi == 0.0
    assert (visc, shear) == (_PREMELT_VISC, _PREMELT_SHEAR)


# =====================================================================================================================
# Henning model
# =====================================================================================================================
def _henning_expected(T):
    """The TidalPy 0.7 three-regime Henning formula."""
    crit, width = 0.5, 0.05
    visc_slope, visc_fall = 13.5, 370.0
    shear_param_1, shear_param_2, shear_fall = 40000.0, 25.0, 700.0
    phi = _melt_fraction(T)
    break_temp = _SOLIDUS + crit * (_LIQUIDUS - _SOLIDUS)
    if phi <= 0.0:
        visc, shear = _PREMELT_VISC, _PREMELT_SHEAR
    elif phi < crit:
        visc = _PREMELT_VISC * math.exp(-visc_slope * phi)
        # Henning et al. (2009) Eq. 20 as published: 25 = 40000 / 1600, the reference solidus.
        shear = _PREMELT_SHEAR * math.exp((shear_param_1 / T) - shear_param_2)
    elif phi <= crit + width:
        visc = _PREMELT_VISC * math.exp(-visc_slope * crit) * math.exp(-visc_fall * (phi - crit))
        shear = (_PREMELT_SHEAR * math.exp((shear_param_1 / break_temp) - shear_param_2)
                 * math.exp(-shear_fall * (phi - crit)))
    else:
        visc, shear = _LIQ_VISC, _LIQ_SHEAR
    return max(_LIQ_VISC, visc), max(_LIQ_SHEAR, shear)


@pytest.mark.parametrize("T", [1500.0, 1700.0, 1810.0, 1900.0, 1950.0])
def test_henning_regimes(T):
    """Each regime, including the liquid floors at 1950 K, matches the reference formula."""
    _, visc, shear = _melt(_model("HenningPartialMelt", liquid_shear=_LIQ_SHEAR), T)
    expected_visc, expected_shear = _henning_expected(T)
    assert visc == pytest.approx(expected_visc, rel=1e-9)
    assert shear == pytest.approx(expected_shear, rel=1e-9)


@pytest.mark.parametrize("solidus", [250.0, 1200.0, 1600.0])
def test_henning_shear_continuous_at_any_solidus(solidus):
    """The sub-critical shear law exp[b1 (1/T - 1/T_sol)] is 1 at the solidus, whatever the solidus."""
    liquidus = 1.25 * solidus
    model = _model("HenningPartialMelt", solidus=solidus, liquidus=liquidus, liquid_shear=_LIQ_SHEAR)
    _, _, shear = _melt(model, solidus * (1.0 + 1.0e-9))
    assert shear == pytest.approx(_PREMELT_SHEAR, rel=1e-4)
    T = solidus + 0.2 * (liquidus - solidus)
    _, _, shear = _melt(model, T)
    assert shear == pytest.approx(_PREMELT_SHEAR * math.exp(40000.0 * (1.0 / T - 1.0 / solidus)), rel=1e-9)


def test_henning_weakens_with_temperature():
    model = _model("HenningPartialMelt", liquid_shear=_LIQ_SHEAR)
    _, visc_cool, shear_cool = _melt(model, 1650.0)
    _, visc_hot, shear_hot = _melt(model, 1750.0)
    assert visc_hot < visc_cool
    assert shear_hot < shear_cool


def test_henning_retired_offset_key_warns_and_is_ignored():
    """An old config with hn_shear_param_2 builds with a warning, and the key has no effect."""
    with pytest.warns(UserWarning, match="hn_shear_param_2"):
        model = partial_melt.make_partial_melt("henning", {"solidus_k": 1200.0, "hn_shear_param_2": 3.0})
    assert "hn_shear_param_2" not in model.get_config_dict()
    reference = partial_melt.make_partial_melt("henning", {"solidus_k": 1200.0})
    assert _melt(model, 1300.0) == _melt(reference, 1300.0)


def test_liquid_viscosity_is_a_model_parameter():
    """The model's liquid viscosity is reported and applied past breakdown."""
    model = _model("HenningPartialMelt", liquid_viscosity=3.0)
    assert model.liquid_viscosity == 3.0
    assert model.get_config_dict()["liquid_viscosity_pas"] == 3.0
    _, visc, _ = _melt(model, 1950.0)
    assert visc == 3.0
    _, visc_partial, _ = _melt(model, 1700.0)
    assert visc_partial == pytest.approx(_PREMELT_VISC * math.exp(-13.5 * 0.25), rel=1e-12)


# =====================================================================================================================
# Bulk-modulus weakening
# =====================================================================================================================
def _hashin_shtrikman(phi, framework_shear):
    return _PREMELT_BULK + phi / (1.0 / (_LIQ_BULK - _PREMELT_BULK)
                                  + (1.0 - phi) / (_PREMELT_BULK + (4.0 / 3.0) * framework_shear))


def test_bulk_weakening_is_off_by_default():
    model = partial_melt.make_partial_melt("henning")
    assert model.bulk_melt_weakening is False
    assert model.calc_bulk_modulus_melt(1900.0, _PREMELT_BULK, 0.0) == _PREMELT_BULK


@pytest.mark.parametrize("T", [1500.0, 1640.0, 1800.0, 2000.0, 2100.0])
def test_bulk_weakening_follows_hashin_shtrikman(T):
    model = _model("HenningPartialMelt", bulk_melt_weakening=True, liquid_bulk_modulus=_LIQ_BULK)
    phi, _, shear = _melt(model, T)
    bulk = model.calc_bulk_modulus_melt(T, _PREMELT_BULK, shear)
    if phi == 0.0:
        assert bulk == _PREMELT_BULK
    else:
        assert bulk == pytest.approx(_hashin_shtrikman(phi, shear), rel=1e-12)
    assert _LIQ_BULK <= bulk <= _PREMELT_BULK


def test_bulk_weakening_limits():
    """Weak while the framework holds, the Reuss average once it collapses, the melt's value at phi = 1."""
    model = _model("OffPartialMelt", bulk_melt_weakening=True, liquid_bulk_modulus=_LIQ_BULK)
    phi = 0.1
    temperature = _SOLIDUS + phi * (_LIQUIDUS - _SOLIDUS)
    framework = model.calc_bulk_modulus_melt(temperature, _PREMELT_BULK, _PREMELT_SHEAR)
    reuss = 1.0 / ((1.0 - phi) / _PREMELT_BULK + phi / _LIQ_BULK)
    assert model.calc_bulk_modulus_melt(temperature, _PREMELT_BULK, 0.0) == pytest.approx(reuss, rel=1e-12)
    assert reuss < framework < _PREMELT_BULK
    # Melt weakens the bulk modulus far less than Henning weakens the shear modulus.
    assert framework / _PREMELT_BULK > 0.8
    assert model.calc_bulk_modulus_melt(_LIQUIDUS, _PREMELT_BULK, 0.0) == pytest.approx(_LIQ_BULK, rel=1e-12)


# =====================================================================================================================
# Factory, config dict, and binary round trip
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
    assert isinstance(partial_melt.make_partial_melt(name), getattr(partial_melt, cls_attr))


def test_factory_unknown_raises():
    with pytest.raises(ValueError):
        partial_melt.make_partial_melt("not_a_model")


def test_factory_config_override():
    config = partial_melt.make_partial_melt("henning", {"solidus_k": 1500.0, "crit_melt_frac": 0.4}).get_config_dict()
    assert config["solidus_k"] == pytest.approx(1500.0)
    assert config["crit_melt_frac"] == pytest.approx(0.4)


def test_config_dict_keys():
    config = partial_melt.HenningPartialMelt().get_config_dict()
    assert {"model", "solidus_k", "liquidus_k", "liquid_shear_pa", "liquid_viscosity_pas", "bulk_melt_weakening",
            "liquid_bulk_modulus_pa", "crit_melt_frac", "hn_visc_slope_1", "hn_shear_param_1_k"} <= set(config)
    assert config["model"] == "henning"


# Constructor keyword (and property) to config key, where the config key carries a unit.
_CONFIG_KEYS = {
    "fs_visc_power_slope":  "fs_visc_power_slope_k",
    "fs_shear_power_slope": "fs_shear_power_slope_k",
    "hn_shear_param_1":     "hn_shear_param_1_k",
}


@pytest.mark.parametrize("cls_name,params", [
    ("SpohnPartialMelt", dict(fs_visc_power_slope=25000.0, fs_visc_log10_at_solidus=15.5,
                              fs_shear_power_slope=80000.0, fs_shear_log10_at_solidus=10.2)),
    ("HenningPartialMelt", dict(crit_melt_frac=0.4, crit_melt_frac_width=0.08, hn_visc_slope_1=12.0,
                                hn_visc_falloff_slope=350.0, hn_shear_param_1=41000.0,
                                hn_shear_falloff_slope=650.0)),
])
def test_model_parameter_properties(cls_name, params):
    """Each model parameter is a property named like its constructor keyword, emitted under its config key."""
    model = _model(cls_name, liquid_shear=_LIQ_SHEAR, **params)
    config = model.get_config_dict()
    for name, value in params.items():
        assert getattr(model, name) == value, name
        assert config[_CONFIG_KEYS.get(name, name)] == value, name


@pytest.mark.parametrize("name", ["off", "spohn", "henning"])
def test_binary_round_trip(name, tmp_path):
    """A binary round trip keeps the config, bulk weakening, and a representative evaluation."""
    model = partial_melt.make_partial_melt(name, {"solidus_k": 1550.0, "liquidus_k": 1950.0,
                                                  "liquid_viscosity_pas": 0.7, "bulk_melt_weakening": True,
                                                  "liquid_bulk_modulus_pa": 1.5e10})
    path = str(tmp_path / f"{name}.tpyb")
    model.save_binary(path)
    reloaded = partial_melt.make_partial_melt(name)
    reloaded.load_binary(path)
    assert reloaded.get_config_dict() == model.get_config_dict()
    assert reloaded.bulk_melt_weakening is True
    assert reloaded.calc_bulk_modulus_melt(1700.0, 1.3e11, 1.0e9) == model.calc_bulk_modulus_melt(1700.0, 1.3e11, 1.0e9)
    assert _melt(reloaded, 1700.0) == _melt(model, 1700.0)


@pytest.mark.parametrize("cls_name", _CLASS_NAMES)
def test_isinstance_chain(cls_name):
    model = getattr(partial_melt, cls_name)()
    assert isinstance(model, partial_melt.PartialMeltBase)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)
