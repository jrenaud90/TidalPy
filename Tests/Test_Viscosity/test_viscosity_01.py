"""Viscosity models (Arrhenius, reference, constant): formulas against TidalPy 0.7, factory, config, and I/O."""

import math

import pytest
from scipy.constants import R as _R   # molar gas constant [J/mol/K], the same source as the C++ config

from TidalPy.Utilities.classes.classes import PhysicsBase, TidalPyBaseClass
from TidalPy.Viscosity import viscosity as viscosity_models

_CLASS_NAMES = ["ArrheniusViscosity", "ReferenceViscosity", "ConstantViscosity"]


def _reference_expected(T, P, params):
    return params["reference_viscosity"] * math.exp(
        (params["molar_activation_energy"] / _R) * ((1.0 / T) - (1.0 / params["reference_temperature"]))
        + (P * params["molar_activation_volume"]) / (_R * T))


def _arrhenius_expected(T, P, params):
    eta = (params["arrhenius_coeff"] * params["stress"] ** (1.0 - params["stress_expo"])
           * params["grain_size"] ** params["grain_size_expo"]
           * math.exp((params["molar_activation_energy"] + P * params["molar_activation_volume"]) / (_R * T)))
    if params["additional_temp_dependence"]:
        eta *= T
    return eta


# =====================================================================================================================
# Constant and reference
# =====================================================================================================================
@pytest.mark.parametrize("T,P", [(300.0, 0.0), (1500.0, 1.0e10), (2500.0, 5.0e11)])
def test_constant(T, P):
    model = viscosity_models.ConstantViscosity(reference_viscosity=1.0e21)
    assert model.calc_viscosity(T, P) == pytest.approx(1.0e21)


@pytest.mark.parametrize("T", [800.0, 1000.0, 1500.0])
def test_reference_pressure_always_raises_viscosity(T):
    """A positive activation volume raises the viscosity at every temperature, including at and above T_ref."""
    model = viscosity_models.ReferenceViscosity(
        reference_viscosity=1.0e22,
        reference_temperature=1000.0,
        molar_activation_energy=3.0e5,
        molar_activation_volume=5.0e-6)
    ratio = model.calc_viscosity(T, 1.0e10) / model.calc_viscosity(T, 0.0)
    assert ratio == pytest.approx(math.exp(5.0e4 / (_R * T)), rel=1e-12)
    assert ratio > 1.0


@pytest.mark.parametrize("T,P", [(1000.0, 0.0), (1500.0, 2.0e10), (1800.0, 1.0e11)])
def test_reference_formula(T, P):
    params = dict(reference_viscosity=1.0e22, reference_temperature=1000.0, molar_activation_energy=3.0e5,
                  molar_activation_volume=1.0e-6)
    model = viscosity_models.ReferenceViscosity(**params)
    assert model.calc_viscosity(T, P) == pytest.approx(_reference_expected(T, P, params), rel=1e-9)


def test_reference_at_reference_temperature():
    model = viscosity_models.ReferenceViscosity(reference_viscosity=5.0e21, reference_temperature=1200.0)
    assert model.calc_viscosity(1200.0, 0.0) == pytest.approx(5.0e21)


@pytest.mark.parametrize("factory,cool_T,hot_T", [
    (lambda: viscosity_models.ReferenceViscosity(
        reference_viscosity=1.0e22, reference_temperature=1000.0, molar_activation_energy=3.0e5), 1100.0, 1500.0),
    (lambda: viscosity_models.ArrheniusViscosity(arrhenius_coeff=1.0, molar_activation_energy=3.0e5), 1500.0, 2000.0),
], ids=["reference", "arrhenius"])
def test_viscosity_decreases_with_temperature(factory, cool_T, hot_T):
    model = factory()
    assert model.calc_viscosity(hot_T, 0.0) < model.calc_viscosity(cool_T, 0.0)


# =====================================================================================================================
# Arrhenius
# =====================================================================================================================
_ICE_PARAMS = dict(arrhenius_coeff=1.1e7, stress=1.0, stress_expo=1.0, grain_size=5.0e-4, grain_size_expo=2.0,
                   molar_activation_energy=59.4e3, molar_activation_volume=0.0, additional_temp_dependence=True)
_ROCK_PARAMS = dict(arrhenius_coeff=2.6e9, stress=1.0, stress_expo=1.0, grain_size=1.0e-3, grain_size_expo=0.0,
                    molar_activation_energy=3.0e5, molar_activation_volume=6.0e-6, additional_temp_dependence=False)


@pytest.mark.parametrize("T,P", [(250.0, 0.0), (273.0, 1.0e9), (300.0, 5.0e9)])
def test_arrhenius_formula(T, P):
    model = viscosity_models.ArrheniusViscosity(**_ICE_PARAMS)
    assert model.calc_viscosity(T, P) == pytest.approx(_arrhenius_expected(T, P, _ICE_PARAMS), rel=1e-9)


@pytest.mark.parametrize("T,P", [(1400.0, 1.0e9), (1600.0, 2.5e10), (2500.0, 1.3e11)])
def test_arrhenius_activation_volume(T, P):
    """With a nonzero activation volume the pressure enters the exponent and raises the viscosity with depth."""
    model = viscosity_models.ArrheniusViscosity(**_ROCK_PARAMS)
    assert model.calc_viscosity(T, P) == pytest.approx(_arrhenius_expected(T, P, _ROCK_PARAMS), rel=1e-9)
    # exp(P Va / R T) is not a small correction at these pressures.
    pressure_factor = math.exp(P * _ROCK_PARAMS["molar_activation_volume"] / (_R * T))
    assert model.calc_viscosity(T, P) == pytest.approx(model.calc_viscosity(T, 0.0) * pressure_factor, rel=1e-9)
    assert model.calc_viscosity(T, P) > 1.5 * model.calc_viscosity(T, 0.0)


def test_arrhenius_no_extra_temp():
    model = viscosity_models.ArrheniusViscosity(
        arrhenius_coeff=2.0, molar_activation_energy=3.0e5, additional_temp_dependence=False)
    assert model.calc_viscosity(1500.0, 0.0) == pytest.approx(2.0 * math.exp(3.0e5 / (_R * 1500.0)), rel=1e-9)


# Constructor keyword (and property) to config key; dimensional config keys carry their unit.
_ARRHENIUS_CONFIG_KEYS = {
    "arrhenius_coeff":            "arrhenius_coeff",
    "stress":                     "stress_pa",
    "stress_expo":                "stress_expo",
    "grain_size":                 "grain_size_m",
    "grain_size_expo":            "grain_size_expo",
    "molar_activation_energy":    "molar_activation_energy_j_mol",
    "molar_activation_volume":    "molar_activation_volume_m3_mol",
    "additional_temp_dependence": "additional_temp_dependence",
}


def test_arrhenius_properties_match_constructor():
    """Every Arrhenius parameter is a property named like its constructor keyword, emitted under its config key."""
    params = dict(arrhenius_coeff=1.1e7, stress=2.0e6, stress_expo=3.5, grain_size=5.0e-4, grain_size_expo=2.0,
                  molar_activation_energy=5.4e5, molar_activation_volume=1.5e-5,
                  additional_temp_dependence=True)
    model = viscosity_models.ArrheniusViscosity(**params)
    config = model.get_config_dict()
    for name, value in params.items():
        assert getattr(model, name) == value, name
        assert config[_ARRHENIUS_CONFIG_KEYS[name]] == value, name


# =====================================================================================================================
# Factory, config dict, and binary round trip
# =====================================================================================================================
@pytest.mark.parametrize("name,cls_attr", [
    ("arrhenius", "ArrheniusViscosity"),
    ("arr", "ArrheniusViscosity"),
    ("reference", "ReferenceViscosity"),
    ("REF", "ReferenceViscosity"),
    ("constant", "ConstantViscosity"),
    ("const", "ConstantViscosity"),
])
def test_factory_names(name, cls_attr):
    assert isinstance(viscosity_models.make_viscosity(name), getattr(viscosity_models, cls_attr))


def test_factory_unknown_raises():
    with pytest.raises(ValueError):
        viscosity_models.make_viscosity("not_a_model")


def test_factory_config_override():
    model = viscosity_models.make_viscosity(
        "reference", {"reference_viscosity_pas": 7.0e21, "reference_temperature_k": 1300.0})
    config = model.get_config_dict()
    assert config["reference_viscosity_pas"] == pytest.approx(7.0e21)
    assert config["reference_temperature_k"] == pytest.approx(1300.0)


def test_config_dict_keys():
    config = viscosity_models.ArrheniusViscosity().get_config_dict()
    assert {"model", "arrhenius_coeff", "stress_pa", "grain_size_m", "molar_activation_energy_j_mol",
            "additional_temp_dependence"} <= set(config)
    assert config["model"] == "arrhenius"


@pytest.mark.parametrize("name,config", [
    ("constant", {"reference_viscosity_pas": 3.3e21}),
    ("reference", {"reference_viscosity_pas": 1.2e22, "reference_temperature_k": 1100.0}),
    ("arrhenius", {"arrhenius_coeff": 5.0, "additional_temp_dependence": True,
                   "grain_size_expo": 2.0}),
])
def test_binary_round_trip(name, config, tmp_path):
    model = viscosity_models.make_viscosity(name, config)
    path = str(tmp_path / f"{name}.tpyb")
    model.save_binary(path)
    reloaded = viscosity_models.make_viscosity(name)
    reloaded.load_binary(path)
    assert reloaded.get_config_dict() == model.get_config_dict()
    assert reloaded.calc_viscosity(1500.0, 1.0e9) == model.calc_viscosity(1500.0, 1.0e9)


@pytest.mark.parametrize("cls_name", _CLASS_NAMES)
def test_isinstance_chain(cls_name):
    model = getattr(viscosity_models, cls_name)()
    assert isinstance(model, viscosity_models.ViscosityBase)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)
