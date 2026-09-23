"""
Tests for TidalPy.viscosity_x.viscosity — viscosity model hierarchy.

Covers the Arrhenius / Reference / Constant models: the viscosity formulas
(cross-checked against the legacy viscosity_models definitions), the name factory
(incl. unknown-name ValueError and rich-subclass return), config dicts, binary
round-trips, and the isinstance chain.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import math
import os
import tempfile

import pytest

from scipy.constants import R as _R   # molar gas constant [J/mol/K] (same source as the C++ config)


def _import():
    from TidalPy.viscosity_x import viscosity as _mod
    return _mod


# =====================================================================================================================
# Constant
# =====================================================================================================================
@pytest.mark.parametrize("T,P", [(300.0, 0.0), (1500.0, 1.0e10), (2500.0, 5.0e11)])
def test_constant(T, P):
    m = _import().ConstantViscosity(reference_viscosity=1.0e21)
    assert m.calc_viscosity(T, P) == pytest.approx(1.0e21)


# =====================================================================================================================
# Reference
# =====================================================================================================================
def _reference_expected(T, P, eta_ref, T_ref, Ea, Va):
    return eta_ref * math.exp((Ea / _R) * ((1.0 / T) - (1.0 / T_ref)) + (P * Va) / (_R * T))


@pytest.mark.parametrize("T", [800.0, 1000.0, 1500.0])
def test_reference_pressure_always_raises_viscosity(T):
    """A positive activation volume raises the viscosity at every temperature, including at and above T_ref."""
    m = _import().ReferenceViscosity(reference_viscosity=1.0e22, reference_temperature=1000.0,
                                     molar_activation_energy=3.0e5, molar_activation_volume=5.0e-6)
    ratio = m.calc_viscosity(T, 1.0e10) / m.calc_viscosity(T, 0.0)
    assert ratio == pytest.approx(math.exp(5.0e4 / (_R * T)), rel=1e-12)
    assert ratio > 1.0


@pytest.mark.parametrize("T,P", [(1000.0, 0.0), (1500.0, 2.0e10), (1800.0, 1.0e11)])
def test_reference_formula(T, P):
    eta_ref, T_ref, Ea, Va = 1.0e22, 1000.0, 3.0e5, 1.0e-6
    m = _import().ReferenceViscosity(
        reference_viscosity=eta_ref, reference_temperature=T_ref,
        molar_activation_energy=Ea, molar_activation_volume=Va)
    assert m.calc_viscosity(T, P) == pytest.approx(_reference_expected(T, P, eta_ref, T_ref, Ea, Va), rel=1e-9)


def test_reference_at_reference_temperature():
    m = _import().ReferenceViscosity(reference_viscosity=5.0e21, reference_temperature=1200.0)
    assert m.calc_viscosity(1200.0, 0.0) == pytest.approx(5.0e21)


def test_reference_decreases_with_temperature():
    m = _import().ReferenceViscosity(reference_viscosity=1.0e22, reference_temperature=1000.0,
                                     molar_activation_energy=3.0e5)
    assert m.calc_viscosity(1500.0, 0.0) < m.calc_viscosity(1100.0, 0.0)


# =====================================================================================================================
# Arrhenius
# =====================================================================================================================
def _arrhenius_expected(T, P, A, sigma, n, d, mexp, Ea, Va, extra_T):
    eta = A * sigma ** (1.0 - n) * d ** mexp * math.exp((Ea + P * Va) / (_R * T))
    if extra_T:
        eta *= T
    return eta


@pytest.mark.parametrize("T,P", [(250.0, 0.0), (273.0, 1.0e9), (300.0, 5.0e9)])
def test_arrhenius_formula(T, P):
    A, sigma, n, d, mexp, Ea, Va = 1.1e7, 1.0, 1.0, 5.0e-4, 2.0, 59.4e3, 0.0
    m = _import().ArrheniusViscosity(
        arrhenius_coeff=A, stress=sigma, stress_expo=n, grain_size=d,
        grain_size_expo=mexp, molar_activation_energy=Ea, molar_activation_volume=Va,
        additional_temp_dependence=True)
    expected = _arrhenius_expected(T, P, A, sigma, n, d, mexp, Ea, Va, True)
    assert m.calc_viscosity(T, P) == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize("T,P", [(1400.0, 1.0e9), (1600.0, 2.5e10), (2500.0, 1.3e11)])
def test_arrhenius_activation_volume(T, P):
    """With a nonzero activation volume the pressure enters the exponent, and the viscosity rises with depth."""
    A, sigma, n, d, mexp, Ea, Va = 2.6e9, 1.0, 1.0, 1.0e-3, 0.0, 3.0e5, 6.0e-6
    m = _import().ArrheniusViscosity(
        arrhenius_coeff=A, stress=sigma, stress_expo=n, grain_size=d,
        grain_size_expo=mexp, molar_activation_energy=Ea, molar_activation_volume=Va,
        additional_temp_dependence=False)
    expected = _arrhenius_expected(T, P, A, sigma, n, d, mexp, Ea, Va, False)
    assert m.calc_viscosity(T, P) == pytest.approx(expected, rel=1e-9)
    # The pressure term is exp(P Va / R T), which is not a small correction at these pressures.
    assert m.calc_viscosity(T, P) == pytest.approx(
        m.calc_viscosity(T, 0.0) * math.exp(P * Va / (_R * T)), rel=1e-9)
    assert m.calc_viscosity(T, P) > 1.5 * m.calc_viscosity(T, 0.0)


def test_arrhenius_no_extra_temp():
    A, Ea = 2.0, 3.0e5
    m = _import().ArrheniusViscosity(arrhenius_coeff=A, molar_activation_energy=Ea,
                                     additional_temp_dependence=False)
    T, P = 1500.0, 0.0
    expected = A * math.exp(Ea / (_R * T))
    assert m.calc_viscosity(T, P) == pytest.approx(expected, rel=1e-9)


def test_arrhenius_decreases_with_temperature():
    m = _import().ArrheniusViscosity(arrhenius_coeff=1.0, molar_activation_energy=3.0e5)
    assert m.calc_viscosity(2000.0, 0.0) < m.calc_viscosity(1500.0, 0.0)


# Constructor keyword (and property) -> config key. Dimensional config keys carry their unit.
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
    m = _import().ArrheniusViscosity(**params)
    config = m.get_config_dict()
    for name, value in params.items():
        assert getattr(m, name) == value, name
        assert config[_ARRHENIUS_CONFIG_KEYS[name]] == value, name


# =====================================================================================================================
# Factory
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
    mod = _import()
    m = mod.make_viscosity(name)
    assert isinstance(m, getattr(mod, cls_attr))


def test_factory_unknown_raises():
    with pytest.raises(ValueError):
        _import().make_viscosity("not_a_model")


def test_factory_config_override():
    m = _import().make_viscosity("reference", {"reference_viscosity_pas": 7.0e21, "reference_temperature_k": 1300.0})
    d = m.get_config_dict()
    assert d["reference_viscosity_pas"] == pytest.approx(7.0e21)
    assert d["reference_temperature_k"] == pytest.approx(1300.0)


# =====================================================================================================================
# Config dict
# =====================================================================================================================
def test_config_dict_keys():
    m = _import().ArrheniusViscosity()
    d = m.get_config_dict()
    for key in ("model", "arrhenius_coeff", "stress_pa", "grain_size_m",
                "molar_activation_energy_j_mol", "additional_temp_dependence"):
        assert key in d
    assert d["model"] == "arrhenius"


# =====================================================================================================================
# Binary round-trip
# =====================================================================================================================
@pytest.mark.parametrize("name,config", [
    ("constant", {"reference_viscosity_pas": 3.3e21}),
    ("reference", {"reference_viscosity_pas": 1.2e22, "reference_temperature_k": 1100.0}),
    ("arrhenius", {"arrhenius_coeff": 5.0, "additional_temp_dependence": True,
                   "grain_size_expo": 2.0}),
])
def test_binary_round_trip(name, config):
    mod = _import()
    m = mod.make_viscosity(name, config)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, f"{name}.tpyb")
        m.save_binary(path)
        reloaded = mod.make_viscosity(name)
        reloaded.load_binary(path)
    assert reloaded.get_config_dict() == m.get_config_dict()
    assert reloaded.calc_viscosity(1500.0, 1.0e9) == m.calc_viscosity(1500.0, 1.0e9)


# =====================================================================================================================
# isinstance chain
# =====================================================================================================================
def test_isinstance_chain():
    from TidalPy.Utilities_x.classes_x.classes import PhysicsBase, TidalPyBaseClass
    mod = _import()
    m = mod.ArrheniusViscosity()
    assert isinstance(m, mod.ViscosityBase)
    assert isinstance(m, PhysicsBase)
    assert isinstance(m, TidalPyBaseClass)
