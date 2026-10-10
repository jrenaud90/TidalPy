"""Equation-of-state laws (TidalPy.Material.laws) against closed forms written here independently.

The generic behavior every spec model shares (config and binary round trips, unknown keys, bounds) is tested in
Tests/Test_Utilities/Test_Classes/test_spec_models_01.py; this file tests the physics.
"""

import math

import numpy as np
import pytest

from TidalPy.Material.laws import make_eos

_K0 = 1.3e11
_K0_PRIME = 4.3
_RHO0 = 3300.0


def _bm_pressure(eta, k0=_K0, kp=_K0_PRIME):
    """3rd-order Birch-Murnaghan P(eta), eta = rho / rho0."""
    return 1.5 * k0 * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0)) * (1.0 + 0.75 * (kp - 4.0) * (eta ** (2.0 / 3.0) - 1.0))


def _vinet_pressure(eta, k0=_K0, kp=_K0_PRIME):
    """Vinet P(eta)."""
    x = eta ** (-1.0 / 3.0)
    return 3.0 * k0 * (1.0 - x) / x ** 2 * math.exp(1.5 * (kp - 1.0) * (1.0 - x))


_PRESSURE_LAWS = {"birch_murnaghan": _bm_pressure, "vinet": _vinet_pressure}
_LAW_CONFIG = {"reference_density_kg_m3": _RHO0, "reference_bulk_modulus_pa": _K0, "bulk_modulus_derivative": _K0_PRIME}


@pytest.mark.parametrize("law", sorted(_PRESSURE_LAWS))
@pytest.mark.parametrize("pressure", [1.0e8, 1.0e10, 6.0e10, 1.35e11, 3.6e11])
def test_pressure_law_inverts_its_own_pressure(law, pressure):
    eos = make_eos(law, _LAW_CONFIG)
    eta = eos.calc_density(pressure) / _RHO0
    assert _PRESSURE_LAWS[law](eta) == pytest.approx(pressure, rel=1e-10)


@pytest.mark.parametrize("law", sorted(_PRESSURE_LAWS))
def test_pressure_law_bulk_modulus_is_eta_dp_deta(law):
    eos = make_eos(law, _LAW_CONFIG)
    pressure = 5.0e10
    eta = eos.calc_density(pressure) / _RHO0
    step = 1e-6 * eta
    slope = (_PRESSURE_LAWS[law](eta + step) - _PRESSURE_LAWS[law](eta - step)) / (2.0 * step)
    assert eos.calc_eos(pressure)["bulk_modulus"] == pytest.approx(eta * slope, rel=1e-6)


@pytest.mark.parametrize("law", sorted(_PRESSURE_LAWS))
def test_thermal_pressure_shifts_the_cold_law(law):
    """Hot density at P is the cold density at P - alpha0 K0 (T - T_ref), and only when thermal is on."""
    config = dict(_LAW_CONFIG, thermal_expansion_1_k=3.0e-5, reference_temperature_k=300.0)
    eos = make_eos(law, config)
    temperature = 1800.0
    thermal_pressure = 3.0e-5 * _K0 * (temperature - 300.0)
    assert eos.calc_density(2.0e10, temperature, thermal=True) == pytest.approx(
        eos.calc_density(2.0e10 - thermal_pressure), rel=1e-12)
    assert eos.calc_density(2.0e10, temperature, thermal=False) == eos.calc_density(2.0e10)


def test_constant_law_thermal_expansion():
    eos = make_eos("constant", {"reference_density_kg_m3": 3000.0, "bulk_modulus_pa": 1.0e11,
                                "thermal_expansion_1_k": 2.0e-5, "reference_temperature_k": 300.0})
    state = eos.calc_eos(1.0e9, 1300.0, thermal=True)
    assert state["density"] == pytest.approx(3000.0 * math.exp(-2.0e-5 * 1000.0), rel=1e-14)
    assert state["bulk_modulus"] == 1.0e11
    assert eos.calc_density(1.0e9, 1300.0) == 3000.0


@pytest.mark.parametrize("law", sorted(_PRESSURE_LAWS) + ["murnaghan"])
@pytest.mark.parametrize("pressure", [1.0e10, 3.0e10, 1.3e11])
def test_pressure_law_expansivity_is_its_density_s_own(law, pressure):
    """The thermal pressure alpha0 K0 (T - T_ref) keeps alpha K_T at alpha0 K0, so the reported expansivity is
    alpha0 K0 / K_T, which is -(1/rho) (d rho / dT)_P of the law's own density."""
    config = dict(_LAW_CONFIG, thermal_expansion_1_k=3.0e-5, reference_temperature_k=300.0)
    eos = make_eos(law, config)
    temperature = 1600.0
    state = eos.calc_eos(pressure, temperature, thermal=True)
    assert state["thermal_expansion"] == pytest.approx(3.0e-5 * _K0 / state["bulk_modulus"], rel=1e-13)
    step = 1.0e-2
    slope = (eos.calc_density(pressure, temperature + step, thermal=True)
             - eos.calc_density(pressure, temperature - step, thermal=True)) / (2.0 * step)
    assert state["thermal_expansion"] == pytest.approx(-slope / state["density"], rel=1e-6)
    # Compression stiffens the law, so the expansivity falls with depth.
    assert eos.calc_eos(2.0 * pressure, temperature, thermal=True)["thermal_expansion"] < state["thermal_expansion"]


@pytest.mark.parametrize("law", sorted(_PRESSURE_LAWS) + ["murnaghan", "constant", "polytrope"])
def test_only_the_modified_polytrope_takes_an_anderson_gruneisen_scaling(law):
    """A law whose density follows the temperature has the expansivity of that density, and a polytrope has no
    reference density to scale from, so they refuse the Anderson-Gruneisen keys."""
    with pytest.raises(ValueError, match="anderson_gruneisen_parameter"):
        make_eos(law, {"anderson_gruneisen_parameter": 5.5})


def test_modified_polytrope_anderson_gruneisen_expansivity():
    """alpha = alpha0 exp[(d0 / k) ((rho0 / rho)^k - 1)], and (rho0 / rho)^d0 for k = 0."""
    config = {"reference_density_kg_m3": 8300.0, "thermal_expansion_1_k": 3.0e-5, "anderson_gruneisen_parameter": 5.5}
    for kappa in (0.0, 1.4):
        eos = make_eos("modified_polytrope", dict(config, anderson_gruneisen_exponent=kappa))
        state = eos.calc_eos(3.0e11)
        ratio = 8300.0 / state["density"]
        expected = 3.0e-5 * (ratio ** 5.5 if kappa == 0.0 else math.exp((5.5 / kappa) * (ratio ** kappa - 1.0)))
        assert state["thermal_expansion"] == pytest.approx(expected, rel=1e-13)


def test_adiabatic_bulk_modulus():
    """K_S = K_T (1 + alpha gamma T), and K_T without a Gruneisen parameter."""
    config = dict(_LAW_CONFIG, thermal_expansion_1_k=3.0e-5, gruneisen_parameter=1.3)
    state = make_eos("birch_murnaghan", config).calc_eos(3.0e10, 2000.0)
    assert state["adiabatic_bulk_modulus"] == pytest.approx(
        state["bulk_modulus"] * (1.0 + state["thermal_expansion"] * 1.3 * 2000.0), rel=1e-14)
    plain = make_eos("birch_murnaghan", _LAW_CONFIG).calc_eos(3.0e10, 2000.0)
    assert plain["adiabatic_bulk_modulus"] == plain["bulk_modulus"]


@pytest.mark.parametrize("pressure", [-1.0e9, 0.0, 1.0e9, 3.0e10])
def test_murnaghan_closed_form(pressure):
    eos = make_eos("murnaghan", {"reference_density_kg_m3": 2750.0, "reference_bulk_modulus_pa": 2.0e10,
                                 "bulk_modulus_derivative": 5.0})
    state = eos.calc_eos(pressure)
    if pressure > 0.0:
        assert state["density"] == pytest.approx(2750.0 * (1.0 + 5.0 * pressure / 2.0e10) ** 0.2, rel=1e-14)
        assert state["bulk_modulus"] == pytest.approx(2.0e10 + 5.0 * pressure, rel=1e-14)
    else:
        assert state["density"] == pytest.approx(2750.0 * math.exp(pressure / 2.0e10), rel=1e-14)
        assert state["bulk_modulus"] == 2.0e10


def test_polytrope_closed_form():
    eos = make_eos("polytrope", {"polytropic_constant": 2.1e5, "polytropic_index": 1.0})
    state = eos.calc_eos(1.0e11)
    assert state["density"] == pytest.approx((1.0e11 / 2.1e5) ** 0.5, rel=1e-14)
    assert state["bulk_modulus"] == pytest.approx(2.0e11, rel=1e-14)
    assert eos.calc_eos(0.0)["density"] == 0.0


def test_modified_polytrope_closed_form():
    eos = make_eos("seager", {"reference_density_kg_m3": 8300.0, "polytrope_coefficient": 0.00349,
                              "polytrope_exponent": 0.528})
    pressure = 2.0e11
    state = eos.calc_eos(pressure)
    density = 8300.0 + 0.00349 * pressure ** 0.528
    assert state["density"] == pytest.approx(density, rel=1e-14)
    assert state["bulk_modulus"] == pytest.approx(density / (0.00349 * 0.528 * pressure ** (0.528 - 1.0)), rel=1e-12)


def test_interpolated_law_reads_its_tables():
    eos = make_eos("interpolate", {"radius_m": [0.0, 1.0e6, 2.0e6], "density_kg_m3": [5000.0, 4000.0, 3000.0],
                                   "bulk_modulus_pa": [2.0e11, 1.5e11, 1.0e11]})
    state = eos.calc_eos(0.0, radius=1.5e6)
    assert state["density"] == pytest.approx(3500.0)
    assert state["bulk_modulus"] == pytest.approx(1.25e11)
    # Held at the end values beyond the table.
    assert eos.calc_density(0.0, radius=3.0e6) == pytest.approx(3000.0)
    assert math.isnan(make_eos("interpolate", {"radius_m": [0.0, 1.0], "density_kg_m3": [1.0, 1.0]})
                      .calc_eos(0.0, radius=0.5)["bulk_modulus"])


def test_interpolated_law_rejects_a_bad_table():
    with pytest.raises(ValueError, match="ascending"):
        make_eos("interpolate", {"radius_m": [1.0, 0.0], "density_kg_m3": [1.0, 1.0]})
    with pytest.raises(ValueError, match="does not match"):
        make_eos("interpolate", {"radius_m": [0.0, 1.0], "density_kg_m3": [1.0]})


def test_vectorized_evaluation_broadcasts():
    eos = make_eos("vinet", _LAW_CONFIG)
    pressure = np.linspace(0.0, 1.0e11, 7)
    state = eos.calc_eos(pressure, 1500.0)
    assert state["density"].shape == (7,)
    for i, value in enumerate(pressure):
        assert state["density"][i] == eos.calc_density(float(value), 1500.0)
    grid = eos.calc_density(pressure[:, None], np.array([300.0, 1500.0])[None, :], thermal=True)
    assert grid.shape == (7, 2)
