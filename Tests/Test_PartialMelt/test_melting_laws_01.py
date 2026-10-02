"""Melting curves, melt weakening, and bulk mixing (TidalPy.PartialMelt.melting) against closed forms written here
independently; the generic spec behavior is tested in Tests/Test_Utilities/Test_Classes/test_spec_models_01.py."""

import math

import numpy as np
import pytest

from TidalPy.PartialMelt.melting import (
    make_bulk_modulus_mixing,
    make_bulk_viscosity_mixing,
    make_melt_weakening,
    make_melting_curve,
)

# Monteux et al. (2016) peridotite solidus, the simon_glatzel_2 defaults.
_SOLIDUS = dict(t0=1661.2, a=1.336e9, c=7.437, transition=20.0e9, t0_high=2081.8, a_high=1.0169e11, c_high=1.226)


def _simon_glatzel(pressure, t0, a, c, reference=0.0):
    base = 1.0 + (max(pressure, reference) - reference) / a
    return t0 * base ** (1.0 / c) if base > 0.0 else 0.0


# =====================================================================================================================
# Melting curves
# =====================================================================================================================
@pytest.mark.parametrize("pressure", [-1.0e9, 0.0, 5.0e9, 19.9e9, 20.1e9, 1.3e11])
def test_two_branch_simon_glatzel(pressure):
    curve = make_melting_curve("simon_glatzel_2", {})
    if pressure > _SOLIDUS["transition"]:
        expected = _simon_glatzel(pressure, _SOLIDUS["t0_high"], _SOLIDUS["a_high"], _SOLIDUS["c_high"])
    else:
        expected = _simon_glatzel(pressure, _SOLIDUS["t0"], _SOLIDUS["a"], _SOLIDUS["c"])
    assert curve.calc_melting_temperature(pressure) == pytest.approx(expected, rel=1e-14)


def test_reference_pressure_and_a_falling_curve():
    """A reference pressure shifts the law; a negative a gives a curve that falls with pressure and reaches 0 K where
    its base does (ice Ih's melting curve falls with pressure)."""
    curve = make_melting_curve(
        "simon_glatzel", {"temperature_k": 273.16, "simon_a_pa": -4.0e8, "simon_c": 9.0, "reference_pressure_pa": 611.0})
    assert curve.calc_melting_temperature(0.0) == pytest.approx(273.16)
    assert curve.calc_melting_temperature(2.0e8) == pytest.approx(
        _simon_glatzel(2.0e8, 273.16, -4.0e8, 9.0, 611.0), rel=1e-14)
    assert curve.calc_melting_temperature(2.0e8) < 273.16
    assert curve.calc_melting_temperature(5.0e8) == 0.0


def test_constant_and_interpolated_curves():
    assert make_melting_curve("constant", {"temperature_k": 1500.0}).calc_melting_temperature(1.0e10) == 1500.0
    curve = make_melting_curve("interp", {"pressure_pa": [0.0, 1.0e10], "temperature_k": [1000.0, 2000.0]})
    assert np.allclose(curve.calc_melting_temperature(np.array([0.0, 5.0e9, 2.0e10])), [1000.0, 1500.0, 2000.0])
    assert math.isnan(curve.calc_melting_temperature(math.nan))


def test_zero_simon_a_is_rejected():
    with pytest.raises(ValueError, match="nonzero"):
        make_melting_curve("simon_glatzel", {"simon_a_pa": 0.0})


# =====================================================================================================================
# Melt weakening
# =====================================================================================================================
_SOLID_SHEAR, _SOLID_VISCOSITY = 6.0e10, 1.0e21
_LIQUID_SHEAR, _LIQUID_VISCOSITY = 1.0e-5, 0.2
_T_SOL, _T_LIQ = 1600.0, 2000.0


def _weaken(model, temperature, **parameters):
    return make_melt_weakening(model, parameters).calc_weakening(
        temperature, _T_SOL, _T_LIQ, _SOLID_SHEAR, _SOLID_VISCOSITY, _LIQUID_SHEAR, _LIQUID_VISCOSITY)


@pytest.mark.parametrize("model", ["none", "spohn", "henning"])
def test_every_law_is_solid_below_the_solidus_and_liquid_above_the_liquidus(model):
    assert _weaken(model, 1500.0) == (_SOLID_SHEAR, _SOLID_VISCOSITY)
    assert _weaken(model, 2100.0) == (_LIQUID_SHEAR, _LIQUID_VISCOSITY)


def test_none_keeps_the_solid_until_fully_molten():
    assert _weaken("none", 1999.0) == (_SOLID_SHEAR, _SOLID_VISCOSITY)


_CRIT, _WIDTH = 0.5, 0.05


def _breakdown(phi, shear, viscosity):
    """The framework's pair floored at the liquid's, blended into the liquid across [crit, crit + width] (the
    viscosity log-linearly, the shear modulus linearly), and the liquid's past the band."""
    if phi >= _CRIT + _WIDTH:
        return _LIQUID_SHEAR, _LIQUID_VISCOSITY
    shear, viscosity = max(shear, _LIQUID_SHEAR), max(viscosity, _LIQUID_VISCOSITY)
    if phi > _CRIT:
        blend = (phi - _CRIT) / _WIDTH
        shear = (1.0 - blend) * shear + blend * _LIQUID_SHEAR
        viscosity = math.exp((1.0 - blend) * math.log(viscosity) + blend * math.log(_LIQUID_VISCOSITY))
    return shear, viscosity


@pytest.mark.parametrize("temperature", [1650.0, 1800.0, 1810.0, 1990.0])
def test_spohn(temperature):
    """By default anchored on the solid's own values at the solidus (here the solid's values passed)."""
    phi = (temperature - _T_SOL) / (_T_LIQ - _T_SOL)
    shift = 1.0 / temperature - 1.0 / _T_SOL
    expected = _breakdown(
        phi, _SOLID_SHEAR * 10.0 ** (82000.0 * shift), _SOLID_VISCOSITY * 10.0 ** (27000.0 * shift))
    shear, viscosity = _weaken("spohn", temperature)
    assert viscosity == pytest.approx(expected[1], rel=1e-12)
    assert shear == pytest.approx(expected[0], rel=1e-12)


def test_spohn_with_absolute_anchors_is_the_published_fit():
    """Finite anchors give Fischer and Spohn's 10^(27000 / T - 1) Pa s and 10^(82000 / T - 40.6) Pa at T_sol =
    1600 K, below the breakdown band; the solid at the solidus is then not read."""
    temperature = 1700.0
    shear, viscosity = _weaken(
        "spohn", temperature, fs_visc_log10_at_solidus=15.875, fs_shear_log10_at_solidus=10.65)
    assert viscosity == pytest.approx(10.0 ** (27000.0 / temperature - 1.0), rel=1e-12)
    assert shear == pytest.approx(10.0 ** (82000.0 / temperature - 40.6), rel=1e-12)


def test_spohn_reads_the_solid_at_the_solidus():
    law = make_melt_weakening("spohn", {})
    shear, viscosity = law.calc_weakening(
        1700.0, _T_SOL, _T_LIQ, _SOLID_SHEAR, _SOLID_VISCOSITY, _LIQUID_SHEAR, _LIQUID_VISCOSITY,
        solid_shear_at_solidus=7.0e10, solid_viscosity_at_solidus=3.0e19)
    shift = 1.0 / 1700.0 - 1.0 / _T_SOL
    assert viscosity == pytest.approx(3.0e19 * 10.0 ** (27000.0 * shift), rel=1e-12)
    assert shear == pytest.approx(7.0e10 * 10.0 ** (82000.0 * shift), rel=1e-12)


@pytest.mark.parametrize("temperature", [1650.0, 1790.0, 1805.0, 1815.0, 1900.0])
def test_henning(temperature):
    phi = (temperature - _T_SOL) / (_T_LIQ - _T_SOL)
    break_temperature = _T_SOL + _CRIT * (_T_LIQ - _T_SOL)
    if phi < _CRIT:
        viscosity = _SOLID_VISCOSITY * math.exp(-13.5 * phi)
        shear = _SOLID_SHEAR * math.exp(40000.0 * (1.0 / temperature - 1.0 / _T_SOL))
    else:
        viscosity = _SOLID_VISCOSITY * math.exp(-13.5 * _CRIT) * math.exp(-370.0 * (phi - _CRIT))
        shear = (_SOLID_SHEAR * math.exp(40000.0 * (1.0 / break_temperature - 1.0 / _T_SOL))
                 * math.exp(-700.0 * (phi - _CRIT)))
    expected_shear, expected_viscosity = _breakdown(phi, shear, viscosity)
    result_shear, result_viscosity = _weaken("henning", temperature)
    assert result_viscosity == pytest.approx(expected_viscosity, rel=1e-12)
    assert result_shear == pytest.approx(expected_shear, rel=1e-12)


@pytest.mark.parametrize("model", ["spohn", "henning"])
def test_weakening_is_continuous_through_the_melting_range(model):
    """No step anywhere from below the solidus to above the liquidus: on a 0.01 K grid the viscosity changes by at
    most a few percent per step (a step at the band's end was a factor of 1e10), and the shear modulus by a small
    fraction of the solid's."""
    temperatures = np.arange(_T_SOL - 1.0, _T_LIQ + 1.0, 0.01)
    pairs = np.array([_weaken(model, temperature) for temperature in temperatures])
    assert np.max(np.abs(np.diff(np.log(pairs[:, 1])))) < 0.05
    assert np.max(np.abs(np.diff(pairs[:, 0]))) < 1.0e-2 * _SOLID_SHEAR


def test_a_zero_width_band_is_a_step_into_the_liquid():
    phi_temperature = _T_SOL + _CRIT * (_T_LIQ - _T_SOL)
    law = {"crit_melt_frac_width": 0.0}
    assert _weaken("henning", phi_temperature + 1.0e-6, **law) == (_LIQUID_SHEAR, _LIQUID_VISCOSITY)
    assert _weaken("henning", phi_temperature - 1.0, **law)[1] > 1.0e3 * _LIQUID_VISCOSITY


# =====================================================================================================================
# Bulk mixing
# =====================================================================================================================
@pytest.mark.parametrize("phi", [0.0, 0.1, 0.5, 1.0])
def test_hashin_shtrikman(phi):
    k_s, k_l, mu = 1.3e11, 2.0e10, 6.0e10
    expected = k_s if phi == 0.0 else k_s + phi / (1.0 / (k_l - k_s) + (1.0 - phi) / (k_s + 4.0 / 3.0 * mu))
    assert make_bulk_modulus_mixing("hs", {}).calc_bulk_modulus(k_s, k_l, mu, phi) == pytest.approx(expected, rel=1e-14)


def test_hashin_shtrikman_reaches_the_liquid_when_fully_molten():
    assert make_bulk_modulus_mixing("hs", {}).calc_bulk_modulus(1.3e11, 2.0e10, 0.0, 1.0) == pytest.approx(2.0e10)


@pytest.mark.parametrize("phi, exponent", [(0.0, 1.0), (0.1, 1.0), (0.3, 0.0)])
def test_compaction(phi, exponent):
    law = make_bulk_viscosity_mixing("compaction", {"coefficient": 1.3, "exponent": exponent})
    expected = 1.0e22 if phi == 0.0 else 1.0 / (1.0 / 1.0e22 + phi ** exponent / (1.3 * 1.0e18))
    assert law.calc_bulk_viscosity(1.0e22, 1.0e18, phi) == pytest.approx(expected, rel=1e-14)


def test_compaction_without_a_solid_dashpot():
    law = make_bulk_viscosity_mixing("mckenzie", {})
    assert law.calc_bulk_viscosity(math.nan, 1.0e18, 0.2) == pytest.approx(1.0e18 / 0.2, rel=1e-14)
