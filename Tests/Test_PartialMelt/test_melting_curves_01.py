"""Pressure-dependent solidus and liquidus: the Simon-Glatzel law, its high-pressure branch, and what reads them."""
import math

import pytest

from TidalPy.PartialMelt import HenningPartialMelt, OffPartialMelt, SpohnPartialMelt, make_partial_melt

# Monteux et al. (2016) fits to the peridotite and chondritic-mantle melting curves of Fiquet et al. (2010) and
# Andrault et al. (2011): T = T0 (1 + P / a)^(1 / c), with a second branch above 20 GPa.
_MONTEUX = {
    "solidus_k": 1661.2, "solidus_simon_a_pa": 1.336e9, "solidus_simon_c": 7.437,
    "solidus_transition_pressure_pa": 20.0e9, "solidus_high_k": 2081.8, "solidus_high_simon_a_pa": 1.0169e11,
    "solidus_high_simon_c": 1.226,
    "liquidus_k": 1982.1, "liquidus_simon_a_pa": 6.594e9, "liquidus_simon_c": 5.374,
    "liquidus_transition_pressure_pa": 20.0e9, "liquidus_high_k": 2006.8, "liquidus_high_simon_a_pa": 3.465e10,
    "liquidus_high_simon_c": 1.844,
}


def _simon(temperature, a, c, pressure):
    return temperature * (1.0 + pressure / a) ** (1.0 / c)


def test_constant_curves_are_the_default():
    """With no pressure law the solidus and liquidus are the same at every pressure, and the config is unchanged."""
    model = make_partial_melt("henning", {"solidus_k": 1600.0, "liquidus_k": 2000.0})
    for pressure in (0.0, 1.0e9, 1.0e11):
        assert model.calc_solidus(pressure) == 1600.0
        assert model.calc_liquidus(pressure) == 2000.0
    assert model.solidus_curve["simon_a"] == 0.0
    assert not any(key.endswith(("_simon_a_pa", "_simon_c", "_high_k")) for key in model.get_config_dict())


def test_simon_glatzel_law_and_its_high_pressure_branch():
    """Each branch follows T0 (1 + P / a)^(1 / c) in the absolute pressure; tension holds the zero-pressure value."""
    model = make_partial_melt("henning", _MONTEUX)
    for pressure in (0.0, 5.0e9, 19.9e9):
        assert model.calc_solidus(pressure) == pytest.approx(_simon(1661.2, 1.336e9, 7.437, pressure), rel=1e-14)
        assert model.calc_liquidus(pressure) == pytest.approx(_simon(1982.1, 6.594e9, 5.374, pressure), rel=1e-14)
    for pressure in (20.1e9, 60.0e9, 135.0e9):
        assert model.calc_solidus(pressure) == pytest.approx(_simon(2081.8, 1.0169e11, 1.226, pressure), rel=1e-14)
        assert model.calc_liquidus(pressure) == pytest.approx(_simon(2006.8, 3.465e10, 1.844, pressure), rel=1e-14)
    assert model.calc_solidus(-1.0e9) == 1661.2
    assert math.isnan(model.calc_solidus(float("nan")))


def test_monteux_branches_meet_and_reach_the_published_cmb_values():
    """The two branches of each fit join at 20 GPa, and at the core-mantle boundary (135 GPa) give about 4150 K for
    the solidus (Andrault et al. 2011) and about 4750 K for the liquidus."""
    model = make_partial_melt("henning", _MONTEUX)
    for curve in (model.calc_solidus, model.calc_liquidus):
        assert curve(20.0e9 * (1.0 - 1e-9)) == pytest.approx(curve(20.0e9 * (1.0 + 1e-9)), rel=2e-3)
    assert model.calc_solidus(135.0e9) == pytest.approx(4150.0, rel=0.01)
    assert model.calc_liquidus(135.0e9) == pytest.approx(4750.0, rel=0.01)
    assert model.calc_liquidus(135.0e9) > model.calc_solidus(135.0e9)


def test_melt_fraction_reads_the_curves_at_the_local_pressure():
    """phi = (T - T_sol(P)) / (T_liq(P) - T_sol(P)): the same temperature melts at low pressure, not at high."""
    model = make_partial_melt("henning", _MONTEUX)
    pressure = 5.0e9
    solidus, liquidus = model.calc_solidus(pressure), model.calc_liquidus(pressure)
    temperature = 0.5 * (solidus + liquidus)
    assert model.calc_melt_fraction(temperature, pressure) == pytest.approx(0.5, rel=1e-14)
    assert model.calc_melt_fraction(temperature, 60.0e9) == 0.0
    # The default pressure is zero.
    assert model.calc_melt_fraction(1800.0) == model.calc_melt_fraction(1800.0, 0.0)


@pytest.mark.parametrize("model_name", ["spohn", "henning"])
def test_weakening_laws_anchor_at_the_local_solidus(model_name):
    """The weakening laws use the solidus at the local pressure: just above it the model weakens, just below it it
    returns the pre-melt pair, at any pressure."""
    model = make_partial_melt(model_name, _MONTEUX)
    for pressure in (0.0, 30.0e9):
        solidus = model.calc_solidus(pressure)
        _, viscosity_below, shear_below = model.calc_partial_melt(solidus - 1.0, 1.0e21, 6.0e10, pressure)
        assert (viscosity_below, shear_below) == (1.0e21, 6.0e10)
        melt, viscosity_above, _ = model.calc_partial_melt(solidus + 50.0, 1.0e21, 6.0e10, pressure)
        assert melt > 0.0
        assert viscosity_above < 1.0e21


def test_pressure_reaches_the_bulk_viscosity_melt_effect():
    """The compaction bulk viscosity uses the melt fraction at the local pressure."""
    config = dict(_MONTEUX, bulk_viscosity_melt_weakening=True)
    model = make_partial_melt("henning", config)
    temperature = 0.5 * (model.calc_solidus(0.0) + model.calc_liquidus(0.0))
    assert model.calc_bulk_viscosity_melt(temperature, 1.0e22, 1.0e18) < 1.0e22
    assert model.calc_bulk_viscosity_melt(temperature, 1.0e22, 1.0e18, 60.0e9) == 1.0e22


@pytest.mark.parametrize("model_class", [OffPartialMelt, SpohnPartialMelt, HenningPartialMelt])
def test_constructors_take_the_curve_parameters(model_class):
    """Each model's constructor takes the twelve curve parameters (SI units, no suffix)."""
    model = model_class(
        solidus=1661.2, liquidus=1982.1, solidus_simon_a=1.336e9, solidus_simon_c=7.437,
        solidus_transition_pressure=20.0e9, solidus_high=2081.8, solidus_high_simon_a=1.0169e11,
        solidus_high_simon_c=1.226, liquidus_simon_a=6.594e9, liquidus_simon_c=5.374,
        liquidus_transition_pressure=20.0e9, liquidus_high=2006.8, liquidus_high_simon_a=3.465e10,
        liquidus_high_simon_c=1.844)
    reference = make_partial_melt("henning", _MONTEUX)
    for pressure in (0.0, 10.0e9, 100.0e9):
        assert model.calc_solidus(pressure) == reference.calc_solidus(pressure)
        assert model.calc_liquidus(pressure) == reference.calc_liquidus(pressure)


def test_curves_round_trip_through_the_config():
    """get_config_dict writes every curve key of a pressure-dependent curve, and the factory rebuilds the model."""
    model = make_partial_melt("henning", _MONTEUX)
    config = model.get_config_dict()
    for key, value in _MONTEUX.items():
        assert config[key] == value, key
    rebuilt = make_partial_melt(config.pop("model"), config)
    assert rebuilt.get_config_dict() == model.get_config_dict()
    assert rebuilt.solidus_curve == model.solidus_curve
    assert rebuilt.liquidus_curve == model.liquidus_curve
