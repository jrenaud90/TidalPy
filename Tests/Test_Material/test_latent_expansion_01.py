"""The melting curves' pressure slopes, and the latent expansion of a melting range whose curves follow the pressure.

Inside a melting range the melt fraction phi = (T - T_sol) / (T_liq - T_sol) changes with pressure as well as with
temperature when the curves follow the pressure. An isentrope then reads (c_p + L dphi/dT) dT = (alpha T / rho -
L dphi/dP) dP, and the material reports the second term as a latent expansion,
alpha_L = rho L [(1 - phi) dT_sol/dP + phi dT_liq/dP] / ((T_liq - T_sol) T), so an adiabat runs at
dT/dr = -(alpha + alpha_L) g T / c_p. Without it, a convecting mantle whose adiabat enters a pressure-dependent
melting range is pushed onto the liquidus from both sides and its thermal solve cannot be integrated.
"""
from pathlib import Path

import numpy as np
import pytest

import TidalPy
from TidalPy.Material import Material
from TidalPy.PartialMelt import make_melting_curve
from TidalPy.Structures import build_world

_PRESSURES = np.array([0.5e9, 5.0e9, 15.0e9, 25.0e9, 60.0e9, 120.0e9])   # [Pa]
_FINITE_DIFFERENCE_STEP = 1.0e3   # [Pa]

# Peridotite melting curves (Monteux et al. 2016), with their 20 GPa branch transition.
_SOLIDUS = {"model": "simon_glatzel_2", "temperature_k": 1661.2, "simon_a_pa": 1.336e9, "simon_c": 7.437,
            "transition_pressure_pa": 2.0e10, "high_temperature_k": 2081.8, "high_simon_a_pa": 1.0169e11,
            "high_simon_c": 1.226}
_LIQUIDUS = {"model": "simon_glatzel_2", "temperature_k": 1982.1, "simon_a_pa": 6.594e9, "simon_c": 5.374,
             "transition_pressure_pa": 2.0e10, "high_temperature_k": 78.74, "high_simon_a_pa": 4.054e6,
             "high_simon_c": 2.44}
_LATENT_HEAT = 4.0e5   # [J kg-1]
_HOT_MANTLE_TEMPERATURE = 1900.0   # [K] the bundled Earth's mantle 300 K warmer, partially molten at depth


def _melting_material():
    """A constant-density silicate with a silicate melt and the peridotite melting curves."""
    thermal = {"thermal_conductivity_w_mk": 3.3, "heat_capacity_j_kgk": 1200.0}
    return Material(config={
        "latent_heat_j_kg": _LATENT_HEAT,
        "solid": {**thermal, "eos": {"model": "constant", "reference_density_kg_m3": 3300.0,
                                     "bulk_modulus_pa": 1.3e11, "thermal_expansion_1_k": 3.0e-5}},
        "liquid": {**thermal, "eos": {"model": "constant", "reference_density_kg_m3": 3000.0,
                                      "bulk_modulus_pa": 3.0e10, "thermal_expansion_1_k": 6.0e-5},
                   "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.1}},
        "melting": {"solidus": dict(_SOLIDUS), "liquidus": dict(_LIQUIDUS), "weakening": {"model": "henning"}},
    })


@pytest.mark.parametrize("config", [
    {"model": "simon_glatzel", "temperature_k": 1600.0, "simon_a_pa": 2.0e9, "simon_c": 4.0},
    {"model": "simon_glatzel", "temperature_k": 1600.0, "simon_a_pa": 2.0e9, "simon_c": 4.0,
     "reference_pressure_pa": 1.0e9},
    _SOLIDUS,
    _LIQUIDUS,
    {"model": "interpolate", "pressure_pa": [0.0, 10.0e9, 30.0e9, 100.0e9],
     "temperature_k": [1500.0, 2200.0, 2900.0, 4100.0]},
], ids=["simon_glatzel", "simon_glatzel_reference", "simon_glatzel_2_solidus", "simon_glatzel_2_liquidus",
        "interpolate"])
def test_melting_slope_matches_finite_difference(config):
    """Away from kinks, each curve's slope is its central difference."""
    curve = make_melting_curve(config["model"], {key: value for key, value in config.items() if key != "model"})
    # Points away from the interpolated table's knots and the 20 GPa branch transition.
    pressures = np.array([3.3e9, 7.7e9, 41.0e9, 87.0e9])
    expected = (curve.calc_melting_temperature(pressures + _FINITE_DIFFERENCE_STEP)
                - curve.calc_melting_temperature(pressures - _FINITE_DIFFERENCE_STEP)) / (2.0 * _FINITE_DIFFERENCE_STEP)
    np.testing.assert_allclose(curve.calc_melting_slope(pressures), expected, rtol=1e-6)
    assert isinstance(curve.calc_melting_slope(5.0e9), float)


def test_melting_slope_where_curves_are_flat():
    """A constant curve has no slope; a Simon-Glatzel curve is held flat below its reference pressure, and a table
    beyond its ends; a NaN pressure has no slope."""
    constant = make_melting_curve("constant", {"temperature_k": 1600.0})
    np.testing.assert_array_equal(constant.calc_melting_slope(_PRESSURES), 0.0)
    referenced = make_melting_curve(
        "simon_glatzel", {"temperature_k": 1600.0, "simon_a_pa": 2.0e9, "simon_c": 4.0, "reference_pressure_pa": 1.0e9})
    assert referenced.calc_melting_slope(0.5e9) == 0.0
    table = make_melting_curve("interpolate", {"pressure_pa": [1.0e9, 2.0e9], "temperature_k": [1500.0, 1700.0]})
    assert table.calc_melting_slope(0.5e9) == 0.0
    assert table.calc_melting_slope(3.0e9) == 0.0
    assert table.calc_melting_slope(1.5e9) == pytest.approx(200.0 / 1.0e9, rel=1e-14)
    assert np.isnan(constant.calc_melting_slope(np.nan))


def test_latent_expansion_inside_a_pressure_dependent_range():
    """Inside the range the latent expansion is rho L [(1 - phi) T_sol' + phi T_liq'] / ((T_liq - T_sol) T)."""
    material = _melting_material()
    solidus = make_melting_curve("simon_glatzel_2", {k: v for k, v in _SOLIDUS.items() if k != "model"})
    liquidus = make_melting_curve("simon_glatzel_2", {k: v for k, v in _LIQUIDUS.items() if k != "model"})
    for melt_fraction in (0.1, 0.5, 0.9):
        temperature = (solidus.calc_melting_temperature(_PRESSURES)
                       + melt_fraction * (liquidus.calc_melting_temperature(_PRESSURES)
                                          - solidus.calc_melting_temperature(_PRESSURES)))
        state = material.calc_state(_PRESSURES, temperature, use_melting=True, use_pressure_melting=True)
        np.testing.assert_allclose(state["melt_fraction"], melt_fraction, rtol=1e-10)
        span = state["liquidus"] - state["solidus"]
        expected = (state["density"] * _LATENT_HEAT
                    * ((1.0 - melt_fraction) * solidus.calc_melting_slope(_PRESSURES)
                       + melt_fraction * liquidus.calc_melting_slope(_PRESSURES))
                    / (span * temperature))
        np.testing.assert_allclose(state["latent_expansion"], expected, rtol=1e-12)
        assert np.all(state["latent_expansion"] > 0.0)


def test_no_latent_expansion_outside_a_pressure_dependent_range():
    """Zero below the solidus, above the liquidus, with melting curves at zero pressure, and with melting off."""
    material = _melting_material()
    pressure = 25.0e9
    solidus, liquidus = material.calc_melting_range(pressure, use_pressure_melting=True)
    for temperature in (solidus - 50.0, liquidus + 50.0):
        state = material.calc_state(pressure, temperature, use_melting=True, use_pressure_melting=True)
        assert state["latent_expansion"] == 0.0
    inside = 0.5 * (solidus + liquidus)
    assert material.calc_state(pressure, inside, use_melting=True)["latent_expansion"] == 0.0
    assert material.calc_state(pressure, inside, use_pressure_melting=True)["latent_expansion"] == 0.0


def test_convecting_mantle_through_a_pressure_dependent_melting_range():
    """The bundled Earth's convecting mantle, 300 K warmer than bundled, is partially molten along much of its
    adiabat; the thermal solve converges, and inside the melting range the solved adiabat runs at
    dT/dr = -(alpha + alpha_L) g T / c_p."""
    world = build_world(str(Path(TidalPy.__file__).parent / "WorldPack" / "earth_simple.toml"))
    assert world.mantle.use_melting and world.mantle.use_pressure_melting
    world.mantle.temperature = _HOT_MANTLE_TEMPERATURE
    result = world.solve_eos(solve_temperature=True)
    assert result["success"], result["message"]
    assert result["thermal_converged"]

    mantle_i = [layer.name for layer in world].index("mantle")
    # Clear of the boundary layers (millimeters thick in a magma ocean) and of the finite difference's reach.
    margin = max(2.0 * result["layer_boundary_thickness"][mantle_i], 1.0e3)
    radii = np.linspace(world.mantle.radius_inner + margin, world.mantle.radius_outer - margin, 400)
    temperature = world.get_temperature(radii)
    state = world.mantle.calc_state(world.get_pressure(radii), temperature)
    partial = (state["melt_fraction"] > 0.05) & (state["melt_fraction"] < 0.95)
    assert np.count_nonzero(partial) > 10
    assert np.all(state["latent_expansion"][partial] > 0.0)

    step = 10.0   # [m]
    for radius in radii[partial][::5]:
        slope = (world.get_temperature(radius + step) - world.get_temperature(radius - step)) / (2.0 * step)
        point = world.mantle.calc_state(world.get_pressure(radius), world.get_temperature(radius))
        expected = (-(point["thermal_expansion"] + point["latent_expansion"]) * world.get_gravity(radius)
                    * world.get_temperature(radius) / point["heat_capacity"])
        assert slope == pytest.approx(expected, rel=1e-5)
