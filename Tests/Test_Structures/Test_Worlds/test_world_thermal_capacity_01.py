"""The heat a layer's profile stores per kelvin of its temperature (layer_thermal_capacity), and the temperature rate
it divides.

The capacity is the integral of rho c_p S over the layer: S is one through an isothermal layer, T(r) / T along a
convecting interior, and the share a conducting stretch moves with the layer's temperature. Inside a melting range c_p
carries the latent heat, so only the part of the layer that is partially molten stores it.
"""
import math
from pathlib import Path

import numpy as np
import pytest

import TidalPy
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld

_EARTH = str(Path(TidalPy.__file__).parent / "WorldPack" / "earth_simple.toml")
_HOT_MANTLE_TEMPERATURE = 1900.0   # [K] the bundled Earth's mantle 300 K warmer, partially molten through much of it
_INTERIOR_POINTS = 40001


def _thermal_earth(mantle_temperature=None):
    world = build_world(_EARTH)
    if mantle_temperature is not None:
        world.mantle.temperature = mantle_temperature
    result = world.solve_eos(solve_temperature=True, surface_temperature=300.0)
    assert result["success"], result["message"]
    return world, result, [layer.name for layer in world].index("mantle")


def _mantle_capacity(world, result, mantle_i, sensible_only=False):
    """The integral of rho c_p S over the mantle by the trapezoid rule on the world's getters, written out here: S is
    T(r) / T along the convecting interior and, across each boundary layer, the steady-conduction share linear in 1/r
    from the interface (which holds) to the interior's end (which moves by 1 at the top and T_base / T at the base).
    With sensible_only, without the latent heat."""
    mantle = world.mantle
    boundary = result["layer_boundary_thickness"][mantle_i]
    layer_temperature = mantle.temperature
    base_scale = result["layer_base_temperature"][mantle_i] / layer_temperature
    interior_base = mantle.radius_inner + boundary
    interior_top = mantle.radius_outer - boundary

    def share(radii, far_radius, held_radius):
        return (1.0 / radii - 1.0 / far_radius) / (1.0 / held_radius - 1.0 / far_radius)

    stretches = (
        (mantle.radius_inner, interior_base, lambda radii, temperature: base_scale * share(
            radii, mantle.radius_inner, interior_base)),
        (interior_base, interior_top, lambda radii, temperature: temperature / layer_temperature),
        (interior_top, mantle.radius_outer, lambda radii, temperature: share(radii, mantle.radius_outer, interior_top)),
    )
    total = 0.0
    for lower, upper, sensitivity in stretches:
        radii = np.linspace(lower, upper, _INTERIOR_POINTS)
        temperature = world.get_temperature(radii)
        state = mantle.calc_state(world.get_pressure(radii), temperature)
        heat_capacity = state["heat_capacity"] - (state["latent_heat_capacity"] if sensible_only else 0.0)
        integrand = 4.0 * math.pi * radii**2 * world.get_density(radii) * heat_capacity * sensitivity(
            radii, temperature)
        total += np.trapezoid(integrand, radii)
    return total


def test_an_isothermal_layer_stores_its_mass_times_its_heat_capacity():
    """A layer with no cooling model is at one temperature, so S = 1 and the capacity is the integral of rho c_p."""
    core = Layer("core", 0, 0.0, 1.0e6, material="simple_rock", temperature=1200.0)
    world = BaseWorld("rock", 1.0e6, 4.0 / 3.0 * math.pi * 3300.0 * 1.0e18)
    world.add_layer(core)
    result = world.solve_eos()
    assert result["success"], result["message"]
    heat_capacity = world.core.calc_state(0.0, 1200.0)["heat_capacity"]
    assert result["layer_thermal_capacity"][0] == pytest.approx(world.core.mass * heat_capacity, rel=1e-9)


def test_a_convecting_mantle_weights_its_adiabat():
    """Along a convecting interior the profile scales with the temperature at its top, so the capacity there is the
    integral of rho c_p T(r) / T; it exceeds M c_p by the adiabat's warming (Stevenson et al. 1983)."""
    world, result, mantle_i = _thermal_earth()
    capacity = result["layer_thermal_capacity"][mantle_i]
    assert capacity == pytest.approx(_mantle_capacity(world, result, mantle_i), rel=1e-5)
    sensible = world.mantle.mass * world.mantle.calc_state(0.0, world.mantle.temperature)["heat_capacity"]
    assert capacity > 1.2 * sensible


def test_the_latent_heat_counts_only_where_the_profile_melts():
    """With the mantle partially molten through much of its interior, the capacity integrates the effective heat
    capacity (latent heat included) over the solved profile, well above its sensible part."""
    world, result, mantle_i = _thermal_earth(_HOT_MANTLE_TEMPERATURE)
    capacity = result["layer_thermal_capacity"][mantle_i]
    assert capacity == pytest.approx(_mantle_capacity(world, result, mantle_i), rel=1e-5)
    assert capacity > 1.1 * _mantle_capacity(world, result, mantle_i, sensible_only=True)


def test_the_temperature_rate_divides_the_heat_budget_by_the_capacity():
    world, result, _ = _thermal_earth()
    for layer_i in range(len(world)):
        budget = (result["layer_heat_flow_in"][layer_i] - result["layer_heat_flow_out"][layer_i]
                  + result["layer_heating"][layer_i])
        capacity = result["layer_thermal_capacity"][layer_i] + result["layer_latent_capacity"][layer_i]
        assert result["layer_temperature_rate"][layer_i] == pytest.approx(budget / capacity, rel=1e-12, abs=1e-30)
