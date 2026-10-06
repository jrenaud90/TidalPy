"""System.evolve: a world's spin tracking its spin-orbit equilibria (TidalPy.Structures.system.evolution).

The root search, the combination of a root's two sides, and the windows are checked on synthetic spin balances; the
driver on earth_thermal about the Sun at 0.15 au, e = 0.2, which captures into the 5:2 equilibrium about 1e-6 below
k/2 (a free LSODA integration at a tolerance of 1e-4 on the spin steps over it).
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G, au, year
from TidalPy.Structures import System
from TidalPy.Structures.configs import build_world
from TidalPy.Structures.system import EvolutionResult
from TidalPy.Structures.system.evolution import (
    SpinProbe, build_window, combine_sides, join_segments, locate_root, refine_root, search_points)

MYR = 1.0e6 * year   # [s]
# The trapezoid rule: np.trapezoid from NumPy 2.0, np.trapz before it.
trapezoid = getattr(np, "trapezoid", None) or np.trapz


class SyntheticBalance:
    """A stand-in for the pair rates: the balance is `balance(s)`, the rates [s, balance], the heating s."""

    orbital_frequency = 1.0

    def __init__(self, balance):
        self.balance = balance
        self.num_probes = 0

    def finest_offset(self, finest):
        return finest

    def probe(self, spin_ratio):
        self.num_probes += 1
        value = self.balance(spin_ratio)
        return SpinProbe(spin_ratio, value, np.array([spin_ratio, value]), spin_ratio)


# ======================================================================================================================
# Roots of synthetic balances
# ======================================================================================================================
def test_search_points_meet_the_commensurability_from_the_start_side():
    points = search_points(2.51, 2.5, False, 1.0e-6)
    assert np.all(np.diff(points) < 0.0)
    assert points[0] < 2.51 and points[-1] == pytest.approx(2.25)
    assert np.any(points > 2.5) and np.any(points < 2.5)
    assert np.min(np.abs(points - 2.5)) == pytest.approx(1.0e-6)


def test_locate_finds_a_smooth_root_ahead_of_a_despinning_spin():
    root = 2.5 - 3.0e-4
    rates = SyntheticBalance(lambda s: -1.0e3 * (s - root))
    lower, upper = locate_root(rates, 2.505, 1.0e-12, 1.0e-10)
    spin_ratio, values, _ = combine_sides(lower, upper)
    assert spin_ratio == pytest.approx(root, abs=1.0e-12)
    assert values[1] == pytest.approx(0.0, abs=1.0e-6)


def test_jump_root_is_the_equal_weight_combination_of_its_sides():
    rates = SyntheticBalance(lambda s: 1.0 if s < 2.0 else -1.0)
    lower, upper = locate_root(rates, 2.005, 1.0e-10, 1.0e-10)
    assert upper.spin_ratio - lower.spin_ratio <= 1.0e-10
    spin_ratio, values, _ = combine_sides(lower, upper)
    assert spin_ratio == pytest.approx(2.0, abs=1.0e-10)
    assert values[1] == pytest.approx(0.0, abs=1.0e-12)   # The two sides' balances cancel
    assert rates.num_probes < 200


def test_no_root_ahead_returns_none():
    rates = SyntheticBalance(lambda s: -1.0)
    assert locate_root(rates, 1.505, 1.0e-10, 1.0e-10) is None


def test_refine_converges_on_a_near_discontinuity():
    rates = SyntheticBalance(lambda s: math.tanh(-(s - 1.5) / 1.0e-9))
    lower, upper = refine_root(rates, rates.probe(1.4), rates.probe(1.6), 1.0e-12)
    assert upper.spin_ratio - lower.spin_ratio <= 1.0e-12
    assert combine_sides(lower, upper)[0] == pytest.approx(1.5, abs=1.0e-11)
    assert rates.num_probes < 200


def test_window_stops_short_of_the_unstable_neighbor():
    stable, unstable = 2.5 - 1.0e-6, 2.5 - 4.0e-6
    rates = SyntheticBalance(lambda s: -1.0e12 * (s - stable) * (s - unstable))
    root = combine_sides(*locate_root(rates, 2.501, 1.0e-12, 1.0e-10))[0]
    assert root == pytest.approx(stable, abs=1.0e-11)
    window = build_window(rates, root, 1.0e-8)
    assert unstable < window.lower < root < window.upper
    assert rates.balance(window.lower) > 0.0 > rates.balance(window.upper)


def test_window_steps_past_the_balance_noise():
    """The noise reaches 1e-7 from the root, ten times the resolution: the edges lie past it."""
    root = 2.0 - 5.0e-4
    rates = SyntheticBalance(lambda s: -1.0e4 * (s - root) + 1.0e-3 * math.sin(1.0e13 * s))
    window = build_window(rates, root, 1.0e-8)
    assert window is not None
    assert window.lower < root - 1.0e-7 and window.upper > root + 1.0e-7


def test_window_needs_a_restoring_side():
    rates = SyntheticBalance(lambda s: -1.0)   # A root with no restoring side below it is no equilibrium
    assert build_window(rates, 2.0, 1.0e-8) is None


def test_join_segments_drops_each_repeated_start():
    times = [np.array([0.0, 1.0, 2.0]), np.array([2.0, 3.0]), np.array([3.0, 4.0, 5.0])]
    np.testing.assert_array_equal(join_segments(times, 1), [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    states = [np.array([[0.0, 1.0], [10.0, 11.0]]), np.array([[1.0, 2.0], [11.0, 12.0]])]
    np.testing.assert_array_equal(join_segments(states, 2), [[0.0, 1.0, 2.0], [10.0, 11.0, 12.0]])
    assert join_segments([], 3).shape == (3, 0)


# ======================================================================================================================
# The driver
# ======================================================================================================================
def build_pair(spin_ratio):
    """earth_thermal about the Sun at 0.15 au, e = 0.2, its EOS solved with its temperature profile."""
    planet = build_world("earth_thermal", {"name": "Planet", "radial_solver": {"rtol": 1.0e-10, "atol": 1.0e-14}})
    sun = build_world("sol")
    sun.set_spin_frequency(2.0 * np.pi / (25.38 * 86400.0))
    system = System("Exoplanet")
    system.add_world(
        sun,
        is_star=True)
    system.add_world(
        planet,
        tidal_host=sun,
        semi_major_axis=0.15 * au,
        eccentricity=0.2)
    system.set_tidal_host(sun, planet)
    planet.solve_eos(solve_temperature=True, surface_temperature=system.calc_equilibrium_temperature(planet))
    planet.set_spin_frequency(spin_ratio * system.calc_orbital_frequency(planet))
    return system, planet


@pytest.fixture(scope="module")
def captured():
    system, planet = build_pair(2.53)
    return system, planet, system.evolve(planet, (0.0, 2.0 * MYR), evolve_thermal=False)


def test_a_despinning_world_is_captured_below_5_2(captured):
    _, _, result = captured
    assert isinstance(result, EvolutionResult)
    assert result.success, result.message
    assert result.time[0] == 0.0 and result.time[-1] == pytest.approx(2.0 * MYR)
    assert np.all(np.diff(result.time) > 0.0)
    assert [segment["mode"] for segment in result.segments][-1] == "tracked"
    assert result.tracked[-1]
    # The cold-mantle equilibrium sits about 1e-6 below 5:2.
    assert 2.5 - 1.0e-5 < result.spin_ratio[-1] < 2.5
    assert np.all(np.isfinite(result.spin_ratio)) and np.all(result.tidal_heating > 0.0)
    assert result.temperature is None


def test_a_held_spin_stays_at_its_equilibrium_while_the_orbit_drains(captured):
    _, _, result = captured
    assert np.all(np.abs(result.spin_ratio[result.tracked] - 2.5) < 1.0e-5)
    assert result.semi_major_axis[-1] < result.semi_major_axis[0]
    assert result.eccentricity[-1] < result.eccentricity[0]


def test_the_system_is_left_at_the_final_state(captured):
    system, planet, result = captured
    sun = system.get_tidal_host(planet)
    assert system.get_semi_major_axis(planet) == pytest.approx(result.semi_major_axis[-1], rel=1e-14)
    assert system.get_eccentricity(planet) == pytest.approx(result.eccentricity[-1], rel=1e-14)
    assert sun.spin_frequency == pytest.approx(result.host_spin_frequency[-1], rel=1e-14)
    assert planet.spin_frequency / system.calc_orbital_frequency(planet) == pytest.approx(result.spin_ratio[-1],
                                                                                         rel=1e-12)


def test_orbital_energy_lost_matches_the_heat_while_tracked(captured):
    """At a held spin the orbit pays for the planet's tidal heat (the Sun's own dissipation is far smaller)."""
    system, planet, result = captured
    tracked = result.tracked
    a = result.semi_major_axis[tracked]
    lost = 0.5 * G * system.get_tidal_host(planet).mass * planet.mass * (1.0 / a[-1] - 1.0 / a[0])
    heat = trapezoid(result.tidal_heating[tracked], result.time[tracked])
    assert lost == pytest.approx(heat, rel=0.05)


def test_a_spin_away_from_the_bands_stays_free():
    system, planet = build_pair(3.25)
    result = system.evolve(planet, (0.0, 1.0e-3 * MYR), evolve_thermal=False)
    assert result.success, result.message
    assert not np.any(result.tracked)
    assert result.spin_ratio[-1] < 3.25   # Despinning
    assert [segment["mode"] for segment in result.segments] == ["free"]


def test_the_thermal_state_evolves_with_its_eos():
    system, planet = build_pair(3.25)
    temperature0 = np.array([layer.temperature for layer in planet])
    result = system.evolve(planet, (0.0, 0.01 * MYR))
    assert result.success, result.message
    assert result.temperature.shape == (len(planet), result.time.size)
    np.testing.assert_allclose(result.temperature[:, 0], temperature0)
    assert result.temperature[-1, -1] > temperature0[-1]   # The tidally heated mantle warms
    assert result.num_eos_solves > 0
    np.testing.assert_allclose([layer.temperature for layer in planet], result.temperature[:, -1], rtol=1e-14)


def test_a_wall_time_cap_returns_what_was_integrated():
    system, planet = build_pair(3.25)
    result = system.evolve(planet, (0.0, 1.0 * MYR), evolve_thermal=False, max_wall_time=0.0)
    assert not result.success
    assert "max_wall_time" in result.message


@pytest.mark.parametrize("kwargs, message", [
    (dict(time_span=(1.0, 1.0)), "time span"),
    (dict(time_span=(0.0, 1.0), capture_margin=0.5), "capture_margin"),
    (dict(time_span=(0.0, 1.0), resolution=1.0e-2), "root_tolerance < resolution"),
])
def test_bad_arguments_raise(kwargs, message):
    system, planet = build_pair(2.6)
    with pytest.raises(ValueError, match=message):
        system.evolve(planet, **kwargs)


def test_a_world_without_a_host_raises():
    system = System("solo")
    planet = build_world("earth_thermal")
    system.add_world(planet)
    planet.set_spin_frequency(1.0e-5)
    with pytest.raises(ValueError, match="tidal host"):
        system.evolve(planet, (0.0, 1.0))


def test_a_thermal_run_needs_layers():
    system, planet = build_pair(2.6)
    star = build_world("sol", {"name": "Companion"})   # A world with no layers
    system.add_world(
        star,
        tidal_host=planet,
        semi_major_axis=1.0e9,
        eccentricity=0.1)
    star.set_spin_frequency(1.0e-5)
    with pytest.raises(ValueError, match="layers"):
        system.evolve(star, (0.0, 1.0))
