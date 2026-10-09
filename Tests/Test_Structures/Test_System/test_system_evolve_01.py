"""System.evolve: a pair of worlds, each spin tracking its spin-orbit equilibria (evolution_.hpp).

The root search, the combination of a root's two sides, and the windows are checked on synthetic spin balances; the
driver on earth_thermal about the Sun at 0.15 au, e = 0.2, which captures into the 5:2 equilibrium about 1e-6 below
k/2 (a free LSODA integration at a tolerance of 1e-4 on the spin steps over it), and on Pluto and Charon, both spins
held synchronous while the orbit circularizes.
"""
import copy
import math
import pickle

import numpy as np
import pytest

from TidalPy.constants import G, au, year
from TidalPy.Structures import System, build_system
from TidalPy.Structures.configs import build_world
from TidalPy.Structures.system import EvolutionResult, PairedEvolutionResult
from TidalPy.Structures.system.system import (
    build_spin_window, locate_spin_root, refine_spin_root, spin_search_points)

MYR = 1.0e6 * year   # [s]
# The trapezoid rule: np.trapezoid from NumPy 2.0, np.trapz before it.
trapezoid = getattr(np, "trapezoid", None) or np.trapz


class Counted:
    """A synthetic balance that counts its calls."""

    def __init__(self, balance):
        self.balance = balance
        self.calls = 0

    def __call__(self, spin_ratio):
        self.calls += 1
        return self.balance(spin_ratio)


# ======================================================================================================================
# Roots of synthetic balances
# ======================================================================================================================
def test_search_points_meet_the_commensurability_from_the_start_side():
    points = spin_search_points(2.51, 2.5, False, 1.0e-6)
    assert np.all(np.diff(points) < 0.0)
    assert points[0] < 2.51 and points[-1] == pytest.approx(2.25)
    assert np.any(points > 2.5) and np.any(points < 2.5)
    assert np.min(np.abs(points - 2.5)) == pytest.approx(1.0e-6)


def test_locate_finds_a_smooth_root_ahead_of_a_despinning_spin():
    root = 2.5 - 3.0e-4
    found = locate_spin_root(lambda s: -1.0e3 * (s - root), 2.505, 1.0e-12, 1.0e-10)
    assert found["spin_ratio"] == pytest.approx(root, abs=1.0e-12)
    assert found["balance"] == pytest.approx(0.0, abs=1.0e-6)


def test_jump_root_is_the_equal_weight_combination_of_its_sides():
    balance = Counted(lambda s: 1.0 if s < 2.0 else -1.0)
    found = locate_spin_root(balance, 2.005, 1.0e-10, 1.0e-10)
    assert found["upper"] - found["lower"] <= 1.0e-10
    assert found["spin_ratio"] == pytest.approx(2.0, abs=1.0e-10)
    assert found["balance"] == pytest.approx(0.0, abs=1.0e-12)   # The two sides' balances cancel
    assert balance.calls == found["num_probes"] < 200


def test_no_root_ahead_returns_none():
    assert locate_spin_root(lambda s: -1.0, 1.505, 1.0e-10, 1.0e-10) is None


def test_refine_converges_on_a_near_discontinuity():
    found = refine_spin_root(lambda s: math.tanh(-(s - 1.5) / 1.0e-9), 1.4, 1.6, 1.0e-12)
    assert found["upper"] - found["lower"] <= 1.0e-12
    assert found["spin_ratio"] == pytest.approx(1.5, abs=1.0e-11)
    assert found["num_probes"] < 200


def test_window_stops_short_of_the_unstable_neighbor():
    stable, unstable = 2.5 - 1.0e-6, 2.5 - 4.0e-6

    def balance(s):
        return -1.0e12 * (s - stable) * (s - unstable)

    root = locate_spin_root(balance, 2.501, 1.0e-12, 1.0e-10)["spin_ratio"]
    assert root == pytest.approx(stable, abs=1.0e-11)
    window = build_spin_window(balance, root, 1.0e-8)
    assert unstable < window["lower"] < root < window["upper"]
    assert balance(window["lower"]) > 0.0 > balance(window["upper"])


def test_window_steps_past_the_balance_noise():
    """The noise reaches 1e-7 from the root, ten times the resolution: the edges lie past it."""
    root = 2.0 - 5.0e-4
    window = build_spin_window(lambda s: -1.0e4 * (s - root) + 1.0e-3 * math.sin(1.0e13 * s), root, 1.0e-8)
    assert window is not None
    assert window["lower"] < root - 1.0e-7 and window["upper"] > root + 1.0e-7


def test_search_points_need_a_positive_finest_offset():
    with pytest.raises(ValueError, match="finest"):
        spin_search_points(2.6, 2.5, False, 0.0)
    np.testing.assert_allclose(spin_search_points(2.6, 2.5, False, 0.3), [2.25, 2.2])   # Python's truncation


def test_window_needs_a_restoring_side():
    assert build_spin_window(lambda s: -1.0, 2.0, 1.0e-8) is None   # No restoring side below: no equilibrium


def test_a_failing_balance_point_is_skipped():
    root = 1.5 - 2.0e-4

    def balance(s):
        if abs(s - 1.5) < 1.0e-6:
            raise RuntimeError("a failed tide solve")
        return -1.0e3 * (s - root)

    assert locate_spin_root(balance, 1.505, 1.0e-12, 1.0e-10)["spin_ratio"] == pytest.approx(root, abs=1.0e-12)


# ======================================================================================================================
# The driver
# ======================================================================================================================
def build_pair(spin_ratio, star_tides=True):
    """earth_thermal about the Sun at 0.15 au, e = 0.2, its EOS solved with its temperature profile."""
    planet = build_world("earth_thermal", {"name": "Planet", "radial_solver": {"rtol": 1.0e-10, "atol": 1.0e-14}})
    sun = build_world("sol")
    sun.set_spin_frequency(2.0 * np.pi / (25.38 * 86400.0))
    if not star_tides:
        sun.set_tide_model(None)
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
    _, _, pair = captured
    assert isinstance(pair, PairedEvolutionResult)
    assert pair.success, pair.message
    assert pair.world_names == ("Planet", "Sol") and list(pair) == ["Planet", "Sol"]
    result = pair["Planet"]
    assert isinstance(result, EvolutionResult)
    assert result.time[0] == 0.0 and result.time[-1] == pytest.approx(2.0 * MYR)
    assert np.all(np.diff(result.time) > 0.0)
    assert [segment["worlds"]["Planet"]["mode"] for segment in pair.segments][-1] == "tracked"
    assert result.tracked[-1]
    # The cold-mantle equilibrium sits about 1e-6 below 5:2.
    assert 2.5 - 1.0e-5 < result.spin_ratio[-1] < 2.5
    assert np.all(np.isfinite(result.spin_ratio)) and np.all(result.tidal_heating > 0.0)
    assert result.temperature is None and pair["Sol"].temperature is None


def test_a_held_spin_stays_at_its_equilibrium_while_the_orbit_drains(captured):
    _, _, pair = captured
    result = pair["Planet"]
    assert np.all(np.abs(result.spin_ratio[result.tracked] - 2.5) < 1.0e-5)
    assert pair.semi_major_axis[-1] < pair.semi_major_axis[0]
    assert pair.eccentricity[-1] < pair.eccentricity[0]


def test_the_pair_rates_are_the_sum_of_the_worlds(captured):
    _, _, pair = captured
    planet, sun = pair["Planet"], pair["Sol"]
    np.testing.assert_allclose(pair.da_dt, planet.da_dt + sun.da_dt, rtol=1e-14)
    np.testing.assert_allclose(pair.de_dt, planet.de_dt + sun.de_dt, rtol=1e-14)
    assert np.all(pair.dn_dt == pytest.approx(-1.5 * (planet.spin_frequency / planet.spin_ratio)
                                              / pair.semi_major_axis * pair.da_dt, rel=1e-9))


def test_the_system_is_left_at_the_final_state(captured):
    system, planet, pair = captured
    sun = system.get_tidal_host(planet)
    assert system.get_semi_major_axis(planet) == pytest.approx(pair.semi_major_axis[-1], rel=1e-14)
    assert system.get_eccentricity(planet) == pytest.approx(pair.eccentricity[-1], rel=1e-14)
    assert sun.spin_frequency == pytest.approx(pair["Sol"].spin_frequency[-1], rel=1e-12)
    assert planet.spin_frequency / system.calc_orbital_frequency(planet) == pytest.approx(
        pair["Planet"].spin_ratio[-1], rel=1e-12)


def test_orbital_energy_lost_matches_the_heat_while_tracked(captured):
    """At a held spin the orbit pays for the planet's tidal heat (the Sun's own dissipation is far smaller)."""
    system, planet, pair = captured
    result = pair["Planet"]
    tracked = result.tracked
    a = pair.semi_major_axis[tracked]
    lost = 0.5 * G * system.get_tidal_host(planet).mass * planet.mass * (1.0 / a[-1] - 1.0 / a[0])
    heat = trapezoid(result.tidal_heating[tracked], pair.time[tracked])
    assert lost == pytest.approx(heat, rel=0.05)


def test_a_spin_away_from_the_bands_stays_free():
    system, planet = build_pair(3.25)
    pair = system.evolve(planet, (0.0, 1.0e-3 * MYR), evolve_thermal=False)
    assert pair.success, pair.message
    assert not np.any(pair["Planet"].tracked)
    assert pair["Planet"].spin_ratio[-1] < 3.25   # Despinning
    assert [segment["worlds"]["Planet"]["mode"] for segment in pair.segments] == ["free"]


def test_a_rigid_partner_keeps_its_spin():
    system, planet = build_pair(3.25, star_tides=False)
    sun_spin = system.get_tidal_host(planet).spin_frequency
    pair = system.evolve(planet, (0.0, 1.0e-3 * MYR), evolve_thermal=False)
    assert pair.success, pair.message
    sun = pair["Sol"]
    assert sun.rigid and sun.num_tide_solves == 0
    assert np.all(sun.spin_frequency == sun_spin) and np.all(sun.tidal_heating == 0.0)
    np.testing.assert_array_equal(pair.da_dt, pair["Planet"].da_dt)
    assert pair.segments[0]["worlds"]["Sol"]["mode"] == "rigid"


def test_the_thermal_state_evolves_with_its_eos():
    system, planet = build_pair(3.25)
    temperature0 = np.array([layer.temperature for layer in planet])
    pair = system.evolve(planet, (0.0, 0.01 * MYR))
    assert pair.success, pair.message
    result = pair["Planet"]
    assert result.temperature.shape == (len(planet), result.time.size)
    np.testing.assert_allclose(result.temperature[:, 0], temperature0)
    assert result.temperature[-1, -1] > temperature0[-1]   # The tidally heated mantle warms
    assert result.num_eos_solves > 0
    assert pair["Sol"].temperature is None and pair["Sol"].num_eos_solves == 0   # The Sun has no layers
    np.testing.assert_allclose([layer.temperature for layer in planet], result.temperature[:, -1], rtol=1e-14)


@pytest.fixture(scope="module")
def pluto_charon():
    system = build_system("pluto_charon_system")
    for world in system:
        if len(world) > 0:
            world.solve_eos()
    return system, system.evolve("charon", (0.0, 0.3 * MYR), evolve_thermal=False)


def test_both_synchronous_spins_are_tracked_while_the_orbit_circularizes(pluto_charon):
    system, pair = pluto_charon
    assert pair.success, pair.message
    assert pair.world_names == ("charon", "pluto")
    for name in pair:
        result = pair[name]
        assert np.all(result.tracked)
        assert np.all(np.abs(result.spin_ratio - 1.0) < 1.0e-6)
        assert result.num_tide_solves > 0
    assert pair.eccentricity[-1] < 0.1 * pair.eccentricity[0]
    assert pair.semi_major_axis[-1] < pair.semi_major_axis[0]
    for name in pair:
        assert system[name].spin_frequency == pytest.approx(pair[name].spin_frequency[-1], rel=1e-12)


def test_results_are_read_only_and_survive_pickling_and_copies(captured):
    _, _, pair = captured
    with pytest.raises(ValueError):
        pair["Planet"].time[0] = -1.0
    for restored in (pickle.loads(pickle.dumps(pair)), copy.deepcopy(pair)):
        assert restored.world_names == pair.world_names
        assert restored.segments == pair.segments
        assert restored.success == pair.success and restored.message == pair.message
        np.testing.assert_array_equal(restored.da_dt, pair.da_dt)
        for name in pair:
            np.testing.assert_array_equal(restored[name].spin_ratio, pair[name].spin_ratio)
            np.testing.assert_array_equal(restored[name].tracked, pair[name].tracked)
            assert restored[name].num_tide_solves == pair[name].num_tide_solves
    single = pickle.loads(pickle.dumps(pair["Planet"]))
    assert single.world_name == "Planet"
    np.testing.assert_array_equal(single.tidal_heating, pair["Planet"].tidal_heating)


def test_either_member_gives_the_same_evolution(pluto_charon):
    """The pair is treated alike from either world."""
    system = build_system("pluto_charon_system")
    for world in system:
        if len(world) > 0:
            world.solve_eos()
    from_pluto = system.evolve("pluto", (0.0, 0.3 * MYR), evolve_thermal=False)
    _, from_charon = pluto_charon
    assert from_pluto.world_names == ("pluto", "charon")
    assert from_pluto.semi_major_axis[-1] == pytest.approx(from_charon.semi_major_axis[-1], rel=1e-12)
    assert from_pluto.eccentricity[-1] == pytest.approx(from_charon.eccentricity[-1], rel=1e-9)
    for name in ("pluto", "charon"):
        assert from_pluto[name].spin_ratio[-1] == pytest.approx(from_charon[name].spin_ratio[-1], rel=1e-12)


def test_a_cap_after_stored_steps_leaves_the_system_at_the_last_step():
    system, planet = build_pair(2.53)
    pair = system.evolve(planet, (0.0, 2.0 * MYR), max_wall_time=0.5)
    assert not pair.success and "max_wall_time" in pair.message
    assert pair.time.size >= 1
    assert system.get_semi_major_axis(planet) == pytest.approx(pair.semi_major_axis[-1], rel=1e-14)
    assert planet.spin_frequency == pytest.approx(pair["Planet"].spin_frequency[-1], rel=1e-12)
    np.testing.assert_allclose([layer.temperature for layer in planet], pair["Planet"].temperature[:, -1],
                               rtol=1e-14)


def test_a_wall_time_cap_returns_what_was_integrated():
    system, planet = build_pair(3.25)
    pair = system.evolve(planet, (0.0, 1.0 * MYR), evolve_thermal=False, max_wall_time=0.0)
    assert not pair.success
    assert "max_wall_time" in pair.message


@pytest.mark.parametrize("kwargs, message", [
    (dict(time_span=(1.0, 1.0)), "time span"),
    (dict(time_span=(0.0, 1.0), capture_margin=0.5), "capture_margin"),
    (dict(time_span=(0.0, 1.0), resolution=1.0e-2), "root_tolerance < resolution"),
    (dict(time_span=(0.0, 1.0), spin_rtol=0.0), "spin_rtol"),
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


def test_a_thermal_world_needs_layer_temperatures():
    system, planet = build_pair(2.6)
    planet[0].temperature = 0.0
    with pytest.raises(ValueError, match="temperature"):
        system.evolve(planet, (0.0, 1.0))
