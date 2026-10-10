"""System.evolve: a pair of worlds, each spin integrated as its offset from a spin-orbit commensurability
(evolution_.hpp).

The driver runs on earth_thermal about the Sun at 0.15 au, e = 0.2, which captures into the 5:2 equilibrium about 1e-6
below k/2, and on Pluto and Charon, both spins held synchronous while the orbit circularizes. The tidal torque passes
smoothly through each commensurability (each world's continuation frequency), which the integration relies on.
"""
import copy
import pickle
import time

import numpy as np
import pytest

import TidalPy
from TidalPy.constants import G, au, year
from TidalPy.Structures import System, build_system
from TidalPy.Structures.configs import build_world
from TidalPy.Structures.system import EvolutionResult, PairedEvolutionResult
from TidalPy.Tides import make_tide

MYR = 1.0e6 * year   # [s]
# The trapezoid rule: np.trapezoid from NumPy 2.0, np.trapz before it.
trapezoid = getattr(np, "trapezoid", None) or np.trapz


def build_pair(spin_ratio, star_tides=True):
    """earth_thermal about the Sun at 0.15 au, e = 0.2, its EOS solved with its temperature profile."""
    planet = build_world("earth_thermal", {"name": "Planet"})
    sun = build_world("sol")
    sun.set_spin_frequency(2.0 * np.pi / (25.38 * 86400.0))
    if not star_tides:
        sun.set_tide_model(None)
    system = System("Exoplanet")
    system.add_world(sun, is_star=True)
    system.add_world(planet, tidal_host=sun, semi_major_axis=0.15 * au, eccentricity=0.2)
    system.set_tidal_host(sun, planet)
    planet.solve_eos(solve_temperature=True, surface_temperature=system.calc_equilibrium_temperature(planet))
    planet.set_spin_frequency(spin_ratio * system.calc_orbital_frequency(planet))
    return system, planet


def build_pluto_charon():
    system = build_system("pluto_charon_system")
    for world in system:
        if len(world) > 0:
            world.solve_eos()
    return system


# ======================================================================================================================
# The tidal torque through a commensurability
# ======================================================================================================================
@pytest.mark.parametrize("commensurability", [1.0, 2.5])
def test_the_torque_passes_smoothly_through_a_commensurability(commensurability):
    """Well below the continuation frequency the resonant mode's dissipation falls linearly to zero, so its torque is
    odd and linear through k/2 with no dead band or jump (at the tight Love-solve tolerance an evolution runs with)."""
    system, planet = build_pair(commensurability)
    planet.set_solver_defaults(radial_solver={"rtol": 1.0e-10, "atol": 1.0e-10})
    n = system.calc_orbital_frequency(planet)

    def torque(offset):
        planet.set_spin_frequency((commensurability + offset) * n)
        return system.calc_world_evolution(planet)["dspin_dt"]

    at_lock = torque(0.0)
    # The resonant mode's part: odd in the offset and linear while its frequency 2 n |offset| is below 1e-16 rad/s.
    resonant = {offset: torque(offset) - at_lock for offset in (-1.0e-11, -1.0e-12, 1.0e-12, 1.0e-11)}
    assert resonant[1.0e-12] < 0.0 < resonant[-1.0e-12]
    assert resonant[1.0e-12] == pytest.approx(-resonant[-1.0e-12], rel=1e-2)
    assert resonant[1.0e-11] == pytest.approx(10.0 * resonant[1.0e-12], rel=1e-2)
    # Restoring further out, through the resonant mode's peak.
    assert torque(-1.0e-8) > at_lock > torque(1.0e-8)


def test_a_warm_ice_shell_raises_the_continuation_frequency():
    """A warm all-solid ice shell is near-fluid far below its Maxwell rate, so its world continues its tides from where
    omega eta falls to the liquid threshold rather than from minimum_frequency; the torque stays smooth and restoring
    through synchrony."""
    minimum = TidalPy.config["numerical"]["minimum_frequency"]
    frequencies = []
    for temperature in (150.0, 240.0):
        system = build_pluto_charon()
        charon = system["charon"]
        charon["hydrosphere"].temperature = temperature
        charon.solve_eos(solve_temperature=True, surface_temperature=system.calc_equilibrium_temperature(charon))
        charon.set_solver_defaults(radial_solver={"rtol": 1.0e-10, "atol": 1.0e-10})
        frequencies.append(charon.calc_continuation_frequency())
        n = system.calc_orbital_frequency("charon")
        torques = []
        for offset in (-1.0e-10, 1.0e-10):
            charon.set_spin_frequency((1.0 + offset) * n)
            torques.append(system.calc_world_evolution(charon)["dspin_dt"])
        assert torques[0] > 0.0 > torques[1]
    assert frequencies[0] == pytest.approx(minimum)
    assert frequencies[1] > 100.0 * minimum
    # An analytic tide model has no Love solve to keep from the near-fluid band.
    charon.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.1], "fixed_q": [10.0]}))
    assert charon.calc_continuation_frequency() == pytest.approx(minimum)


def test_the_torque_has_no_kink_at_the_continuation_frequency():
    """The continuation is smooth: a near-synchronous torque keeps its slope across the frequency where the resonant
    mode reaches the continuation frequency, even for warm ice, whose Im k2 there is far from linear in frequency (a
    cut at that frequency gave slopes 50 times apart on its two sides, which stalled the implicit integrator)."""
    system = build_pluto_charon()
    charon = system["charon"]
    charon["hydrosphere"].temperature = 268.9
    charon.solve_eos(solve_temperature=True, surface_temperature=system.calc_equilibrium_temperature(charon))
    charon.set_solver_defaults(radial_solver={"rtol": 1.0e-10, "atol": 1.0e-10})
    n = system.calc_orbital_frequency("charon")
    edge = charon.calc_continuation_frequency() / (2.0 * n)   # The offset where 2 (s - 1) n reaches it

    def torque(offset):
        charon.set_spin_frequency((1.0 + offset) * n)
        return system.calc_world_evolution(charon)["dspin_dt"]

    below, at, above = (torque(factor * edge) for factor in (0.8, 1.0, 1.2))
    assert 0.5 < (above - at) / (at - below) < 2.0


# ======================================================================================================================
# The driver
# ======================================================================================================================
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
    # The cold-mantle equilibrium sits about 1e-6 below 5:2 (a root search on the torque gives 2.4999992348).
    assert result.spin_ratio[-1] == pytest.approx(2.4999992348, abs=1.0e-8)
    assert pair.segments[-1]["worlds"]["Planet"]["reference_ratio"] == 2.5
    assert np.all(np.isfinite(result.spin_ratio)) and np.all(result.tidal_heating > 0.0)
    assert result.temperature is None and pair["Sol"].temperature is None
    # Held on the equilibrium with long steps once captured.
    assert pair.time.size < 200 and pair.num_jacobians > 0


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


def test_orbital_energy_lost_matches_the_heat_while_held(captured):
    """At a held spin the orbit pays for the planet's tidal heat (the Sun's own dissipation is far smaller)."""
    system, planet, pair = captured
    held = np.abs(pair["Planet"].spin_ratio - 2.5) < 1.0e-5
    a = pair.semi_major_axis[held]
    lost = 0.5 * G * system.get_tidal_host(planet).mass * planet.mass * (1.0 / a[-1] - 1.0 / a[0])
    heat = trapezoid(pair["Planet"].tidal_heating[held], pair.time[held])
    assert lost == pytest.approx(heat, rel=0.05)


def test_the_result_converges_with_the_tolerances(captured):
    """A run at tolerances ten times tighter ends within the looser run's error."""
    _, _, pair = captured
    system, planet = build_pair(2.53)
    tight = system.evolve(planet, (0.0, 2.0 * MYR), evolve_thermal=False, semi_major_axis_rtol=1.0e-6,
                          eccentricity_rtol=1.0e-6, spin_rtol=1.0e-5)
    assert tight.success, tight.message
    assert tight.semi_major_axis[-1] == pytest.approx(pair.semi_major_axis[-1], rel=1.0e-8)
    assert tight.eccentricity[-1] == pytest.approx(pair.eccentricity[-1], rel=1.0e-5)
    assert tight["Planet"].spin_ratio[-1] == pytest.approx(pair["Planet"].spin_ratio[-1], abs=1.0e-9)


def test_a_run_is_deterministic(captured):
    _, _, pair = captured
    system, planet = build_pair(2.53)
    again = system.evolve(planet, (0.0, 2.0 * MYR), evolve_thermal=False)
    np.testing.assert_array_equal(again.time, pair.time)
    np.testing.assert_array_equal(again["Planet"].spin_ratio, pair["Planet"].spin_ratio)


def test_a_spin_away_from_a_lock_despins():
    system, planet = build_pair(3.2)
    pair = system.evolve(planet, (0.0, 1.0e-3 * MYR), evolve_thermal=False)
    assert pair.success, pair.message
    assert pair["Planet"].spin_ratio[-1] < 3.2
    assert pair.segments[0]["worlds"]["Planet"]["reference_ratio"] == 3.0


def test_an_event_near_the_end_of_the_span_finishes_the_run():
    """A segment that restarts just before the span's end starts with a step that fits in what is left."""
    system, planet = build_pair(4.6)
    pair = system.evolve(planet, (0.0, 0.05 * MYR), evolve_thermal=False)
    assert pair.success, pair.message
    assert pair.time[-1] == pytest.approx(0.05 * MYR)


def test_a_capped_despin_leaves_the_system_at_its_last_step():
    """Wherever a cap stops a despin (its spin's reference moving at each half-integer), the system is left at the
    spin of the last stored step."""
    system, planet = build_pair(4.6)
    pair = system.evolve(planet, (0.0, 1.0 * MYR), evolve_thermal=False, max_wall_time=2.0)
    assert pair.time.size >= 1
    assert planet.spin_frequency / system.calc_orbital_frequency(planet) == pytest.approx(
        pair["Planet"].spin_ratio[-1], rel=1e-12)


def test_degree_3_references_are_its_commensurabilities():
    """With degree-3 tides a spin's references are the ratios j / m with m up to 3, the spins at which a mode is
    resonant."""
    system, planet = build_pair(4.6)
    planet.set_tide_config(max_degree_l=3)
    pair = system.evolve(planet, (0.0, 0.3 * MYR), evolve_thermal=False)
    assert pair.success, pair.message
    references = {segment["worlds"]["Planet"]["reference_ratio"] for segment in pair.segments}
    assert {14.0 / 3.0, 4.5} <= references
    for reference in references:
        assert min(abs(reference * order - round(reference * order)) for order in (1, 2, 3)) < 1.0e-12


@pytest.mark.parametrize("method", ["BDF", "LSODA"])
def test_every_method_holds_the_same_equilibrium(captured, method):
    _, _, pair = captured
    system, planet = build_pair(2.53)
    other = system.evolve(planet, (0.0, 2.0 * MYR), evolve_thermal=False, method=method)
    assert other.success, other.message
    assert other["Planet"].spin_ratio[-1] == pytest.approx(pair["Planet"].spin_ratio[-1], abs=1.0e-10)
    assert other.eccentricity[-1] == pytest.approx(pair.eccentricity[-1], rel=1.0e-5)


def test_a_rigid_partner_keeps_its_spin():
    system, planet = build_pair(3.25, star_tides=False)
    sun_spin = system.get_tidal_host(planet).spin_frequency
    pair = system.evolve(planet, (0.0, 1.0e-3 * MYR), evolve_thermal=False)
    assert pair.success, pair.message
    sun = pair["Sol"]
    assert sun.rigid and sun.num_tide_solves == 0
    assert np.all(sun.spin_frequency == sun_spin) and np.all(sun.tidal_heating == 0.0)
    np.testing.assert_array_equal(pair.da_dt, pair["Planet"].da_dt)
    assert pair.segments[0]["worlds"]["Sol"] == {"reference_ratio": None, "rigid": True, "crossing_armed": False}


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
    system = build_pluto_charon()
    return system, system.evolve("charon", (0.0, 0.3 * MYR), evolve_thermal=False)


def test_both_synchronous_spins_hold_while_the_orbit_circularizes(pluto_charon):
    system, pair = pluto_charon
    assert pair.success, pair.message
    assert pair.world_names == ("charon", "pluto")
    for name in pair:
        result = pair[name]
        # Pluto starts on synchrony but off its equilibrium, which the eccentricity shifts; it relaxes back.
        assert np.all(np.abs(result.spin_ratio - 1.0) < 1.0e-7)
        assert abs(result.spin_ratio[-1] - 1.0) < 1.0e-12
        assert result.num_tide_solves > 0
    assert pair.eccentricity[-1] < 0.1 * pair.eccentricity[0]
    assert pair.semi_major_axis[-1] < pair.semi_major_axis[0]
    for name in pair:
        assert system[name].spin_frequency == pytest.approx(pair[name].spin_frequency[-1], rel=1e-12)


def test_either_member_gives_the_same_evolution(pluto_charon):
    """The pair is treated alike from either world."""
    from_pluto = build_pluto_charon().evolve("pluto", (0.0, 0.3 * MYR), evolve_thermal=False)
    _, from_charon = pluto_charon
    assert from_pluto.world_names == ("pluto", "charon")
    assert from_pluto.semi_major_axis[-1] == pytest.approx(from_charon.semi_major_axis[-1], rel=1e-12)
    assert from_pluto.eccentricity[-1] == pytest.approx(from_charon.eccentricity[-1], rel=1e-6)
    for name in ("pluto", "charon"):
        assert from_pluto[name].spin_ratio[-1] == pytest.approx(from_charon[name].spin_ratio[-1], rel=1e-12)


def test_results_are_read_only_and_survive_pickling_and_copies(captured):
    _, _, pair = captured
    with pytest.raises(ValueError):
        pair["Planet"].time[0] = -1.0
    for restored in (pickle.loads(pickle.dumps(pair)), copy.deepcopy(pair)):
        assert restored.world_names == pair.world_names
        assert restored.segments == pair.segments
        assert restored.success == pair.success and restored.message == pair.message
        assert restored.num_rhs_calls == pair.num_rhs_calls and restored.num_jacobians == pair.num_jacobians
        np.testing.assert_array_equal(restored.da_dt, pair.da_dt)
        for name in pair:
            np.testing.assert_array_equal(restored[name].spin_ratio, pair[name].spin_ratio)
            assert restored[name].num_tide_solves == pair[name].num_tide_solves
    single = pickle.loads(pickle.dumps(pair["Planet"]))
    assert single.world_name == "Planet"
    np.testing.assert_array_equal(single.tidal_heating, pair["Planet"].tidal_heating)


# ======================================================================================================================
# Independent checks: analytic tides, conservation, and cost
# ======================================================================================================================
def build_io(model, e0, spin_ratio):
    """Io about a rigid Jupiter under an analytic degree-2 tide (k2 = 0.3, Q = 100 or a 600 s time lag)."""
    moon = build_world("io", {"name": "Moon"})
    if model == "fixed_q":
        moon.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]}))
    else:
        moon.set_tide_model(make_tide("fixed_dt", {"fixed_k": [0.3], "fixed_dt_s": [600.0]}))
    host = build_world("jupiter", {"name": "Host"})
    host.set_tide_model(None)
    system = System("Io")
    system.add_world(host)
    system.add_world(moon, tidal_host=host, semi_major_axis=4.217e8, eccentricity=e0)
    moon.set_spin_frequency(spin_ratio * system.calc_orbital_frequency(moon))
    return system, moon, host


def test_constant_q_circularization_matches_the_analytic_decay():
    """A synchronous body under a constant-Q tide circularizes as de/dt = -e / tau, tau = 1 / ((21/2) (k2/Q) (M/m)
    (R/a)^5 n), with tau following a^6.5 as the orbit shrinks (Goldreich and Soter 1966)."""
    system, moon, host = build_io("fixed_q", 0.01, 1.0)
    a0, n0 = system.get_semi_major_axis(moon), system.calc_orbital_frequency(moon)
    tau = 1.0 / (10.5 * (0.3 / 100.0) * (host.mass / moon.mass) * (moon.radius / a0) ** 5 * n0)
    pair = system.evolve(moon, (0.0, 5.0 * tau), evolve_thermal=False)
    assert pair.success, pair.message
    rate = (pair.semi_major_axis / a0) ** -6.5 / tau
    expected = 0.01 * np.exp(-np.concatenate(([0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(pair.time)))))
    np.testing.assert_allclose(pair.eccentricity, expected, rtol=1e-3)
    assert np.all(np.abs(pair["Moon"].spin_ratio - 1.0) < 1e-6)   # Held synchronous


def test_constant_lag_spin_is_pseudo_synchronous():
    """Under a constant time lag the spin settles at Hut's (1981) pseudo-synchronous rate for the current e."""
    system, moon, _ = build_io("fixed_dt", 0.065, 1.2)
    pair = system.evolve(moon, (0.0, 0.2 * MYR), evolve_thermal=False)
    assert pair.success, pair.message
    e2 = pair.eccentricity[-1] ** 2
    hut = (1.0 + 7.5 * e2 + 5.625 * e2**2 + 0.3125 * e2**3) / ((1.0 + 3.0 * e2 + 0.375 * e2**2) * (1.0 - e2) ** 1.5)
    assert pair["Moon"].spin_ratio[-1] == pytest.approx(hut, rel=1e-3)


def test_angular_momentum_and_energy_balance_at_every_step(captured):
    """Total angular momentum is conserved and the orbit and spins pay for the tidal heat, step by step."""
    system, planet, pair = captured
    sun = system.get_tidal_host(planet)
    m1, m2 = planet.mass, sun.mass
    reduced = m1 * m2 / (m1 + m2)
    a, e = pair.semi_major_axis, pair.eccentricity
    root = np.sqrt(G * (m1 + m2) * a * (1.0 - e**2))
    orbit = reduced * root
    dorbit = reduced * G * (m1 + m2) * ((1.0 - e**2) * pair.da_dt - 2.0 * a * e * pair.de_dt) / (2.0 * root)
    dspins = sum(world.get_moment_of_inertia() * pair[world.name].dspin_dt for world in (planet, sun))
    assert np.max(np.abs(dorbit + dspins) / np.abs(dorbit)) < 1e-6
    heating = pair["Planet"].tidal_heating + pair["Sol"].tidal_heating
    dorbit_energy = G * m1 * m2 / (2.0 * a**2) * pair.da_dt
    dspin_energy = sum(world.get_moment_of_inertia() * pair[world.name].spin_frequency * pair[world.name].dspin_dt
                       for world in (planet, sun))
    np.testing.assert_allclose(heating, -(dorbit_energy + dspin_energy), rtol=1e-6)


def test_a_thermal_pair_evolves_for_millions_of_years():
    """Thermal Pluto and Charon over 10 Myr at the present epoch: both spins stay synchronous and the layers' heat
    flows set their temperatures."""
    system = build_pluto_charon()
    temperatures = {name: np.array([layer.temperature for layer in system[name]]) for name in ("pluto", "charon")}
    pair = system.evolve("charon", (4500.0 * MYR, 4510.0 * MYR))
    assert pair.success, pair.message
    for name in ("pluto", "charon"):
        result = pair[name]
        assert result.temperature.shape == (len(temperatures[name]), pair.time.size)
        np.testing.assert_allclose(result.temperature[:, 0], temperatures[name])
        assert np.all(np.isfinite(result.temperature)) and np.all(result.temperature > 0.0)
        assert abs(result.spin_ratio[-1] - 1.0) < 1.0e-9
        assert result.num_eos_solves > 0


def test_a_cold_gyr_runs_in_seconds():
    """Cold Pluto-Charon over 1 Gyr: the performance target is about a minute per Gyr; this case takes a few
    seconds, so the bound only catches a large regression."""
    system = build_pluto_charon()
    start = time.perf_counter()
    pair = system.evolve("charon", (0.0, 1000.0 * MYR), evolve_thermal=False)
    assert pair.success, pair.message
    assert time.perf_counter() - start < 60.0
    assert pair.time.size < 500


# ======================================================================================================================
# Settings, caps, and failures
# ======================================================================================================================
def test_settings_come_from_the_system_file_then_the_configuration():
    config = build_pluto_charon().get_config_dict()
    config["evolution"] = {"max_wall_time": 0.0}

    def built():
        system = build_system(config)
        for world in system:
            if len(world) > 0:
                world.solve_eos()
        return system

    # The system file's cap stops the run at once; an argument wins over it.
    pair = built().evolve("charon", (0.0, 1.0), evolve_thermal=False)
    assert not pair.success and "max_wall_time" in pair.message
    pair = built().evolve("charon", (0.0, 1.0), evolve_thermal=False, max_wall_time=np.inf)
    assert pair.success, pair.message
    assert TidalPy.config["evolution"]["max_wall_time"] == np.inf


def test_a_system_file_rejects_an_unknown_evolution_key():
    config = build_pluto_charon().get_config_dict()
    config["evolution"] = {"spin_rotl": 1.0e-4}
    with pytest.raises(ValueError, match="spin_rtol"):
        build_system(config)


@pytest.mark.parametrize("setting", [{"method": "nonsense"}, {"spin_rtol": -1.0}, {"evolve_thermal": "false"},
                                     {"max_wall_time": -1.0}])
def test_a_system_file_rejects_an_unusable_evolution_value(setting):
    config = build_pluto_charon().get_config_dict()
    config["evolution"] = setting
    with pytest.raises(ValueError, match=next(iter(setting))):
        build_system(config)


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
    (dict(time_span=(0.0, 1.0), spin_rtol=0.0), "spin_rtol"),
    (dict(time_span=(0.0, 1.0), method="RK45"), "implicit"),
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
