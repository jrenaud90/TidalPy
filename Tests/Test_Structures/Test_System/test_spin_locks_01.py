"""Spin-orbit resonance locks of the System evolution calls: the hold test, the Filippov weight and combined rates, the
release, equilibria narrower than the tolerance, the probe leaving the world untouched, pair mode, and the arguments."""
import math

import numpy as np
import pytest

import TidalPy
from TidalPy.configurations import find_invalid_config_values, find_unknown_config_keys, get_packaged_config
from TidalPy.Dynamics import OrbitSolver
from TidalPy.Structures import build_world
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.base import BaseWorld

# An Io-like moon about a Jupiter-mass host.
HOST_MASS = 1.898e27
SEMI_MAJOR_AXIS = 4.2e8
TOLERANCE = 1.0e-5
RADIUS = 1.6e6
DENSITY = 3500.0
LAYER_MATERIAL = {"solid": {
    "eos": {"model": "constant", "reference_density_kg_m3": DENSITY, "bulk_modulus_pa": 2.0e11},
    "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10}}}
# The stable zeros of the Maxwell moon's spin balance below 3:2 (spin ratio minus 1.5), found by root search: a warm one
# at viscosity 1e17 Pa s and e = 0.01 (restoring band down to -5.9e-4), and a cold one at 1e19 Pa s and e = 0.0008 whose
# restoring band, from -3.58e-6 to -3.37e-7, is narrower than TOLERANCE.
WARM_ROOT = -2.0399e-5
COLD_ROOT = -3.3711e-7
COLD_UNSTABLE_ROOT = -3.5828e-6


def cpl_moon(name="moon", eccentricity_truncation=10):
    """A two-layer moon with a constant-phase-lag tide model, whose torque jumps at every half-integer spin ratio."""
    return build_world({
        "name": name,
        "type": "terrestrial",
        "radius_m": 1.6e6,
        "mass_kg": 8.9e22,
        "tides": {
            "global_tidal_model": "cpl",
            "max_degree_l": 2,
            "eccentricity_trunc_lvl": eccentricity_truncation,
            "obliquity_trunc_lvl": "off",
            "fixed_k": [0.3],
            "fixed_q": [50.0],
        },
        "layers": {
            "core": {"radius_fraction": 0.5, "material": LAYER_MATERIAL},
            "mantle": {"radius_fraction": 1.0, "material": LAYER_MATERIAL},
        },
    })


def maxwell_moon(viscosity):
    """A homogeneous Maxwell moon with quasi-homogeneous Love numbers, its EOS solved."""
    material = {"solid": {**LAYER_MATERIAL["solid"],
                          "shear_viscosity": {"model": "constant", "reference_viscosity_pas": viscosity},
                          "shear_rheology": {"model": "maxwell"}}}
    world = build_world({
        "name": "moon",
        "type": "terrestrial",
        "radius_m": RADIUS,
        "mass_kg": (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY,
        "tides": {
            "global_tidal_model": "rheology",
            "love_method": "homogeneous",
            "max_degree_l": 2,
            "eccentricity_trunc_lvl": 10,
            "obliquity_trunc_lvl": "off",
        },
        "eos_solver": {"solve_temperature": False},
        "layers": {"mantle": {"radius_fraction": 1.0, "material": material}},
    })
    world.solve_eos()
    return world


def system_with(world, eccentricity, spin_ratio, host=None):
    """``world`` about a host (a rigid Jupiter-mass point unless given), spinning at ``spin_ratio`` times n."""
    system = System("pair")
    system.add_world(host if host is not None else BaseWorld("host", 7.0e7, HOST_MASS))
    system.add_world(world, tidal_host=0, semi_major_axis=SEMI_MAJOR_AXIS, eccentricity=eccentricity)
    world.set_spin_frequency(spin_ratio * system.calc_orbital_frequency(world))
    return system


def balance(row, dn_dt_other=0.0):
    """The spin balance dspin/dt - s dn/dt of an evolution row [rad s-2]."""
    spin_ratio = row["spin_frequency"] / row["orbital_frequency"]
    return row["dspin_dt"] - spin_ratio * (row["dn_dt"] + dn_dt_other)


def free_row(system, world, spin_ratio):
    """The locks-off row at another spin ratio; the world's own spin is put back afterwards."""
    own_spin = world.spin_frequency
    world.set_spin_frequency(spin_ratio * system.calc_orbital_frequency(world))
    row = system.calc_world_evolution(world, use_locks=False)
    world.set_spin_frequency(own_spin)
    return row


def expected_hold(system, world, held):
    """The two sides of a held single-body row, solved with locks off: (current row, probe row, b0, b1, B, w)."""
    spin_ratio = held["spin_frequency"] / held["orbital_frequency"]
    current = free_row(system, world, spin_ratio)
    probe = free_row(system, world, held["spin_lock_bracket"])
    # Both balances are taken at the held state's own spin ratio.
    balance_current = balance(current)
    balance_probe = probe["dspin_dt"] - spin_ratio * probe["dn_dt"]
    balance_held = balance_current**2 / (balance_current - balance_probe)
    weight = (balance_held - balance_probe) / (balance_current - balance_probe)
    return current, probe, balance_current, balance_probe, balance_held, weight


def assert_held_spin_is_consistent(row, balance_held, rate_scale, dn_dt_other=0.0):
    """The spin rate is the spin model's from the combined dU/dO, and equals s dn/dt + B to round-off."""
    assert row["dspin_dt"] == row["host_mass"] * row["dU_dO"] / row["moment_of_inertia"]
    spin_ratio = row["spin_frequency"] / row["orbital_frequency"]
    # The combination cancels rates of either sign, so round-off is measured against each side's spin rate.
    assert abs(row["dspin_dt"] - (spin_ratio * (row["dn_dt"] + dn_dt_other) + balance_held)) <= 1.0e-12 * rate_scale


# =====================================================================================================================
# Locks off
# =====================================================================================================================
def test_locks_off_is_the_plain_solve():
    world = cpl_moon()
    system = system_with(world, 0.3, 1.5 - 0.3 * TOLERANCE)
    default = system.calc_world_evolution(world)
    off = system.calc_world_evolution(world, use_locks=False)
    assert default == pytest.approx(off, nan_ok=True, rel=0.0, abs=0.0)
    assert not off["spin_locked"]
    assert math.isnan(off["spin_lock_weight"]) and math.isnan(off["spin_lock_bracket"])
    # The rates come from the world's own solve, unchanged.
    assert off["tidal_heating"] == world.get_tidal_heating()
    assert off["dU_dM_minus_dw"] == world.get_tidal_dU_dM_minus_dw()
    rates = OrbitSolver().calc_derivatives(
        off["orbital_frequency"], off["semi_major_axis"], off["eccentricity"], off["target_mass"], off["host_mass"],
        off["dU_dM"], off["dU_dw"], off["dU_dM_minus_dw"])
    assert (off["da_dt"], off["de_dt"], off["dn_dt"]) == (rates["da_dt"], rates["de_dt"], rates["dn_dt"])
    assert off["dspin_dt"] == world.calc_spin_derivative(HOST_MASS)


def test_locks_on_away_from_any_equilibrium_is_unchanged():
    world = cpl_moon()
    system = system_with(world, 0.3, 2.23)
    off = system.calc_world_evolution(world, use_locks=False)
    on = system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)
    assert on == pytest.approx(off, nan_ok=True, rel=0.0, abs=0.0)


# =====================================================================================================================
# Held at a torque jump (constant phase lag)
# =====================================================================================================================
@pytest.mark.parametrize("offset, probe_side", [(-0.3 * TOLERANCE, 1.0), (0.3 * TOLERANCE, -1.0)])
def test_cpl_spin_is_held_at_the_jump(offset, probe_side):
    world = cpl_moon()
    spin_ratio = 1.5 + offset
    system = system_with(world, 0.3, spin_ratio)
    held = system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)
    assert held["spin_locked"]
    assert held["spin_lock_bracket"] == pytest.approx(spin_ratio + probe_side * TOLERANCE, rel=1.0e-15)

    # The held balance, the weight, and every combined value from the two sides, solved with locks off.
    current, probe, balance_current, balance_probe, balance_held, weight = expected_hold(system, world, held)
    # The restoring straddle: positive below, negative above.
    assert probe_side * balance_current > 0.0 > probe_side * balance_probe
    # B lies between the two balances, on the current side of zero.
    assert min(balance_current, balance_probe) < balance_held < max(balance_current, balance_probe)
    assert balance_held * balance_current > 0.0
    assert held["spin_lock_weight"] == pytest.approx(weight, rel=1.0e-12)
    for key in ("tidal_heating", "dU_dM", "dU_dw", "dU_dO", "dU_dM_minus_dw", "da_dt", "de_dt", "dn_dt"):
        assert held[key] == pytest.approx(weight * current[key] + (1.0 - weight) * probe[key], rel=1.0e-12)
    rate_scale = max(abs(current["dspin_dt"]), abs(probe["dspin_dt"]))
    assert_held_spin_is_consistent(held, balance_held, rate_scale)
    # Each side conserves energy, so the combination leaves a residual of the order of the tolerance.
    assert abs(held["energy_residual"]) < 10.0 * TOLERANCE * held["tidal_heating"]


def test_probe_commits_nothing_to_the_world():
    world = cpl_moon()
    system = system_with(world, 0.3, 1.5 - 0.3 * TOLERANCE)
    own = system.calc_world_evolution(world, use_locks=False)
    layer_heating = [world.get_layer_tidal_heating(i) for i in range(2)]
    held = system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)
    assert held["spin_locked"]
    assert held["tidal_heating"] != own["tidal_heating"]
    # The world keeps the solve at its own spin: its result and its layers' heating.
    assert world.get_tidal_heating() == own["tidal_heating"]
    assert world.get_tidal_potential_derivatives()["dU_dO"] == own["dU_dO"]
    assert [world.get_layer_tidal_heating(i) for i in range(2)] == layer_heating
    assert [layer.get_tidal_heating() for layer in world] == layer_heating


def test_lower_eccentricity_releases_the_spin():
    # Below e of about 0.26 the jump at 3:2 no longer outweighs the rest of the constant-phase-lag torque.
    world = cpl_moon()
    system = system_with(world, 0.2, 1.5 - 0.3 * TOLERANCE)
    on = system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)
    off = system.calc_world_evolution(world, use_locks=False)
    assert not on["spin_locked"]
    assert on == pytest.approx(off, nan_ok=True, rel=0.0, abs=0.0)


# =====================================================================================================================
# Smooth equilibria of a Maxwell moon
# =====================================================================================================================
@pytest.fixture(scope="module")
def warm_moon():
    return maxwell_moon(1.0e17)


@pytest.mark.parametrize("offset", [0.5 * TOLERANCE, -0.5 * TOLERANCE, 0.9 * TOLERANCE])
def test_warm_equilibrium_off_the_commensurability_is_held_within_the_tolerance(warm_moon, offset):
    system = system_with(warm_moon, 0.01, 1.5 + WARM_ROOT + offset)
    held = system.calc_world_evolution(warm_moon, use_locks=True, lock_tolerance=TOLERANCE)
    assert held["spin_locked"]
    assert 0.0 < held["spin_lock_weight"] < 1.0
    current, probe, balance_current, _, balance_held, weight = expected_hold(system, warm_moon, held)
    assert held["spin_lock_weight"] == pytest.approx(weight, rel=1.0e-10)
    # The held spin still moves toward the root, more slowly than the free spin.
    assert 0.0 < balance_held / balance_current < 1.0
    rate_scale = max(abs(current["dspin_dt"]), abs(probe["dspin_dt"]))
    assert_held_spin_is_consistent(held, balance_held, rate_scale)


def test_held_spin_rate_is_continuous_at_the_hold_edge_and_vanishes_at_the_root(warm_moon):
    # Just inside the outer edge of the hold, the held spin balance is the free one.
    system = system_with(warm_moon, 0.01, 1.5 + WARM_ROOT + (1.0 - 1.0e-3) * TOLERANCE)
    held = system.calc_world_evolution(warm_moon, use_locks=True, lock_tolerance=TOLERANCE)
    free = system.calc_world_evolution(warm_moon, use_locks=False)
    assert held["spin_locked"]
    assert balance(held) == pytest.approx(balance(free), rel=1.0e-2)
    # Near the root it vanishes faster than the free balance: B / b0 = b0 / (b0 - b1) goes to zero.
    system = system_with(warm_moon, 0.01, 1.5 + WARM_ROOT + 1.0e-3 * TOLERANCE)
    held = system.calc_world_evolution(warm_moon, use_locks=True, lock_tolerance=TOLERANCE)
    free = system.calc_world_evolution(warm_moon, use_locks=False)
    assert held["spin_locked"]
    assert abs(balance(held)) < 1.0e-2 * abs(balance(free))


@pytest.mark.parametrize("offset", [3.0 * TOLERANCE, -3.0 * TOLERANCE])
def test_warm_equilibrium_farther_than_the_tolerance_is_not_held(warm_moon, offset):
    system = system_with(warm_moon, 0.01, 1.5 + WARM_ROOT + offset)
    assert not system.calc_world_evolution(warm_moon, use_locks=True, lock_tolerance=TOLERANCE)["spin_locked"]


def test_warm_equilibrium_vanishes_at_lower_eccentricity(warm_moon):
    # At e = 0.004 the moon has no zero near 3:2, so the spin is free.
    system = system_with(warm_moon, 0.004, 1.5 + WARM_ROOT + 0.5 * TOLERANCE)
    assert not system.calc_world_evolution(warm_moon, use_locks=True, lock_tolerance=TOLERANCE)["spin_locked"]


def test_equilibrium_narrower_than_the_tolerance():
    world = maxwell_moon(1.0e19)
    band = COLD_ROOT - COLD_UNSTABLE_ROOT
    assert band < 0.5 * TOLERANCE
    # Half a tolerance above the root the probe lands below the whole restoring band: free.
    system = system_with(world, 0.0008, 1.5 + COLD_ROOT + 0.5 * TOLERANCE)
    assert not system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)["spin_locked"]
    # A probe that lands inside the band holds it, one tolerance above; a tolerance below the band's width holds it at
    # the root.
    system = system_with(world, 0.0008, 1.5 + COLD_ROOT + TOLERANCE - 0.5 * band)
    assert system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)["spin_locked"]
    system = system_with(world, 0.0008, 1.5 + COLD_ROOT + 0.1 * band)
    assert system.calc_world_evolution(world, use_locks=True, lock_tolerance=0.5 * band)["spin_locked"]


# =====================================================================================================================
# Warnings
# =====================================================================================================================
def test_probes_of_a_synchronous_spin_do_not_warn(spdlog_text):
    world = cpl_moon("synchronous")
    system = system_with(world, 0.1, 1.0)
    held = system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)
    assert held["spin_locked"]
    assert "spins at" not in spdlog_text()
    # The same probe as a solve of the world's own does warn.
    world.set_spin_frequency(held["spin_lock_bracket"] * held["orbital_frequency"])
    system.calc_world_evolution(world)
    assert "world 'synchronous' spins at" in spdlog_text()


# =====================================================================================================================
# Pair mode
# =====================================================================================================================
def test_pair_with_a_rigid_host_holds_as_a_single_body():
    world = cpl_moon()
    system = system_with(world, 0.3, 1.5 - 0.3 * TOLERANCE)
    pair = system.calc_pair_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)
    single = system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)
    assert pair["world"] == pytest.approx(single, nan_ok=True, rel=0.0, abs=0.0)
    assert pair["world"]["spin_locked"] and not pair["host"]["spin_locked"]
    assert pair["da_dt"] == pair["world"]["da_dt"]


def test_pair_holds_both_dissipating_bodies():
    host = cpl_moon("planet")
    world = cpl_moon()
    system = system_with(world, 0.3, 1.5 - 0.3 * TOLERANCE, host=host)
    host.set_spin_frequency((1.5 + 0.3 * TOLERANCE) * system.calc_orbital_frequency(world))
    free = system.calc_pair_evolution(world, use_locks=False)
    held = system.calc_pair_evolution(world, use_locks=True, lock_tolerance=TOLERANCE)
    assert held["world"]["spin_locked"] and held["host"]["spin_locked"]
    # Each body's hold takes the other body's contribution to dn/dt as solved at its own spin.
    for body, other in (("world", "host"), ("host", "world")):
        row = held[body]
        assert row["dspin_dt"] == row["host_mass"] * row["dU_dO"] / row["moment_of_inertia"]
        # The held balance lies between zero and the free one, on the same side.
        balance_free = balance(free[body], free[other]["dn_dt"])
        balance_held = balance(row, free[other]["dn_dt"])
        assert 0.0 < balance_held / balance_free < 1.0
    assert held["dn_dt"] == held["world"]["dn_dt"] + held["host"]["dn_dt"]
    assert held["tidal_heating_total"] == held["world"]["tidal_heating"] + held["host"]["tidal_heating"]
    assert abs(held["energy_residual"]) < 10.0 * TOLERANCE * held["tidal_heating_total"]


def test_system_evolution_passes_the_locks():
    world = cpl_moon()
    system = system_with(world, 0.3, 1.5 - 0.3 * TOLERANCE)
    rows = system.calc_system_evolution(use_locks=True, lock_tolerance=TOLERANCE)
    assert rows[1]["spin_locked"]
    assert rows[1] == pytest.approx(
        system.calc_world_evolution(world, use_locks=True, lock_tolerance=TOLERANCE), nan_ok=True, rel=0.0, abs=0.0)


# =====================================================================================================================
# Arguments and configuration
# =====================================================================================================================
@pytest.mark.parametrize("use_locks", [True, False])
@pytest.mark.parametrize("lock_tolerance", [0.0, -1.0e-5, 0.25, 0.3, math.inf, math.nan])
def test_lock_tolerance_is_validated(use_locks, lock_tolerance):
    world = cpl_moon()
    system = system_with(world, 0.3, 1.2)
    if math.isnan(lock_tolerance) and not use_locks:
        # NaN leaves the tolerance unset, which locks off allow.
        system.calc_world_evolution(world, use_locks=use_locks, lock_tolerance=lock_tolerance)
        return
    for call in (lambda: system.calc_world_evolution(world, use_locks=use_locks, lock_tolerance=lock_tolerance),
                 lambda: system.calc_pair_evolution(world, use_locks=use_locks, lock_tolerance=lock_tolerance),
                 lambda: system.calc_system_evolution(use_locks=use_locks, lock_tolerance=lock_tolerance)):
        with pytest.raises(ValueError, match="lock_tolerance must be finite and in"):
            call()


def test_defaults_come_from_the_dynamics_config(restore_config):
    assert TidalPy.config["dynamics"] == {"use_spin_locks": False, "spin_lock_tolerance": 1.0e-5}
    world = cpl_moon()
    system = system_with(world, 0.3, 1.5 - 0.3 * TOLERANCE)
    assert not system.calc_world_evolution(world)["spin_locked"]
    TidalPy.config["dynamics"]["use_spin_locks"] = True
    held = system.calc_world_evolution(world)
    assert held["spin_locked"]
    assert held["spin_lock_bracket"] == pytest.approx(1.5 + 0.7 * TOLERANCE, rel=1.0e-15)
    TidalPy.config["dynamics"]["spin_lock_tolerance"] = 0.1 * TOLERANCE
    assert not system.calc_world_evolution(world)["spin_locked"]
    # A call's own arguments win over the configuration.
    assert system.calc_world_evolution(world, lock_tolerance=TOLERANCE)["spin_locked"]
    assert not system.calc_world_evolution(world, use_locks=False, lock_tolerance=TOLERANCE)["spin_locked"]


def test_dynamics_config_keys_are_known_and_checked():
    packaged = get_packaged_config()
    assert find_unknown_config_keys({"dynamics": {"use_spin_locks": True, "spin_lock_tolerance": 1.0e-4}}, packaged) \
        == []
    assert find_unknown_config_keys({"dynamics": {"spin_lock_tolerence": 1.0e-4}}, packaged) == \
        ["dynamics.spin_lock_tolerence"]
    assert find_invalid_config_values({"dynamics": {"use_spin_locks": True, "spin_lock_tolerance": 1.0e-4}},
                                      packaged) == []
    for value in (0.0, 0.25, -1.0, math.inf):
        invalid = find_invalid_config_values({"dynamics": {"spin_lock_tolerance": value}}, packaged)
        assert [key for key, _ in invalid] == ["dynamics.spin_lock_tolerance"]
    assert [key for key, _ in find_invalid_config_values({"dynamics": {"use_spin_locks": 1}}, packaged)] == \
        ["dynamics.use_spin_locks"]
