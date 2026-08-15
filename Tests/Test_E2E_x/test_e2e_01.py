"""End-to-end integration tests over the config -> build -> physics pipeline.

Each test starts from a schema-0.2.0 TOML fixture, builds the world or system through the public
config API, attaches a tide model, and runs a full calculation, checking that the integrated workflow
produces physically sensible results. These are coarse-grained integration checks; the fine-grained
numerics are covered by the per-module unit tests.
"""
import math
from pathlib import Path

import numpy as np
import pytest

from TidalPy.structures_x.configs import build_world, build_system, load_toml
from TidalPy.structures_x.system import System
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.worlds.gasgiant import GasGiantWorld
from TidalPy.Tides_x.classes import make_tide

FIXTURES = Path(__file__).parent / "fixtures"


def _world(filename):
    return build_world(load_toml(str(FIXTURES / filename)))


def _configure_tide(world, model, config):
    world.set_tide_model(make_tide(model, config))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)


# =====================================================================================================================
# Terrestrial world with an analytic constant-phase-lag tide
# =====================================================================================================================
def test_terrestrial_fixedq_e2e():
    world = _world("terrestrial.toml")
    assert isinstance(world, LayeredWorld)
    assert world.world_type == "terrestrial"

    _configure_tide(world, "fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]})
    n = 2.0e-5
    world.set_spin_frequency(n)
    world.calc_tides(n, n, 0.05, 0.0, 4.0e8, 6.0e24)
    heating = world.get_tidal_heating()
    assert np.isfinite(heating) and heating > 0.0

    # Leading-order dissipation scales as the square of the eccentricity.
    world.calc_tides(n, n, 0.10, 0.0, 4.0e8, 6.0e24)
    assert math.isclose(world.get_tidal_heating() / heating, 4.0, rel_tol=1e-6)


# =====================================================================================================================
# Terrestrial world with a Maxwell rheology solved through the radial solver
# =====================================================================================================================
def test_terrestrial_rheology_love_e2e():
    world = _world("terrestrial.toml")
    world.solve_eos()
    assert world.eos_solved

    _configure_tide(world, "rheology", None)
    world.solve_love_numbers(frequency_rad_s=1.0e-4)
    assert world.love_solved

    k2 = world.love_number_k
    assert 0.0 < k2.real < 1.5           # a physical potential Love number
    assert k2.imag <= 0.0                # dissipative (negative imaginary part)


# =====================================================================================================================
# Gas giant with an analytic constant-time-lag tide
# =====================================================================================================================
def test_gasgiant_fixedlag_e2e():
    world = _world("gasgiant.toml")
    assert isinstance(world, GasGiantWorld)
    assert world.world_type == "gasgiant"

    _configure_tide(world, "fixed_dt", {"fixed_k": [0.5], "fixed_dt": [0.5]})
    n = 1.2e-5
    world.set_spin_frequency(n)
    world.calc_tides(n, n, 0.05, 0.0, 7.0e9, 1.0e30)
    heating = world.get_tidal_heating()
    assert np.isfinite(heating) and heating > 0.0


# =====================================================================================================================
# Star-planet system: one orbital + spin evolution step
# =====================================================================================================================
def test_star_planet_system_e2e():
    system = build_system(str(FIXTURES / "system.toml"))
    assert isinstance(system, System)
    assert system.num_worlds == 2
    assert system.host.name == "sun" and system.star.name == "sun"

    planet = system["planet"]
    _configure_tide(planet, "fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]})
    planet.set_spin_frequency(system.calc_orbital_frequency(planet))

    ev = system.calc_world_evolution(planet)
    assert ev["da_dt"] < 0.0                                  # tidal decay
    assert ev["de_dt"] < 0.0                                  # circularization
    assert np.isfinite(ev["tidal_heating"]) and ev["tidal_heating"] > 0.0
    # Energy is conserved: heating balances the orbital (and spin) energy loss.
    assert abs(ev["energy_residual"]) < 1.0e-6 * abs(ev["dE_orbit_dt"])


def test_star_planet_system_binary_roundtrip_e2e(tmp_path):
    system = build_system(str(FIXTURES / "system.toml"))
    path = str(tmp_path / "e2e_system.tpyb")
    system.save_binary(path)

    loaded = System()
    loaded.load_binary(path)
    assert loaded.name == system.name
    assert [w.name for w in loaded] == [w.name for w in system]
    assert math.isclose(loaded.get_semi_major_axis("planet"), system.get_semi_major_axis("planet"))
    assert math.isclose(loaded.calc_insolation_flux("planet"), system.calc_insolation_flux("planet"),
                        rel_tol=1e-12)
