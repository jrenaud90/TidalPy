"""End-to-end: schema-0.2.0 TOML fixture to built world or system to a full tidal calculation."""
import math
from pathlib import Path

import numpy as np
import pytest

from TidalPy.Structures.configs import build_world, build_system, load_toml
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.worlds.gasgiant import GasGiantWorld
from TidalPy.Tides.classes import make_tide

FIXTURES = Path(__file__).parent / "fixtures"


def _world(filename):
    return build_world(load_toml(str(FIXTURES / filename)))


def _configure_tide(world, model, config):
    world.set_tide_model(make_tide(model, config))
    # Truncation 2 keeps the heating through e^2, so a synchronous body's dissipation is exactly the e^2 term.
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)


@pytest.mark.parametrize(
    "filename, world_class, world_type, model, config, mean_motion, semi_major_axis, host_mass",
    [
        pytest.param(
            "terrestrial.toml",
            BaseWorld,
            "terrestrial",
            "fixed_q",
            {"fixed_k": [0.3], "fixed_q": [100.0]},
            2.0e-5,
            4.0e8,
            6.0e24,
            id="terrestrial_fixed_q"),
        pytest.param(
            "gasgiant.toml",
            GasGiantWorld,
            "gasgiant",
            "fixed_dt",
            {"fixed_k": [0.5], "fixed_dt_s": [0.5]},
            1.2e-5,
            7.0e9,
            1.0e30,
            id="gasgiant_fixed_dt"),
    ])
def test_analytic_tide_e2e(
        filename,
        world_class,
        world_type,
        model,
        config,
        mean_motion,
        semi_major_axis,
        host_mass):
    """A fixture world with an analytic tide heats, and the heating scales as e^2."""
    world = _world(filename)
    assert isinstance(world, world_class)
    assert world.world_type == world_type

    _configure_tide(world, model, config)
    world.set_spin_frequency(mean_motion)

    def heating_at(eccentricity):
        world.calc_tides(
            mean_motion,
            mean_motion,
            eccentricity,
            0.0,
            semi_major_axis,
            host_mass)
        return world.get_tidal_heating()

    heating = heating_at(0.05)
    assert np.isfinite(heating) and heating > 0.0
    assert math.isclose(heating_at(0.10) / heating, 4.0, rel_tol=1e-12)


def test_terrestrial_rheology_love_e2e():
    """A terrestrial fixture solved through the radial solver gives a physical, dissipative k2."""
    world = _world("terrestrial.toml")
    world.solve_eos()
    assert world.eos_solved

    _configure_tide(world, "rheology", None)
    world.solve_love_numbers(frequency=1.0e-4)
    assert world.love_solved

    k2 = world.love_number_k
    assert 0.0 < k2.real < 1.5
    assert k2.imag <= 0.0


def test_star_planet_system_e2e():
    """One evolution step of a star-planet system decays and circularizes the orbit and conserves energy."""
    system = build_system(str(FIXTURES / "system.toml"))
    assert isinstance(system, System)
    assert system.num_worlds == 2
    assert system.get_tidal_host("planet").name == "sun" and system.star.name == "sun"

    planet = system["planet"]
    _configure_tide(planet, "fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]})
    planet.set_spin_frequency(system.calc_orbital_frequency(planet))

    evolution = system.calc_world_evolution(planet)
    assert evolution["da_dt"] < 0.0
    assert evolution["de_dt"] < 0.0
    assert np.isfinite(evolution["tidal_heating"]) and evolution["tidal_heating"] > 0.0
    assert abs(evolution["energy_residual"]) < 1.0e-6 * abs(evolution["dE_orbit_dt"])


def test_star_planet_system_binary_roundtrip_e2e(tmp_path):
    """A system saved to binary reloads with the same worlds and orbit."""
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
