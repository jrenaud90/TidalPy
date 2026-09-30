"""System orbital and spin evolution: single-body and pair dissipation, energy balance, and skipped worlds."""
import gc
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Utilities.conversions import orbital_motion2semi_a
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Viscosity import make_viscosity
from TidalPy.Rheology.rheology import Maxwell, Elastic
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Dynamics import Spin, OrbitSolver

_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_SMA = orbital_motion2semi_a(_N, _HOST, _MASS)


def _mass(radius):
    return (4.0 / 3.0) * math.pi * radius ** 3 * _DENSITY


def _layered(name, radius, spin_frequency):
    """A homogeneous Maxwell body that dissipates tidally and carries a spin model."""
    mass = _mass(radius)
    world = LayeredWorld(name, radius, mass)
    layer = BaseLayer("mantle", 0, 0.0, radius, mass)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_SHEAR, bulk_modulus_static=_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    world.set_spin_model(Spin())
    world.solve_eos(G_to_use=G)
    world.set_spin_frequency(spin_frequency)
    return world


def _host():
    """A point-mass tidal host (its structure is irrelevant here)."""
    return StarWorld("host", 7.0e8, _HOST)


def _moon(spin_factor=1.5):
    return _layered("moon", _R, spin_factor * _N)


def _system(spin_factor=1.5, eccentricity=_ECC):
    system = System("test")
    system.add_world(_host())
    system.add_world(_moon(spin_factor), tidal_host=0, semi_major_axis=_SMA, eccentricity=eccentricity)
    return system


def _dual_system(
        sma=60.0 * _R,
        host_radius=2.0 * _R,
        world_radius=_R,
        host_spin=0.7,
        world_spin=1.5,
        eccentricity=_ECC):
    """Both the host and the orbiter dissipate; spins are multiples of the shared mean motion."""
    mean_motion = math.sqrt(G * (_mass(host_radius) + _mass(world_radius)) / sma ** 3)
    system = System("dual")
    system.add_world(_layered("host", host_radius, host_spin * mean_motion))
    system.add_world(
        _layered("orbiter", world_radius, world_spin * mean_motion),
        tidal_host=0,
        semi_major_axis=sma,
        eccentricity=eccentricity)
    return system


def _hostless_system():
    system = System("hostless")
    system.add_world(_moon(), semi_major_axis=_SMA, eccentricity=_ECC)
    return system


def test_world_evolution_evolved_flags():
    """An orbiting world evolves using the system's mean motion, orbit, and host mass."""
    ev = _system().calc_world_evolution("moon")
    assert ev["evolved"] is True
    assert ev["has_spin"] is True
    assert ev["world_index"] == 1
    assert math.isclose(ev["orbital_frequency"], _N, rel_tol=1e-9)
    assert math.isclose(ev["semi_major_axis"], _SMA, rel_tol=1e-12)
    assert math.isclose(ev["eccentricity"], _ECC, rel_tol=1e-12)
    assert math.isclose(ev["host_mass"], _HOST, rel_tol=1e-12)
    assert math.isclose(ev["target_mass"], _MASS, rel_tol=1e-12)
    assert ev["tidal_heating"] > 0.0


def test_world_evolution_drives_the_solve():
    """calc_world_evolution runs the world's tidal solve."""
    system = _system()
    moon = system["moon"]
    assert moon.tides_solved is False
    system.calc_world_evolution(moon)
    assert moon.tides_solved is True


def test_world_evolution_matches_standalone_engines():
    """The system's orbital and spin rates match the standalone engines fed the same tidal solve."""
    ev = _system(spin_factor=1.5).calc_world_evolution("moon")

    moon = _moon(spin_factor=1.5)
    moon.calc_tides(
        orbital_frequency=_N,
        spin_frequency=1.5 * _N,
        eccentricity=_ECC,
        obliquity=0.0,
        semi_major_axis=_SMA,
        host_mass=_HOST)
    dU_dM, dU_dw, _ = moon.get_tidal_potential_derivatives()
    orbit = OrbitSolver()
    da_ref = orbit.calc_da_dt(_N, _SMA, _ECC, _MASS, _HOST, dU_dM)
    de_ref = orbit.calc_de_dt(_N, _SMA, _ECC, _MASS, _HOST, dU_dM, dU_dw)
    dn_ref = orbit.calc_dn_dt(_N, _SMA, da_ref)
    dspin_ref = moon.calc_spin_derivative(_HOST)

    assert math.isclose(ev["da_dt"], da_ref, rel_tol=1e-9)
    assert math.isclose(ev["de_dt"], de_ref, rel_tol=1e-9)
    assert math.isclose(ev["dn_dt"], dn_ref, rel_tol=1e-9)
    assert math.isclose(ev["dspin_dt"], dspin_ref, rel_tol=1e-9)


@pytest.mark.parametrize("spin_factor", [1.5, 1.37, 0.5, 2.0])
def test_energy_balance(spin_factor):
    """Heating equals -(dE_orbit/dt + dE_spin/dt), with E_orbit = -G M m / 2a and E_spin = I spin^2 / 2."""
    ev = _system(spin_factor=spin_factor).calc_world_evolution("moon")
    assert ev["has_spin"] is True
    assert abs(ev["energy_residual"]) <= 1e-6 * abs(ev["tidal_heating"])
    dE_orbit = G * _MASS * _HOST / (2.0 * _SMA ** 2) * ev["da_dt"]
    dE_spin = ev["moment_of_inertia"] * ev["spin_frequency"] * ev["dspin_dt"]
    assert math.isclose(ev["dE_orbit_dt"], dE_orbit, rel_tol=1e-9)
    assert math.isclose(ev["dE_spin_dt"], dE_spin, rel_tol=1e-9)
    assert math.isclose(ev["tidal_heating"], -(dE_orbit + dE_spin), rel_tol=1e-6)


def test_circular_orbit_zero_de_dt():
    ev = _system(spin_factor=1.5, eccentricity=0.0).calc_world_evolution("moon")
    assert ev["evolved"] is True
    assert ev["de_dt"] == 0.0


def test_host_entry_not_evolved():
    ev = _system().calc_world_evolution("host")
    assert ev["evolved"] is False
    assert ev["da_dt"] == 0.0
    assert ev["de_dt"] == 0.0
    assert ev["has_spin"] is False


@pytest.mark.parametrize("orbit", [
    pytest.param(dict(tidal_host="host", semi_major_axis=None), id="no_semi_major_axis"),
    pytest.param(dict(semi_major_axis=_SMA, eccentricity=_ECC), id="no_tidal_host"),
])
def test_world_without_orbit_or_host_not_evolved(orbit):
    """A world needs both a tidal host and a semi-major axis about it to evolve."""
    system = System("skipped")
    system.add_world(_host())
    system.add_world(_moon(), **orbit)
    assert system.calc_world_evolution("moon")["evolved"] is False


def test_system_evolution_sweep():
    """The sweep returns one row per world: the host is skipped and the moon evolves."""
    system = _system()
    results = system.calc_system_evolution()
    assert isinstance(results, list)
    assert len(results) == system.num_worlds
    assert results[0]["world_index"] == 0 and results[0]["evolved"] is False
    assert results[1]["world_index"] == 1 and results[1]["evolved"] is True
    assert results[1]["has_spin"] is True
    assert abs(results[1]["energy_residual"]) <= 1e-6 * abs(results[1]["tidal_heating"])


def test_layerless_world_evolves_without_spin():
    """A layerless fixed-Q world evolves its orbit but, without a spin model, its spin terms are NaN."""
    companion_mass = 1.898e27
    sma = 1.0e10
    orbital_frequency = math.sqrt(G * (_HOST + companion_mass) / sma ** 3)

    system = System("layerless")
    system.add_world(_host())
    companion = StarWorld("companion", 5.0e8, companion_mass)
    companion.set_tide_model(make_tide("cpl", {"fixed_k": [0.03], "fixed_q": [1.0e6]}))
    companion.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    companion.set_spin_frequency(orbital_frequency)
    system.add_world(companion, tidal_host=0, semi_major_axis=sma, eccentricity=_ECC)

    ev = system.calc_world_evolution("companion")
    assert ev["evolved"] is True
    assert ev["has_spin"] is False
    # Torqued, but with no moment of inertia to turn the torque into a rate: unknown, not zero.
    assert math.isnan(ev["dspin_dt"])
    assert math.isnan(ev["dE_spin_dt"])
    assert math.isnan(ev["energy_residual"])
    assert ev["tidal_heating"] > 0.0
    assert np.isfinite(ev["da_dt"])
    assert np.isfinite(ev["de_dt"])
    assert np.isfinite(ev["dE_orbit_dt"])


@pytest.mark.parametrize("world_spin, host_spin", [(1.5, 0.7), (1.37, 1.2), (0.5, 2.0)])
def test_pair_evolution_energy_balance(world_spin, host_spin):
    """Both bodies dissipate; the combined and per-body energy balances hold and the rates add."""
    pair = _dual_system(world_spin=world_spin, host_spin=host_spin).calc_pair_evolution("orbiter")
    assert pair["evolved"] is True
    assert pair["world"]["tidal_heating"] > 0.0
    assert pair["host"]["tidal_heating"] > 0.0
    assert abs(pair["energy_residual"]) <= 1e-6 * abs(pair["tidal_heating_total"])
    assert abs(pair["world"]["energy_residual"]) <= 1e-6 * abs(pair["world"]["tidal_heating"])
    assert abs(pair["host"]["energy_residual"]) <= 1e-6 * abs(pair["host"]["tidal_heating"])
    assert math.isclose(pair["da_dt"], pair["world"]["da_dt"] + pair["host"]["da_dt"], rel_tol=1e-12)
    assert math.isclose(pair["de_dt"], pair["world"]["de_dt"] + pair["host"]["de_dt"], rel_tol=1e-12)
    assert math.isclose(
        pair["tidal_heating_total"], pair["world"]["tidal_heating"] + pair["host"]["tidal_heating"], rel_tol=1e-12)


def test_pair_world_matches_single_body():
    """The pair's orbiting-world contribution equals calc_world_evolution."""
    system = _dual_system()
    pair = system.calc_pair_evolution("orbiter")
    single = system.calc_world_evolution("orbiter")
    for key in ("da_dt", "de_dt", "dn_dt", "dspin_dt", "tidal_heating", "energy_residual"):
        assert math.isclose(pair["world"][key], single[key], rel_tol=1e-12, abs_tol=1e-30)


def test_pair_symmetry_identical_bodies():
    """Two identical bodies contribute equally, so the combined rate is twice each."""
    pair = _dual_system(host_radius=_R, world_radius=_R, host_spin=1.5, world_spin=1.5).calc_pair_evolution(
        "orbiter")
    assert math.isclose(pair["world"]["da_dt"], pair["host"]["da_dt"], rel_tol=1e-9)
    assert math.isclose(pair["world"]["tidal_heating"], pair["host"]["tidal_heating"], rel_tol=1e-9)
    assert math.isclose(pair["world"]["dspin_dt"], pair["host"]["dspin_dt"], rel_tol=1e-9)
    assert math.isclose(pair["da_dt"], 2.0 * pair["world"]["da_dt"], rel_tol=1e-12)


def test_pair_rigid_host_reduces_to_single_body():
    """A host with no tide model is rigid, so the pair reduces to the single-body result."""
    system = _system()
    pair = system.calc_pair_evolution("moon")
    single = system.calc_world_evolution("moon")
    assert pair["host"]["tidal_heating"] == 0.0
    assert pair["host"]["da_dt"] == 0.0
    assert pair["host"]["has_spin"] is False
    assert math.isclose(pair["da_dt"], single["da_dt"], rel_tol=1e-12)
    assert math.isclose(pair["tidal_heating_total"], single["tidal_heating"], rel_tol=1e-12)
    assert abs(pair["energy_residual"]) <= 1e-6 * abs(pair["tidal_heating_total"])


@pytest.mark.parametrize("build_system, world", [
    pytest.param(_dual_system, "host", id="host_entry"),
    pytest.param(_hostless_system, "moon", id="no_host"),
])
def test_pair_not_evolved(build_system, world):
    """A pair needs a world with a tidal host."""
    assert build_system().calc_pair_evolution(world)["evolved"] is False


def test_mutual_pair_rows_sum_to_the_pair_evolution():
    """Two worlds hosting each other each get a sweep row, and the rows add up to the pair's rates."""
    system = _dual_system()
    system.set_tidal_host("host", "orbiter")
    rows = system.calc_system_evolution()
    assert rows[0]["evolved"] is True and rows[1]["evolved"] is True
    pair = system.calc_pair_evolution("orbiter")
    assert math.isclose(rows[0]["da_dt"] + rows[1]["da_dt"], pair["da_dt"], rel_tol=1e-12)
    assert math.isclose(
        rows[0]["tidal_heating"] + rows[1]["tidal_heating"], pair["tidal_heating_total"], rel_tol=1e-12)
    # Seen from either member it is the same pair with the roles swapped.
    mirrored = system.calc_pair_evolution("host")
    assert mirrored["evolved"] is True
    assert mirrored["host_index"] == 1
    assert math.isclose(mirrored["da_dt"], pair["da_dt"], rel_tol=1e-12)
    assert math.isclose(mirrored["world"]["tidal_heating"], pair["host"]["tidal_heating"], rel_tol=1e-12)


def test_world_gets_its_tide_state_from_its_system():
    """get_tide_state is None until the world has a system, a host, and an orbit, then matches the evolution."""
    moon = _moon()
    assert moon.get_tide_state() is None
    system = System("test")
    system.add_world(_host())
    system.add_world(moon)
    assert moon.get_tide_state() is None
    system.set_tidal_host(moon, "host")
    assert moon.get_tide_state() is None
    system.set_semi_major_axis(moon, _SMA)
    system.set_eccentricity(moon, _ECC)

    state = moon.get_tide_state()
    assert state["orbital_frequency"] == system.calc_orbital_frequency(moon)
    assert state["semi_major_axis"] == _SMA and state["eccentricity"] == _ECC
    assert state["host_mass"] == _HOST
    assert state["spin_frequency"] == moon.spin_frequency
    assert state["obliquity"] == moon.obliquity
    evolution = system.calc_world_evolution(moon)
    moon.calc_tides(**state)
    assert moon.get_tidal_heating() == evolution["tidal_heating"]


def test_world_stops_asking_a_system_that_is_gone():
    moon = _moon()
    system = System("short_lived")
    system.add_world(_host())
    system.add_world(moon, tidal_host=0, semi_major_axis=_SMA, eccentricity=_ECC)
    assert moon.get_tide_state() is not None
    del system
    gc.collect()
    assert moon.get_tide_state() is None
