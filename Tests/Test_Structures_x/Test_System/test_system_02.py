"""System orbital + spin evolution (single-body tidal dissipation).

The system drives each orbiting world's global tidal solve from its own state (mean motion from
Kepler's third law, spin + obliquity from the world, eccentricity + semi-major axis from the orbit
about the host, host mass from the host world), then turns the tidal-potential derivatives into the
orbital rates (da/dt, de/dt, dn/dt) via the orbital rate engine and the world's spin rate (dspin/dt)
via its attached spin model. Only the orbiting world dissipates; the host is a point mass.

These tests check the evolved/skipped bookkeeping, agreement with the standalone rate engines, the
system-level energy balance ``heating = -(dE_orbit/dt + dE_spin/dt)``, the circular-orbit limit, the
no-spin (layerless) path, and the whole-system sweep.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.utilities.conversions import orbital_motion2semi_a
from TidalPy.structures_x.system import System
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.rheology_x.rheology import Maxwell, Elastic
from TidalPy.Tides_x.classes.tide import make_tide
from TidalPy.dynamics_x import Spin, OrbitSolver


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


def _host():
    """A point-mass tidal host (a star with the host mass; its own structure is irrelevant here)."""
    return StarWorld("host", 7.0e8, _HOST)


def _moon(spin_factor=1.5, eccentricity=_ECC):
    """A homogeneous Maxwell moon that dissipates tidally and carries a spin model."""
    moon = LayeredWorld("moon", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS,
                         shear_modulus_static_pa=_SHEAR, bulk_modulus_static_pa=_BULK)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    moon.add_layer(layer)
    moon.set_tide_model(make_tide("rheology"))
    moon.set_tide_config(min_degree_l=2, max_degree_l=2,
                         eccentricity_truncation=3, obliquity_truncation=0)
    moon.set_spin_model(Spin())
    moon.solve_eos(G_to_use=G)
    moon.set_spin_frequency(spin_factor * _N)
    return moon


def _system(spin_factor=1.5, eccentricity=_ECC):
    system = System("test")
    system.add_world(_host(), is_host=True)
    system.add_world(_moon(spin_factor, eccentricity), semi_major_axis=_SMA, eccentricity=eccentricity)
    return system


# =====================================================================================================================
# Core: evolved bookkeeping + state
# =====================================================================================================================
def test_world_evolution_evolved_flags():
    system = _system()
    ev = system.calc_world_evolution("moon")
    assert ev["evolved"] is True
    assert ev["has_spin"] is True
    assert ev["world_index"] == 1
    # The state used comes from the system: mean motion from Kepler, a/e from the orbit, host mass.
    assert math.isclose(ev["orbital_frequency"], _N, rel_tol=1e-9)
    assert math.isclose(ev["semi_major_axis"], _SMA, rel_tol=1e-12)
    assert math.isclose(ev["eccentricity"], _ECC, rel_tol=1e-12)
    assert math.isclose(ev["host_mass"], _HOST, rel_tol=1e-12)
    assert math.isclose(ev["target_mass"], _MASS, rel_tol=1e-12)
    assert ev["tidal_heating"] > 0.0


def test_world_evolution_drives_the_solve():
    """calc_world_evolution runs the world's tidal solve (the world starts unsolved)."""
    system = _system()
    moon = system["moon"]
    assert moon.tides_solved is False
    system.calc_world_evolution(moon)
    assert moon.tides_solved is True


def test_world_evolution_matches_standalone_engines():
    """The system's orbital + spin rates match the standalone engines fed the same tidal solve."""
    system = _system(spin_factor=1.5)
    ev = system.calc_world_evolution("moon")

    moon = _moon(spin_factor=1.5)
    moon.calc_tides(orbital_frequency=_N, spin_frequency=1.5 * _N, eccentricity=_ECC,
                    obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST)
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
    """System-level energy conservation: heating = -(dE_orbit/dt + dE_spin/dt)."""
    system = _system(spin_factor=spin_factor)
    ev = system.calc_world_evolution("moon")
    assert ev["has_spin"] is True
    # The residual (heating + dE_orbit/dt + dE_spin/dt) is ~0 relative to the heating.
    assert abs(ev["energy_residual"]) <= 1e-6 * abs(ev["tidal_heating"])
    # Cross-check the reported energy terms against E_orbit = -G M m / 2a, E_spin = 1/2 I spin^2.
    dE_orbit = G * _MASS * _HOST / (2.0 * _SMA ** 2) * ev["da_dt"]
    dE_spin = ev["moment_of_inertia"] * ev["spin_frequency"] * ev["dspin_dt"]
    assert math.isclose(ev["dE_orbit_dt"], dE_orbit, rel_tol=1e-9)
    assert math.isclose(ev["dE_spin_dt"], dE_spin, rel_tol=1e-9)
    assert math.isclose(ev["tidal_heating"], -(dE_orbit + dE_spin), rel_tol=1e-6)


def test_circular_orbit_zero_de_dt():
    system = _system(spin_factor=1.5, eccentricity=0.0)
    ev = system.calc_world_evolution("moon")
    assert ev["evolved"] is True
    assert ev["de_dt"] == 0.0


# =====================================================================================================================
# Skipped worlds (host, no orbit, no host)
# =====================================================================================================================
def test_host_entry_not_evolved():
    system = _system()
    ev = system.calc_world_evolution("host")
    assert ev["evolved"] is False
    assert ev["da_dt"] == 0.0
    assert ev["de_dt"] == 0.0
    assert ev["has_spin"] is False


def test_world_without_orbit_not_evolved():
    """A world whose semi-major axis about the host is unset cannot be evolved."""
    system = System("no_orbit")
    system.add_world(_host(), is_host=True)
    system.add_world(_moon(), semi_major_axis=None)   # orbit left unset
    ev = system.calc_world_evolution("moon")
    assert ev["evolved"] is False


def test_no_host_system_not_evolved():
    """With no host designated, a world cannot be evolved (no tidal companion)."""
    system = System("hostless")
    system.add_world(_moon(), semi_major_axis=_SMA, eccentricity=_ECC)
    ev = system.calc_world_evolution("moon")
    assert ev["evolved"] is False


# =====================================================================================================================
# Whole-system sweep
# =====================================================================================================================
def test_system_evolution_sweep():
    system = _system()
    results = system.calc_system_evolution()
    assert isinstance(results, list)
    assert len(results) == system.num_worlds
    # Host entry (index 0) is skipped; the moon (index 1) evolves.
    assert results[0]["world_index"] == 0 and results[0]["evolved"] is False
    assert results[1]["world_index"] == 1 and results[1]["evolved"] is True
    assert results[1]["has_spin"] is True
    assert abs(results[1]["energy_residual"]) <= 1e-6 * abs(results[1]["tidal_heating"])


# =====================================================================================================================
# No-spin path (a layerless orbiting world dissipating analytically)
# =====================================================================================================================
def test_layerless_world_evolves_without_spin():
    """A layerless world (a star with a fixed-Q tide) evolves its orbit but carries no spin model."""
    companion_mass = 1.898e27
    sma = 1.0e10
    orbital_frequency = math.sqrt(G * (_HOST + companion_mass) / sma ** 3)

    system = System("layerless")
    system.add_world(_host(), is_host=True)
    companion = StarWorld("companion", 5.0e8, companion_mass)
    companion.set_tide_model(make_tide("cpl", {"fixed_k": [0.03], "fixed_q": [1.0e6]}))
    companion.set_tide_config(min_degree_l=2, max_degree_l=2,
                              eccentricity_truncation=2, obliquity_truncation=0)
    companion.set_spin_frequency(orbital_frequency)
    system.add_world(companion, semi_major_axis=sma, eccentricity=_ECC)

    ev = system.calc_world_evolution("companion")
    assert ev["evolved"] is True
    assert ev["has_spin"] is False
    assert ev["dspin_dt"] == 0.0
    assert ev["dE_spin_dt"] == 0.0
    assert ev["tidal_heating"] > 0.0
    assert np.isfinite(ev["da_dt"])
    assert np.isfinite(ev["de_dt"])
    # The orbital-energy term is populated even without a spin model; the spin term is not.
    assert np.isfinite(ev["dE_orbit_dt"])
