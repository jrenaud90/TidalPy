"""The world's heat sources beside radiogenics: the tidal heat source (the last calc_tides heating, spread over the
heated layers by every later EOS solve) and prescribed heating (a power or a specific rate per layer)."""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Material import Material, Phase
from TidalPy.Rheology.rheology import Elastic, Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Utilities.conversions import orbital_motion2semi_a

_RADIUS = 1.0e6
_CORE_RADIUS = 4.0e5
_DENSITY = 3300.0
_CONDUCTIVITY = 4.0
_HEAT_CAPACITY = 1200.0
_ORBITAL_MOTION = 2.0 * math.pi / 86400.0
_ECCENTRICITY = 0.05


def _rock(viscosity=1.0e19):
    """Incompressible rock of one density, with a Maxwell-relaxing shear modulus."""
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": _DENSITY, "bulk_modulus_pa": 1.0e11},
        shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": viscosity},
        thermal_conductivity_w_mk=_CONDUCTIVITY,
        heat_capacity_j_kgk=_HEAT_CAPACITY))


def _sphere(use_heating=True):
    """One conducting layer of uniform density."""
    world = BaseWorld("sphere", _RADIUS, 4.0 / 3.0 * math.pi * _DENSITY * _RADIUS**3)
    world.add_layer(Layer("body", 0, 0.0, _RADIUS, material=_rock(), temperature=1500.0, cooling="conduction",
                          use_heating=use_heating))
    return world


def _two_layers():
    """A conducting core under a conducting Maxwell mantle that is heated."""
    mass = 4.0 / 3.0 * math.pi * _DENSITY * _RADIUS**3
    world = BaseWorld("two", _RADIUS, mass)
    world.add_layer(Layer("core", 0, 0.0, _CORE_RADIUS, material=_rock(1.0e24), temperature=1800.0,
                          cooling="conduction", shear_rheology=Elastic()))
    world.add_layer(Layer("mantle", 1, _CORE_RADIUS, _RADIUS, material=_rock(), temperature=1500.0,
                          cooling="conduction", use_heating=True, is_static=False, shear_rheology=Maxwell(),
                          bulk_rheology=Elastic()))
    return world


def _thermal_solve(world):
    result = world.solve_eos(G_to_use=G, solve_temperature=True, surface_temperature=300.0)
    assert result["success"], result["message"]
    assert result["thermal_converged"]
    return result


# ======================================================================================================================
# Prescribed heating
# ======================================================================================================================
def test_a_prescribed_specific_rate_heats_like_a_uniform_source():
    """A uniform sphere heated at s W/kg carries L(r) = (4/3) pi r^3 rho s up its conducting lower half."""
    rate = 5.0e-9
    world = _sphere()
    world.set_prescribed_heating("body", specific_rate=rate)
    assert world.prescribed_heating == {"body": {"specific_rate": rate}}
    result = _thermal_solve(world)
    mass = world.planet_mass_eos
    assert result["layer_heating_prescribed"][0] == pytest.approx(rate * mass, rel=1.0e-12)
    assert result["layer_heating"][0] == pytest.approx(rate * mass, rel=1.0e-12)
    assert result["layer_heating_radiogenic"] == [0.0]
    assert result["layer_heating_tidal"] == [0.0]
    for radius in (0.1 * _RADIUS, 0.3 * _RADIUS, 0.45 * _RADIUS):
        assert world.get_heat_flow(radius) == pytest.approx(4.0 / 3.0 * math.pi * radius**3 * _DENSITY * rate,
                                                            rel=1.0e-8)
    assert world.get_heating(0.5 * _RADIUS) == pytest.approx(_DENSITY * rate, rel=1.0e-12)
    np.testing.assert_allclose(world.get_heating(np.array([0.2, 0.8]) * _RADIUS), _DENSITY * rate, rtol=1.0e-12)


def test_a_prescribed_power_is_what_the_layer_receives():
    """A power is spread by mass, and the layer receives it whatever mass the solve gives the layer; the temperature
    rate takes it with the heat flows."""
    power = 2.0e12
    world = _two_layers()
    world.set_prescribed_heating(1, power=power)
    result = _thermal_solve(world)
    assert result["layer_heating_prescribed"] == pytest.approx([0.0, power], rel=1.0e-8)
    mantle_mass = world.mantle.mass
    rate = result["layer_temperature_rate"][1]
    expected = ((result["layer_heat_flow_in"][1] - result["layer_heat_flow_out"][1] + result["layer_heating"][1])
                / (mantle_mass * _HEAT_CAPACITY + result["layer_latent_capacity"][1]))
    assert rate == pytest.approx(expected, rel=1.0e-12)


def test_prescribed_heating_settings():
    world = _two_layers()
    assert _thermal_solve(world)["layer_heating_prescribed"] == [0.0, 0.0]
    world.set_prescribed_heating("mantle", power=1.0e12)
    # The EOS solve reads it.
    assert not world.eos_solved
    assert world.prescribed_heating == {"mantle": {"power": 1.0e12}}
    with pytest.raises(ValueError, match="not both"):
        world.set_prescribed_heating("mantle", power=1.0, specific_rate=1.0)
    with pytest.raises(ValueError, match="no layer named"):
        world.set_prescribed_heating("crust", power=1.0)
    with pytest.raises(ValueError, match="no layer at index"):
        world.set_prescribed_heating(5, power=1.0)
    with pytest.raises(ValueError, match="integer index"):
        world.set_prescribed_heating(1.7, power=1.0)
    with pytest.raises(ValueError, match="finite"):
        world.set_prescribed_heating("mantle", power=math.inf)
    world.set_prescribed_heating("mantle")
    assert world.prescribed_heating == {}


def test_a_layer_without_use_heating_takes_no_heat():
    world = _sphere(use_heating=False)
    world.set_prescribed_heating("body", specific_rate=5.0e-9)
    result = world.solve_eos(G_to_use=G, solve_temperature=True, surface_temperature=300.0)
    assert result["success"], result["message"]
    assert result["layer_heating"] == [0.0]
    assert world.get_heating(0.5 * _RADIUS) == 0.0


# ======================================================================================================================
# The tidal heat source
# ======================================================================================================================
def _tidal_world():
    world = _two_layers()
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    _thermal_solve(world)
    semi_major_axis = orbital_motion2semi_a(_ORBITAL_MOTION, mass_trap1, world.mass)
    tides = (_ORBITAL_MOTION, _ORBITAL_MOTION, _ECCENTRICITY, 0.0, semi_major_axis, mass_trap1)
    world.calc_tides(*tides)
    return world, tides


def test_the_tidal_heat_source_takes_the_last_tides():
    world, _ = _tidal_world()
    mantle_heating = world.get_layer_tidal_heating(1)
    assert mantle_heating > 0.0
    assert world.tidal_heat_source == {"core": world.get_layer_tidal_heating(0), "mantle": mantle_heating}
    result = _thermal_solve(world)
    # Only the heated mantle takes its tides, and it takes exactly the heating calc_tides gave it.
    assert result["layer_heating_tidal"][0] == 0.0
    assert result["layer_heating_tidal"][1] == pytest.approx(mantle_heating, rel=1.0e-12)
    # The re-solve keeps the source, which is the outer loop of an evolution.
    assert world.tidal_heat_source["mantle"] == mantle_heating
    world.clear_tidal_heating()
    assert world.tidal_heat_source == {}
    assert _thermal_solve(world)["layer_heating_tidal"] == [0.0, 0.0]


def test_the_tidal_heating_follows_the_radial_profile():
    """With a radial-solver Love method the heating follows the shell power dP/dr that calc_tides integrated: at its
    Gauss-Legendre nodes 4 pi r^2 h(r) is dP/dr times one factor for the whole layer."""
    world, tides = _tidal_world()
    mantle_heating = world.get_layer_tidal_heating(1)
    nodes, _ = np.polynomial.legendre.leggauss(16)
    radii = 0.5 * (_RADIUS + _CORE_RADIUS) + 0.5 * (_RADIUS - _CORE_RADIUS) * nodes
    shell_power = world.calc_3d_tides(*tides, radii=radii, latitude_summed=True, longitude_summed=True)["heating"]
    _thermal_solve(world)
    ratio = 4.0 * math.pi * radii**2 * world.get_heating(radii) / shell_power
    np.testing.assert_allclose(ratio, ratio[0], rtol=1.0e-9)
    assert ratio[0] == pytest.approx(1.0, rel=1.0e-2)
    # The shell power is far from uniform across the mantle, so the profile is what the test checks.
    assert shell_power.max() > 2.0 * shell_power.min()


def test_the_temperature_rate_takes_the_latest_tides():
    """A step of an evolution is solve_eos, calc_tides, then the rate: the rate takes the tides just computed, with no
    second solve, and the first step already has them."""
    world, _ = _tidal_world()
    result = world._build_eos_result()
    mantle_heating = world.get_layer_tidal_heating(1)
    assert result["layer_heating_tidal"] == [0.0, 0.0]
    expected = ((result["layer_heat_flow_in"][1] - result["layer_heat_flow_out"][1] + mantle_heating)
                / (world.mantle.mass * _HEAT_CAPACITY + result["layer_latent_capacity"][1]))
    assert world.calc_layer_temperature_rate("mantle") == pytest.approx(expected, rel=1.0e-12)
    assert result["layer_temperature_rate"][1] == pytest.approx(expected, rel=1.0e-12)


def test_a_failed_calc_tides_forgets_the_tidal_heat_source():
    world, tides = _tidal_world()
    assert world.tidal_heat_source
    world.set_prescribed_heating("mantle", power=1.0)   # forgets the solve
    with pytest.raises(RuntimeError, match="EOS solved first"):
        world.calc_tides(*tides)
    assert world.tidal_heat_source == {}


def test_a_layer_that_is_not_tidal_takes_no_tidal_heat():
    """A layer with use_tides off responds elastically: it dissipates nothing, so its heat leaves the total instead of
    moving into the tidal layers."""
    world, tides = _tidal_world()
    total = world.get_tidal_heating()
    world.core.use_tides = False
    world.calc_tides(*tides)
    assert world.get_layer_tidal_heating(0) == 0.0
    # This core is elastic, so it had no dissipation to lose: the total is unchanged and the mantle carries all of it.
    # (test_world_use_tides_radial_01 turns off a dissipating layer.)
    assert world.get_tidal_heating() == pytest.approx(total, rel=1.0e-12)
    assert world.get_layer_tidal_heating(1) == pytest.approx(total, rel=1.0e-12)
    assert world.tidal_heat_source == {"core": 0.0, "mantle": world.get_layer_tidal_heating(1)}


def test_without_a_per_layer_split_the_tides_are_spread_by_mass():
    """A radial-solver tide with layer_tidal_heating off has no per-layer split, so the tidal heat source spreads the
    total over the tidal layers by mass instead of losing it."""
    world, tides = _tidal_world()
    world.set_tide_config(layer_tidal_heating=False)
    world.calc_tides(*tides)
    assert math.isnan(world.get_layer_tidal_heating(1))
    total = world.get_tidal_heating()
    masses = [world.core.mass, world.mantle.mass]
    source = world.tidal_heat_source
    assert source["core"] == pytest.approx(total * masses[0] / sum(masses), rel=1.0e-12)
    assert source["mantle"] == pytest.approx(total * masses[1] / sum(masses), rel=1.0e-12)
    result = _thermal_solve(world)
    assert result["layer_heating_tidal"][1] == pytest.approx(source["mantle"], rel=1.0e-8)


def test_the_per_layer_tidal_heating_matches_the_3d_volume_integral():
    """calc_tides takes each layer's share from the same integral calc_3d_tides gives with every axis summed."""
    world, tides = _tidal_world()
    per_layer = np.asarray(world.calc_3d_tides(
        *tides, latitude_summed=True, longitude_summed=True, radial_summed=True)["per_layer"], dtype=float)
    shares = np.array([world.get_layer_tidal_heating(0), world.get_layer_tidal_heating(1)]) / world.get_tidal_heating()
    np.testing.assert_allclose(shares, per_layer / per_layer.sum(), rtol=1.0e-10)


def test_tidal_heating_in_a_layer_reaching_the_center():
    """In a layer that reaches the center the heating density stays finite there, and the heat flow the solve
    integrates is the volume integral of the reported heating."""
    world = BaseWorld("tidal_sphere", _RADIUS, 4.0 / 3.0 * math.pi * _DENSITY * _RADIUS**3)
    world.add_layer(Layer("body", 0, 0.0, _RADIUS, material=_rock(), temperature=1500.0, cooling="conduction",
                          use_heating=True, is_static=False, shear_rheology=Maxwell(), bulk_rheology=Elastic()))
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    _thermal_solve(world)
    semi_major_axis = orbital_motion2semi_a(_ORBITAL_MOTION, mass_trap1, world.mass)
    world.calc_tides(_ORBITAL_MOTION, _ORBITAL_MOTION, _ECCENTRICITY, 0.0, semi_major_axis, mass_trap1)
    tidal_heating = world.get_layer_tidal_heating(0)
    result = _thermal_solve(world)
    assert result["layer_heating_tidal"][0] == pytest.approx(tidal_heating, rel=1.0e-12)
    # Inside the first Gauss-Legendre node (x = 0.0053, about 5.3 km) the density holds that node's value.
    near_center = world.get_heating(np.array([1.0, 1.0e3, 5.0e3]))
    assert np.all(np.isfinite(near_center))
    np.testing.assert_allclose(near_center, near_center[0], rtol=1.0e-12)
    radii = np.linspace(0.0, 0.4 * _RADIUS, 4001)
    integrand = 4.0 * math.pi * radii**2 * world.get_heating(radii)
    integral = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(radii))
    assert world.get_heat_flow(0.4 * _RADIUS) == pytest.approx(integral, rel=1.0e-6)


def test_get_heating_is_nan_outside_a_solved_world():
    world = _sphere()
    assert math.isnan(world.get_heating(0.5 * _RADIUS))
    _thermal_solve(world)
    assert math.isnan(world.get_heating(2.0 * _RADIUS))
    assert math.isnan(world.get_heating(-1.0))
