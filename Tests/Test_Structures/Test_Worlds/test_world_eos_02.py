"""Robustness of the world EOS solve and its readouts near and outside layer ends, and refusal of stale radial
solution releases."""

import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material import Material, Phase
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.layers import Layer
from TidalPy.Rheology.rheology import Maxwell, Elastic

PLANET_RADIUS = 1.0e6       # [m]
CORE_RADIUS   = 0.5e6       # [m]
CORE_DENSITY  = 6000.0      # [kg m-3]
MANTLE_DENSITY = 5000.0     # [kg m-3]
FORCING_FREQUENCY = 2.0 * math.pi / 86400.0  # [rad s-1]


def _material(density):
    """A constant-density Maxwell solid."""
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 1.0e11},
        shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e19},
        bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e19}))


def _two_layer_world(world_radius=PLANET_RADIUS):
    mass = (4.0 / 3.0) * math.pi * PLANET_RADIUS ** 3 * MANTLE_DENSITY
    world = BaseWorld("two_layer", world_radius, mass)
    for index, (radius_inner, radius_outer, name, density) in enumerate(
            [(0.0, CORE_RADIUS, "core", CORE_DENSITY), (CORE_RADIUS, PLANET_RADIUS, "mantle", MANTLE_DENSITY)]):
        world.add_layer(Layer(name, index, radius_inner, radius_outer, 0.0, _material(density),
                              shear_rheology=Maxwell(), bulk_rheology=Elastic()))
    return world


@pytest.fixture(scope="module")
def solved_world():
    world = _two_layer_world()
    world.solve_eos(G_to_use=G)
    return world


@pytest.mark.parametrize("relative_offset", (2.0e-16, 1.0e-12, 1.0e-9))
def test_radius_a_rounding_step_past_the_surface_reads_the_surface(solved_world, relative_offset):
    surface_gravity = solved_world.get_gravity(PLANET_RADIUS)
    assert np.isfinite(surface_gravity)
    assert solved_world.get_gravity(PLANET_RADIUS * (1.0 + relative_offset)) == surface_gravity
    assert solved_world.get_pressure(PLANET_RADIUS * (1.0 + relative_offset)) == \
        solved_world.get_pressure(PLANET_RADIUS)


@pytest.mark.parametrize("radius", (1.01 * PLANET_RADIUS, 1.5 * PLANET_RADIUS, -1.0))
def test_radius_outside_the_body_reads_nan(solved_world, radius):
    assert math.isnan(solved_world.get_gravity(radius))
    assert math.isnan(solved_world.get_pressure(radius))
    assert math.isnan(solved_world.get_density(radius))


def test_layer_read_past_its_own_top(solved_world):
    """A layer read just above its top gives its top value; well above it gives NaN."""
    core = solved_world.core
    at_top = core.get_gravity(CORE_RADIUS)
    solved_world.get_gravity(0.8 * PLANET_RADIUS)  # a different read first, which a stale buffer would echo
    assert core.get_gravity(CORE_RADIUS * (1.0 + 1.0e-12)) == at_top
    assert math.isnan(core.get_gravity(1.2 * CORE_RADIUS))


@pytest.mark.parametrize("bulk_modulus, radius", ((1.0e11, 1.5e7), (1.0e10, 6.4e6)))
def test_no_hydrostatic_solution_is_a_failure(bulk_modulus, radius):
    """A Vinet sphere this large has no hydrostatic solution, and the solve reports failure."""
    world = BaseWorld("unsolvable", radius, 1.0e24)
    material = Material(solid=Phase(eos={
        "model": "vinet",
        "reference_density_kg_m3": 4000.0,
        "reference_bulk_modulus_pa": bulk_modulus,
        "bulk_modulus_derivative": 4.0}))
    world.add_layer(Layer("L", 0, 0.0, radius, 0.0, material))
    result = world.solve_eos()
    assert result["success"] is False
    assert "no hydrostatic structure" in result["message"]
    assert world.eos_solved is False


def test_layers_must_reach_the_world_radius():
    world = _two_layer_world(world_radius=1.1 * PLANET_RADIUS)
    with pytest.raises(ValueError, match="layers must fill the world"):
        world.solve_eos(G_to_use=G)


def test_release_after_a_new_eos_solve_is_refused():
    world = _two_layer_world()
    world.solve_eos(G_to_use=G)
    world.solve_love_numbers(frequency=FORCING_FREQUENCY, degree_l=2, warnings=False)
    world.solve_eos(G_to_use=G)
    with pytest.raises(RuntimeError):
        world.release_radial_solution()
    world.solve_love_numbers(frequency=FORCING_FREQUENCY, degree_l=2, warnings=False)
    solution = world.release_radial_solution()
    assert solution.success


def test_release_after_an_analytic_solve_is_refused():
    world = _two_layer_world()
    world.solve_eos(G_to_use=G)
    world.solve_love_numbers(frequency=FORCING_FREQUENCY, degree_l=2, warnings=False)
    world.solve_love_numbers(frequency=FORCING_FREQUENCY, degree_l=2, love_method="homogeneous", warnings=False)
    with pytest.raises(RuntimeError):
        world.release_radial_solution()
