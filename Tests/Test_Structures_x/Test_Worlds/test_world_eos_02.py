"""Robustness of the world EOS solve and its dense readouts.

Pins four behaviors: a radius a rounding step past a layer end reads that end rather than unwritten memory, a
radius clearly outside the body reads NaN, a structure with no hydrostatic solution reports failure, and a world
whose layers do not reach its radius is refused. Also checks that a radial solution released after the structure
changed is refused.
"""

import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos import make_material_eos
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.rheology_x.rheology import Maxwell, Elastic

PLANET_RADIUS = 1.0e6       # [m]
CORE_RADIUS   = 0.5e6       # [m]
CORE_DENSITY  = 6000.0      # [kg m-3]
MANTLE_DENSITY = 5000.0     # [kg m-3]
FORCING_FREQUENCY = 2.0 * math.pi / 86400.0  # [rad s-1]


def _two_layer_world(world_radius=PLANET_RADIUS):
    mass = (4.0 / 3.0) * math.pi * PLANET_RADIUS ** 3 * MANTLE_DENSITY
    world = LayeredWorld("two_layer", world_radius, mass)
    for index, (radius_inner, radius_outer, name, density) in enumerate(
            [(0.0, CORE_RADIUS, "core", CORE_DENSITY), (CORE_RADIUS, PLANET_RADIUS, "mantle", MANTLE_DENSITY)]):
        layer = PhysicsLayer(name, index, radius_inner, radius_outer, 0.0)
        layer.set_eos(ConstantDensityEOS(
            reference_density=density,
            shear_modulus_static=5.0e10,
            bulk_modulus_static=1.0e11))
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e19}))
        layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e19}))
        layer.set_shear_rheology(Maxwell())
        layer.set_bulk_rheology(Elastic())
        world.add_layer(layer)
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
    """A layer asked just above its top reads its top; well above it reads NaN, never another radius's value."""
    core = solved_world.core
    at_top = core.get_gravity(CORE_RADIUS)
    solved_world.get_gravity(0.8 * PLANET_RADIUS)  # a different read first, which a stale buffer would echo
    assert core.get_gravity(CORE_RADIUS * (1.0 + 1.0e-12)) == at_top
    assert math.isnan(core.get_gravity(1.2 * CORE_RADIUS))


@pytest.mark.parametrize("bulk_modulus, radius", ((1.0e11, 1.5e7), (1.0e10, 6.4e6)))
def test_no_hydrostatic_solution_is_a_failure(bulk_modulus, radius):
    """A Vinet sphere this large has no hydrostatic solution; the capped solve must say so, not succeed."""
    world = LayeredWorld("unsolvable", radius, 1.0e24)
    layer = PhysicsLayer("L", 0, 0.0, radius, 0.0)
    layer.set_eos(make_material_eos("vinet", {
        "reference_density_kg_m3": 4000.0,
        "reference_bulk_modulus_pa": bulk_modulus,
        "bulk_modulus_derivative": 4.0}))
    world.add_layer(layer)
    result = world.solve_eos()
    assert result["success"] is False
    assert result["max_iters_hit"] is True
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
