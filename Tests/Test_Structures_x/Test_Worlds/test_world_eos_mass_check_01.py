"""An EOS solve fails when the structure it finds is far from the world's stated mass.

Layers with no hydrostatic structure near the stated mass can still meet the surface-pressure target on a collapsed
branch at an absurd central pressure; `[numerical] maximum_eos_mass_ratio` turns that into a failed solve.
"""
import math

import pytest

import TidalPy
from TidalPy.Material_x.eos import make_material_eos
from TidalPy.structures_x.configs import build_world
from TidalPy.structures_x.layers import SolidLiquidLayer
from TidalPy.structures_x.worlds import LayeredWorld

# A 0.69 Earth-mass planet whose core holds 23.6% of the mass inside 0.19 of the radius (about 1.7e5 kg m-3), with
# Birch-Murnaghan layers at 80% of those densities: the only surface-pressure root holds about 60 times the mass.
MASS = 0.692 * 5.972e24
RADIUS = 0.92 * 6.371e6
CORE_RADIUS = 0.19 * RADIUS
CORE_MASS = 0.236 * MASS


def _birch_murnaghan(density):
    return make_material_eos("birch_murnaghan", {
        "reference_density_kg_m3": density,
        "reference_bulk_modulus_pa": 1.0e11,
        "bulk_modulus_derivative": 4.0,
        "shear_modulus_static_pa": 3.7e10})


def _collapsing_world():
    core_density = CORE_MASS / ((4.0 / 3.0) * math.pi * CORE_RADIUS ** 3)
    mantle_density = (MASS - CORE_MASS) / ((4.0 / 3.0) * math.pi * (RADIUS ** 3 - CORE_RADIUS ** 3))
    world = LayeredWorld("collapsing", RADIUS, MASS)
    core = SolidLiquidLayer("core", 0, 0.0, CORE_RADIUS, 0.0, is_solid=False)
    core.set_eos(_birch_murnaghan(0.8 * core_density))
    world.add_layer(core)
    mantle = SolidLiquidLayer("mantle", 1, CORE_RADIUS, RADIUS, 0.0)
    mantle.set_eos(_birch_murnaghan(0.8 * mantle_density))
    world.add_layer(mantle)
    return world


def test_structure_far_from_the_stated_mass_fails():
    world = _collapsing_world()
    result = world.solve_eos()
    assert not result["success"]
    assert not world.eos_solved
    assert "maximum_eos_mass_ratio" in result["message"]
    assert "times its stated mass" in result["message"]


def test_the_limit_comes_from_the_config():
    numerical = TidalPy.config_x["numerical"]
    original = numerical["maximum_eos_mass_ratio"]
    try:
        numerical["maximum_eos_mass_ratio"] = 1.0e6
        TidalPy.constants.update_constants_x()
        assert _collapsing_world().solve_eos()["success"]
    finally:
        numerical["maximum_eos_mass_ratio"] = original
        TidalPy.constants.update_constants_x()


@pytest.mark.parametrize("name", ["io", "earth_simple", "jupiter_simple"])
def test_bundled_worlds_are_within_the_limit(name):
    world = build_world(name)
    assert world.solve_eos()["success"]
    ratio = world.planet_mass_eos / world.mass
    assert 1.0 / TidalPy.constants.maximum_eos_mass_ratio < ratio < TidalPy.constants.maximum_eos_mass_ratio
