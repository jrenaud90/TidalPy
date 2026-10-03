"""Results from an earlier structure or a failed solve never read as the current one's.

The Love diagnostics of a solve on an earlier structure, a mass a layer took from a rejected solve, profiles of a failed
solve read through the calc_* getters, and a layer moved to impossible radii.
"""
import math

import pytest

from TidalPy.exceptions import SolutionFailedError
from TidalPy.Material import Material, Phase
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld

RADIUS = 1.0e6
DENSITY = 3000.0
# Far outside the factor of 10 that [numerical] maximum_eos_mass_ratio lets a solved structure differ from its world's
# stated mass, so the solve is rejected.
WRONG_MASS_FACTOR = 1.0e3
NO_SOLVE_MESSAGE = "No love-number solve has been run."


def _rock():
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": DENSITY, "bulk_modulus_pa": 1.0e11},
        shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e19}))


def _sphere(mass_factor=1.0):
    mass = 4.0 / 3.0 * math.pi * DENSITY * RADIUS**3
    world = BaseWorld("sphere", RADIUS, mass_factor * mass)
    world.add_layer(Layer("mantle", 0, 0.0, RADIUS, material=_rock()))
    return world


def test_love_diagnostics_do_not_outlive_a_new_structure():
    world = build_world("io")
    world.solve_eos()
    world.solve_love_numbers(frequency=4.1e-5)
    assert world.love_success and world.love_error_code == 0
    assert world.love_surface_amplification > 0.0
    world.solve_eos()
    assert not world.love_success
    assert world.love_error_code == -100
    assert world.love_message == NO_SOLVE_MESSAGE
    assert world.love_surface_amplification == 0.0
    assert math.isnan(world.love_surface_rcond)


def test_the_calc_getters_raise_for_a_failed_solve():
    world = _sphere(WRONG_MASS_FACTOR)
    with pytest.raises(SolutionFailedError, match="stated mass"):
        world.calc_density(0.5 * RADIUS)


def test_the_calc_getters_still_solve_a_good_world():
    assert _sphere().calc_density(0.5 * RADIUS) == pytest.approx(DENSITY)


def test_a_layer_refuses_inverted_radii():
    world = _sphere()
    layer = world.mantle
    with pytest.raises(ValueError, match="radius_inner <= radius_outer"):
        layer.set_radii(5.0e5, 3.0e5)
    assert (layer.radius_inner, layer.radius_outer) == (0.0, RADIUS)
    with pytest.raises(ValueError, match="finite radii"):
        layer.set_radii(-1.0, RADIUS)
