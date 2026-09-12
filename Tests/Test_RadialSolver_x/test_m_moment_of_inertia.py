"""Tests for the two moment-of-inertia measures on RadialSolverSolution.

`moi_factor` is the conventional factor moi / (M R^2), which is exactly 0.4 for a uniform sphere. The separate
`moi_sphere_ratio` measures the same moment of inertia against a uniform sphere of equal mass and radius, so it
is exactly 1 there and falls below 1 as mass concentrates toward the center.
"""
import math

import numpy as np
import pytest

from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.RadialSolver_x import homogeneous_love_numbers, radial_solver
from TidalPy.rheology_x import Elastic, Maxwell
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers

RADIUS = 6.0e6
DENSITY = 4000.0
FREQUENCY = 1.0e-5


@pytest.fixture(scope="module")
def uniform():
    solution = homogeneous_love_numbers(RADIUS, DENSITY, 6.0e10 + 0.0j, FREQUENCY, num_slices=80)
    assert solution.success, solution.message
    return solution


def test_uniform_sphere_values(uniform):
    """A uniform sphere has moi = 0.4 M R^2 exactly."""
    expected_mass = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY
    assert uniform.mass == pytest.approx(expected_mass, rel=1e-6)
    assert uniform.moi == pytest.approx(0.4 * expected_mass * RADIUS**2, rel=1e-6)
    assert uniform.moi_factor == pytest.approx(0.4, rel=1e-6)
    assert uniform.moi_sphere_ratio == pytest.approx(1.0, rel=1e-6)


def test_the_two_measures_are_consistent(uniform):
    """The sphere ratio is 2.5 times the factor, by definition."""
    assert uniform.moi_sphere_ratio == pytest.approx(2.5 * uniform.moi_factor, rel=1e-12)
    assert uniform.moi_factor == pytest.approx(uniform.moi / (uniform.mass * uniform.radius**2), rel=1e-12)


def test_dense_core_lowers_both_measures():
    """Concentrating mass toward the center drops the factor below 0.4 and the ratio below 1."""
    build_data = build_rs_input_homogeneous_layers(
        RADIUS, FREQUENCY,
        density_tuple=(9000.0, 3000.0),
        static_bulk_modulus_tuple=(2.0e11, 1.0e11),
        static_shear_modulus_tuple=(1.0e11, 5.0e10),
        bulk_viscosity_tuple=(1.0e30, 1.0e30),
        shear_viscosity_tuple=(1.0e22, 1.0e21),
        layer_type_tuple=("solid", "solid"),
        layer_is_static_tuple=(False, False),
        layer_is_incompressible_tuple=(False, False),
        shear_rheology_model_tuple=Maxwell(),
        bulk_rheology_model_tuple=Elastic(),
        radius_fraction_tuple=(0.5, 1.0),
        slices_tuple=(30, 30))
    solution = radial_solver(*build_data, degree_l=2, solve_for=("tidal",))
    assert solution.success, solution.message
    assert 0.0 < solution.moi_factor < 0.4
    assert 0.0 < solution.moi_sphere_ratio < 1.0
    assert solution.moi_sphere_ratio == pytest.approx(2.5 * solution.moi_factor, rel=1e-12)
