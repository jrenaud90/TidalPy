"""The surface-solve amplification and rcond diagnostics and the poorly conditioned warning (standalone and world)."""
import math

import numpy as np
import pytest

from TidalPy.Rheology import Maxwell
from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.RadialSolver.rs_solution import SEVERE_SURFACE_AMPLIFICATION
from TidalPy.Structures import build_world
from TidalPy.exceptions import SolutionFailedError

frequency = 2.0 * np.pi / (86400. * 1.0)
N = 10
radius_array = np.linspace(0.0, 6000.e3, N)
bulk_density = 3500.
density_array = bulk_density * np.ones_like(radius_array)
bulk_modulus_array = 1.0e11 * np.ones(N, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)
complex_shear_modulus_array = Maxwell().calc_complex_modulus_vectorize_modulus(shear_array, viscosity_array, frequency)
upper_radius_by_layer = np.asarray((radius_array[-1],))

WARNING_TEXT = "poorly conditioned"


def _solve(**kwargs):
    """A 1-layer dynamic incompressible solid with Kamata starts."""
    return radial_solver(
        radius_array,
        density_array,
        bulk_modulus_array,
        complex_shear_modulus_array,
        frequency,
        bulk_density,
        ('solid',),
        (False,),
        (True,),
        upper_radius_by_layer,
        starting_method="kamata",
        raise_on_fail=True,
        **kwargs)


def _run(starting_radius):
    """Degree-3 solve; a starting radius of a meter or less makes the surface constants cancel catastrophically."""
    return _solve(
        degree_l=3,
        solve_for=('tidal',),
        integration_method='RK45',
        integration_rtol=1.0e-7,
        integration_atol=1.0e-10,
        scale_rtols_bylayer_type=False,
        max_num_steps=5_000_000,
        expected_size=250,
        max_step=0,
        verbose=False,
        nondimensionalize=False,
        starting_radius=starting_radius,
        warnings=True)


def _solved_earth(**kwargs):
    world = build_world("earth_simple")
    world.solve_eos()
    return world, world.solve_love_numbers(frequency=1.0e-5, **kwargs)


def test_pathological_solve_records_amplification_and_warns(spdlog_text):
    """A 1 m starting radius at degree 3 (rcond near 1e-9) reports severe amplification and logs the warning."""
    solution = _run(starting_radius=1.0)
    assert solution.success
    assert solution.surface_solve_amplification > SEVERE_SURFACE_AMPLIFICATION
    assert WARNING_TEXT in spdlog_text()


def test_singular_start_fails():
    """A 0.1 m starting radius drops the rcond near 1e-13, below the default floor: the solve fails (k was off by
    2e-4 at rtol 1e-7 when it was accepted)."""
    with pytest.raises(SolutionFailedError, match="singular to working precision"):
        _run(starting_radius=0.1)


def test_healthy_solve_is_silent(spdlog_text):
    """The automatic starting radius stays below the warning thresholds and logs nothing."""
    solution = _run(starting_radius=0.0)
    assert solution.success
    assert 0.0 < solution.surface_solve_amplification < SEVERE_SURFACE_AMPLIFICATION
    assert WARNING_TEXT not in spdlog_text()


def test_world_level_amplification_property(spdlog_text):
    """A layered world records a healthy amplification after a solve and stays silent."""
    world, result = _solved_earth()
    assert result['success']
    assert 0.0 < world.love_surface_amplification < SEVERE_SURFACE_AMPLIFICATION
    assert WARNING_TEXT not in spdlog_text()


@pytest.mark.parametrize('warnings_flag', (True, False))
def test_amplification_is_recorded_whatever_the_warnings_flag(warnings_flag):
    """The diagnostic is always computed; the warnings flag only decides whether it is reported."""
    solution = _solve(degree_l=2, warnings=warnings_flag)
    assert solution.success
    assert 0.0 < solution.surface_solve_amplification < SEVERE_SURFACE_AMPLIFICATION

    world, result = _solved_earth(warnings=warnings_flag)
    assert result['success']
    assert 0.0 < world.love_surface_amplification < SEVERE_SURFACE_AMPLIFICATION


def test_world_level_rcond_property():
    """The world's rcond is NaN before a solve, in (1e-8, 1] after, and matches its released solution's."""
    world = build_world("earth_simple")
    world.solve_eos()
    assert math.isnan(world.love_surface_rcond)
    world.solve_love_numbers(frequency=1.0e-5)
    rcond = world.love_surface_rcond
    assert 1.0e-8 < rcond <= 1.0
    assert world.release_radial_solution().surface_solve_rcond == rcond


def test_world_level_rcond_is_nan_for_the_analytic_methods():
    """The homogeneous analytic method has no surface system, so its rcond is NaN."""
    world, _ = _solved_earth()
    assert world.love_surface_rcond > 0.0
    world.solve_love_numbers(frequency=1.0e-5, love_method="homogeneous")
    assert math.isnan(world.love_surface_rcond)
