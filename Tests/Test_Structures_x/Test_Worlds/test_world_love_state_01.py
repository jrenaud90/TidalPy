"""What a world reports after its Love solves.

A solution released from a world reports its planet scalars and layer radii in SI, whatever units the cached solve
ran in, and the averaged modulus of an analytic Love method describes the last analytic solve only.
"""
import math

import numpy as np

from TidalPy.structures_x.configs import build_world

_FREQUENCY = 2.0 * math.pi / (12.42 * 3600.0)   # [rad s-1] M2


def _earth():
    world = build_world("earth_simple")
    world.solve_eos()
    return world


def test_released_solution_reports_si_scalars():
    world = _earth()
    world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2)
    solution = world.release_radial_solution()
    assert math.isclose(solution.radius, world.radius, rel_tol=1.0e-12)
    assert math.isclose(solution.mass, world.planet_mass_eos, rel_tol=1.0e-9)
    assert math.isclose(solution.surface_gravity, world.surface_gravity_eos, rel_tol=1.0e-9)
    assert math.isclose(solution.moi, world.planet_moi_eos, rel_tol=1.0e-9)
    assert 3.0e11 < solution.central_pressure < 4.0e11          # [Pa]
    assert math.isclose(solution.moi_factor, 0.3307, rel_tol=1.0e-3)
    upper_radii = solution.layer_upper_radius_array
    assert math.isclose(upper_radii[-1], world.radius, rel_tol=1.0e-12)
    assert np.all(np.diff(upper_radii) > 0.0)


def test_effective_modulus_clears_after_a_radial_solve():
    world = _earth()
    world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2, love_method="homogeneous")
    assert np.isfinite(world.love_effective_shear_modulus.real)
    assert np.isfinite(world.love_tidal_volume)
    world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2, love_method="radial_solver")
    assert world.love_success
    assert np.isnan(world.love_effective_shear_modulus.real)
    assert np.isnan(world.love_tidal_volume)
