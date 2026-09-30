"""Propagation matrix: the solution continued between and below its grid radii, and the core seeds it accepts."""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.RadialSolver import radial_solver
from TidalPy.Rheology.rheology import Elastic, Maxwell
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Viscosity import make_viscosity

_RADIUS    = 6.0e6                        # [m]
_DENSITY   = 5500.0                       # [kg m-3]
_SHEAR     = 5.0e10 + 1.0e9j              # [Pa]
_FREQUENCY = 2.0 * np.pi / 86400.         # [rad s-1]


def _solve(num_slices, degree_l=2, solve_for=('tidal',), **kwargs):
    """A uniform static incompressible solid sphere with the propagation matrix."""
    radius_array = np.linspace(0.0, _RADIUS, num_slices)
    return radial_solver(
        radius_array,
        np.full(num_slices, _DENSITY),
        np.full(num_slices, 1.0e14 + 0.0j),
        np.full(num_slices, _SHEAR),
        _FREQUENCY,
        _DENSITY,
        ('solid',),
        (True,),
        (True,),
        np.asarray((_RADIUS,)),
        degree_l=degree_l,
        solve_for=solve_for,
        love_method='propagation_matrix',
        warnings=False,
        **kwargs)


def _analytic_k(degree_l):
    """k_l of a uniform incompressible sphere."""
    gravity = (4.0 / 3.0) * np.pi * G * _DENSITY * _RADIUS
    effective_rigidity = (2.0 * degree_l**2 + 4.0 * degree_l + 3.0) * _SHEAR / (degree_l * _DENSITY * gravity * _RADIUS)
    return (3.0 / (2.0 * (degree_l - 1.0))) / (1.0 + effective_rigidity)


# ======================================================================================================================
# The solution between and below the grid radii
# ======================================================================================================================
@pytest.mark.parametrize('degree_l', (2, 3, 5))
@pytest.mark.parametrize('solve_for', (('tidal',), ('loading',), ('tidal', 'loading')))
def test_off_grid_solution_does_not_depend_on_slices(degree_l, solve_for):
    """For a uniform body the continued solution is exact, so 12 and 401 slices agree between the coarse radii.

    With linear interpolation between the coarse radii y5 was off by about 2% halfway between two of 12 slices, and
    the 3D heating converged only as the inverse square of the slice count.
    """
    coarse = _solve(12, degree_l=degree_l, solve_for=solve_for)
    fine = _solve(401, degree_l=degree_l, solve_for=solve_for)
    assert coarse.success and fine.success
    # Off the coarse grid (its spacing is R / 11), including inside its first slice.
    radii = np.array([0.013, 0.05, 0.123, 0.318, 0.5, 0.777, 0.96, 1.0]) * _RADIUS
    for ytype_index in range(len(solve_for)):
        y_coarse = coarse.get_radial_solution_array(radii, ytype_index)
        y_fine = fine.get_radial_solution_array(radii, ytype_index)
        assert np.all(np.isfinite(y_coarse))
        scale = np.abs(y_fine).max(axis=0)
        np.testing.assert_allclose(y_coarse, y_fine, rtol=0.0, atol=1.0e-8 * scale.max())
        for y_index in range(6):
            np.testing.assert_allclose(
                y_coarse[:, y_index], y_fine[:, y_index], rtol=0.0, atol=1.0e-8 * scale[y_index])


@pytest.mark.parametrize('degree_l', (2, 3))
def test_solution_reaches_the_center(degree_l):
    """Below the first propagated slice the regular solution continues to r = 0: y1 ~ r^(l-1), y5 ~ r^l."""
    solution = _solve(25, degree_l=degree_l)
    assert solution.success
    assert np.all(np.isfinite(solution.result))
    y_center = solution.get_radial_solution(0.0)
    assert np.all(np.isfinite(y_center))
    y_small = solution.get_radial_solution(1.0e-3 * _RADIUS)
    y_double = solution.get_radial_solution(2.0e-3 * _RADIUS)
    # Leading power laws of the regular solution of a uniform sphere.
    assert abs(y_double[0]) / abs(y_small[0]) == pytest.approx(2.0**(degree_l - 1), rel=1.0e-4)
    assert abs(y_double[4]) / abs(y_small[4]) == pytest.approx(2.0**degree_l, rel=1.0e-4)
    assert abs(y_center[4]) == 0.0


def test_grid_values_are_unchanged():
    """The continued solution reproduces the grid at the grid radii and the surface Love number is exact."""
    solution = _solve(25)
    radii = np.linspace(0.0, _RADIUS, 25)
    y_dense = solution.get_radial_solution_array(radii[2:], 0)
    y_grid = solution.result[:, 2:].T
    np.testing.assert_allclose(y_dense, y_grid, rtol=1.0e-10, atol=1.0e-13 * np.abs(y_grid).max())
    assert complex(solution.k) == pytest.approx(_analytic_k(2), rel=1.0e-10)


# ======================================================================================================================
# Core seeds other than the regular solution
# ======================================================================================================================
@pytest.mark.parametrize('core_model', (1, 2, 3, 4))
def test_irregular_core_seed_with_manual_start_is_refused(core_model):
    """Core models 1 to 4 changed k2 by -11% to +8% at a manual 0.3 R start while reporting success."""
    solution = _solve(100, core_model=core_model, starting_radius=0.3 * _RADIUS)
    assert not solution.success
    assert solution.error_code == -22
    assert 'core_model' in solution.message


@pytest.mark.parametrize('core_model', (1, 2, 3, 4))
def test_irregular_core_seed_with_automatic_start_solves(core_model):
    """At the automatic starting radius their effect on k2 is a few 1e-6, and the solution below the start is NaN."""
    solution = _solve(100, core_model=core_model)
    assert solution.success, solution.message
    assert complex(solution.k) == pytest.approx(_analytic_k(2), rel=1.0e-4)
    assert np.all(np.isnan(solution.get_radial_solution(0.0)))


@pytest.mark.parametrize('starting_radius_fraction', (0.05, 0.3))
def test_regular_core_seed_is_exact_at_any_start(starting_radius_fraction):
    solution = _solve(100, core_model=0, starting_radius=starting_radius_fraction * _RADIUS)
    assert solution.success, solution.message
    assert complex(solution.k) == pytest.approx(_analytic_k(2), rel=1.0e-10)
    assert np.all(np.isfinite(solution.get_radial_solution(0.0)))


# ======================================================================================================================
# World path: 3D heating no longer depends on the slice count
# ======================================================================================================================
_WORLD_RADIUS = 1.8216e6
_WORLD_MASS   = 8.9319e22
_HOST         = 1.898e27
_SMA          = 4.217e8
_MEAN_MOTION  = math.sqrt(G * (_HOST + _WORLD_MASS) / _SMA**3)
_STATE        = (_MEAN_MOTION, 0.3 * _MEAN_MOTION, 0.0, 0.0, _SMA, _HOST)


def _world(slices_per_layer):
    """A uniform static incompressible Maxwell body solved with the propagation matrix."""
    density = _WORLD_MASS / ((4.0 / 3.0) * math.pi * _WORLD_RADIUS**3)
    world = BaseWorld("homogeneous", _WORLD_RADIUS, _WORLD_MASS)
    layer = BaseLayer("mantle", 0, 0.0, _WORLD_RADIUS, _WORLD_MASS)
    layer.set_eos(ConstantDensityEOS(reference_density=density, shear_modulus_static=5.0e10))
    layer.is_static = True
    layer.is_incompressible = True
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e17}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, love_method="propagation_matrix")
    world.solve_eos(G_to_use=G, slices_per_layer=slices_per_layer)
    return world


@pytest.mark.parametrize('slices_per_layer', (10, 100))
def test_world_3d_heating_matches_1d_at_any_slice_count(slices_per_layer):
    """The volume integral of the 3D heating equals the 1D heating (it was 4.8e-2 low at 10 slices, 9.9e-5 at 100)."""
    world = _world(slices_per_layer)
    world.calc_tides(*_STATE)
    total = world.calc_3d_tides(
        *_STATE, latitude_summed=True, longitude_summed=True, radial_summed=True, num_threads=1)["total"]
    assert math.isclose(total, world.get_tidal_heating(), rel_tol=1.0e-6)


def test_world_radial_functions_reach_the_center():
    """get_love_radial_y is finite below the first grid radius and the same at 10 and 100 slices."""
    frequency = 1.4 * _MEAN_MOTION
    values = []
    for slices_per_layer in (10, 100):
        world = _world(slices_per_layer)
        result = world.solve_love_numbers(frequency=frequency, love_method="propagation_matrix")
        assert result['success'], result['message']
        values.append(np.array([
            world.get_love_radial_y(fraction * _WORLD_RADIUS, 0, y_index)
            for fraction in (0.0, 0.01, 0.37, 1.0) for y_index in range(6)]))
    assert np.all(np.isfinite(values[0]))
    np.testing.assert_allclose(values[0], values[1], rtol=1.0e-8, atol=1.0e-8 * np.abs(values[1]).max())
