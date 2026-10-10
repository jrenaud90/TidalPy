"""love_only: a shooting solve that keeps only what its Love numbers need, for a world and the standalone solver.

Without dense output each layer's solutions are read at its top from the last step CyRK stores, which ends on the
layer's upper radius, so the Love numbers and surface values match the dense solve's and the steps are the same. The
radial functions below the surface are not kept, and asking for them raises.
"""
import numpy as np
import pytest

from TidalPy.RadialSolver import radial_solver
from TidalPy.Rheology import Maxwell
from TidalPy.Structures import build_world

# Worlds with a solid-only interior, a dynamic liquid ocean, a static liquid core, and several solid layers.
_WORLDS = {"io": 4.1106e-5, "europa_dynamic": 2.0477e-5, "earth_prem": 1.4052e-4, "luna": 2.6617e-6}
_SURFACE_RADIUS = 1.0e6
_CORE_RADIUS = 0.4 * _SURFACE_RADIUS
_FREQUENCY = 2.0e-5


def _solved_world(name):
    world = build_world(name)
    world.solve_eos(raise_on_fail=True)
    return world


def _love(world, frequency, **kwargs):
    result = world.solve_love_numbers(frequency, 2, warnings=False, raise_on_fail=True, **kwargs)
    return np.array([result["love_number_k"], result["love_number_h"], result["love_number_l"]])


@pytest.mark.parametrize("name", tuple(_WORLDS))
def test_a_love_only_world_solve_gives_the_dense_solves_love_numbers_and_surface(name):
    world = _solved_world(name)
    frequency = _WORLDS[name]
    dense = _love(world, frequency, solve_for=("tidal", "loading"))
    dense_surface = [world.get_love_surface_y(ytype, y) for ytype in (0, 1) for y in range(6)]
    love_only = _love(world, frequency, solve_for=("tidal", "loading"), love_only=True)
    np.testing.assert_allclose(love_only, dense, rtol=1.0e-14, atol=0.0)
    np.testing.assert_allclose([world.get_love_surface_y(ytype, y) for ytype in (0, 1) for y in range(6)],
                               dense_surface, rtol=1.0e-14, atol=0.0, equal_nan=True)


def test_a_love_only_world_solve_keeps_no_radial_functions_until_a_dense_solve():
    world = _solved_world("io")
    _love(world, _WORLDS["io"], love_only=True)
    with pytest.raises(ValueError, match="love_only"):
        world.get_love_radial_y(0.5 * world.radius)
    released = world.release_radial_solution()
    assert released.love_only is True
    for read in (lambda: released.result, lambda: released.get_radial_solution(0.5 * world.radius),
                 lambda: released.get_radial_solution_array(np.array([0.5 * world.radius])),
                 lambda: released.get_result_by_ytype_name("tidal")):
        with pytest.raises(ValueError, match="love_only"):
            read()
    # Its interior and Love numbers are still there.
    assert np.isfinite(released.get_density(0.5 * world.radius))
    assert np.isfinite(released.k)

    _love(world, _WORLDS["io"])
    assert np.isfinite(world.get_love_radial_y(0.5 * world.radius))


def _standalone(love_only, **kwargs):
    radius = np.concatenate((np.linspace(0.0, _CORE_RADIUS, 40), np.linspace(_CORE_RADIUS, _SURFACE_RADIUS, 60)))
    density = np.concatenate((np.full(40, 7000.0), np.full(60, 3300.0)))
    bulk = np.full(radius.size, 1.0e11, dtype=np.complex128)
    shear = np.concatenate((np.zeros(40, dtype=np.complex128), Maxwell().calc_complex_modulus_vectorize_modulus(
        np.full(60, 6.0e10), np.full(60, 1.0e18), _FREQUENCY)))
    return radial_solver(radius, density, bulk, shear, _FREQUENCY, 4000.0, ("liquid", "solid"), (False, False),
                         (False, False), np.array([_CORE_RADIUS, _SURFACE_RADIUS]), warnings=False,
                         love_only=love_only, **kwargs)


def test_a_love_only_standalone_solve_gives_the_dense_solves_love_numbers_and_steps():
    dense = _standalone(False)
    love_only = _standalone(True)
    assert dense.success and love_only.success, love_only.message
    assert dense.love_only is False and love_only.love_only is True
    np.testing.assert_allclose(love_only.love, dense.love, rtol=1.0e-14, atol=0.0)
    np.testing.assert_array_equal(love_only.steps_taken, dense.steps_taken)
    with pytest.raises(ValueError, match="love_only"):
        love_only.result
    with pytest.raises(ValueError, match="love_only"):
        love_only.plot_ys(show_plot=False)
    assert love_only.get_density(0.5 * _SURFACE_RADIUS) == pytest.approx(dense.get_density(0.5 * _SURFACE_RADIUS))


def test_love_only_leaves_the_propagation_matrix_its_grid():
    """The propagation matrix fills its grid as it solves, so love_only changes nothing there."""
    radius = np.linspace(0.0, _SURFACE_RADIUS, 200)
    density = np.full(radius.size, 3300.0)
    bulk = np.full(radius.size, 1.0e11, dtype=np.complex128)
    shear = np.full(radius.size, 5.0e10, dtype=np.complex128)
    solutions = [radial_solver(radius, density, bulk, shear, _FREQUENCY, 3300.0, ("solid",), (True,), (True,),
                               np.array([_SURFACE_RADIUS]), love_method="propagation_matrix", warnings=False,
                               love_only=love_only) for love_only in (False, True)]
    assert solutions[1].love_only is False
    np.testing.assert_array_equal(solutions[1].result, solutions[0].result)
    np.testing.assert_array_equal(solutions[1].love, solutions[0].love)
