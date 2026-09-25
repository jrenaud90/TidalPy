"""A radial solution's ``result`` grid is sampled where the caller asked, and pairs with ``sample_radii``.

The standalone ``radial_solver`` fills ``result`` at the caller's own radius array, however unevenly spaced, and a
repeated interface radius takes the lower layer at its first copy and the upper layer at its second. A world's
released solution fills it on the solve's grid rather than leaving it empty. A radius above the surface has no
solution.
"""
import math

import numpy as np
import pytest

from TidalPy.RadialSolver_x import radial_solver
from TidalPy.rheology_x import Maxwell
from TidalPy.structures_x import build_world

_R = 1.0e6
_R_CORE = 0.4 * _R
_FREQUENCY = 2.0e-5


def _uneven_two_layer_inputs(core_is_liquid):
    """A core and mantle on unevenly spaced radii with unequal point counts; the interface radius appears twice."""
    core = _R_CORE * np.sort(np.concatenate(([0.0, 1.0], np.random.default_rng(3).uniform(0.05, 0.95, 38))))
    mantle = _R_CORE + (_R - _R_CORE) * np.linspace(0.0, 1.0, 31) ** 1.7
    radius = np.concatenate((core, mantle))
    n_core, n_mantle = core.size, mantle.size
    density = np.concatenate((np.full(n_core, 7000.0), np.full(n_mantle, 3300.0)))
    bulk = np.concatenate((np.full(n_core, 1.0e11), np.full(n_mantle, 1.0e11)))
    shear_core = np.full(n_core, 0.0 if core_is_liquid else 8.0e10, dtype=np.complex128)
    shear = np.concatenate((shear_core, Maxwell().calc_complex_modulus_vectorize_modulus(
        np.full(n_mantle, 6.0e10), np.full(n_mantle, 1.0e18), _FREQUENCY)))
    bulk = bulk.astype(np.complex128)
    layer_types = ("liquid" if core_is_liquid else "solid", "solid")
    return radius, density, _FREQUENCY, _R_CORE, _R, shear, bulk, layer_types


@pytest.mark.parametrize("core_is_liquid", [False, True])
def test_result_is_sampled_at_the_callers_radii(core_is_liquid):
    radius, density, frequency, r_core, r_surface, shear, bulk, layer_types = _uneven_two_layer_inputs(core_is_liquid)
    solution = radial_solver(
        radius, density, bulk, shear, frequency, 4000.0, layer_types, (True, False), (False, False),
        np.array([r_core, r_surface]), degree_l=2, solve_for=("tidal",))
    assert solution.success, solution.message
    result = solution.result
    assert result.shape == (6, radius.size)
    np.testing.assert_array_equal(solution.sample_radii(), radius)
    interface = int(np.flatnonzero(radius == r_core)[0])
    for i, r in enumerate(radius):
        if i in (interface, interface + 1) or r == 0.0:
            continue
        expected = solution.get_radial_solution(r)
        np.testing.assert_allclose(result[:, i], expected, rtol=1e-12, equal_nan=True)
    # The second copy of the interface radius is the mantle's: its y1..y4 are defined even over a static liquid core.
    assert np.all(np.isfinite(result[:4, interface + 1]))
    if core_is_liquid:
        assert np.all(np.isnan(result[:4, interface]))


def test_a_radius_above_the_surface_has_no_solution():
    radius, density, frequency, r_core, r_surface, shear, bulk, layer_types = _uneven_two_layer_inputs(False)
    solution = radial_solver(
        radius, density, bulk, shear, frequency, 4000.0, layer_types, (True, False), (False, False),
        np.array([r_core, r_surface]), degree_l=2, solve_for=("tidal",))
    for factor in (1.0000001, 1.5, 3.0):
        assert np.all(np.isnan(solution.get_radial_solution(factor * r_surface)))
    assert np.all(np.isfinite(solution.get_radial_solution(r_surface)))


def test_a_released_world_solution_has_a_filled_result():
    io = build_world("io")
    io.solve_eos()
    io.solve_love_numbers(frequency=4.11e-5, degree_l=2)
    released = io.release_radial_solution()
    result = released.result
    radii = released.sample_radii()
    assert result.shape[1] == radii.size
    assert np.count_nonzero(result) > 0.9 * result.size
    middle = radii.size // 2
    np.testing.assert_allclose(result[:, middle], released.get_radial_solution(radii[middle]), rtol=1e-12)
