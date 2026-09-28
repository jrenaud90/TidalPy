"""Where a radial solution's `result` grid is sampled, and that it pairs with `sample_radii`."""
import numpy as np
import pytest

from TidalPy.RadialSolver import radial_solver
from TidalPy.Rheology import Maxwell
from TidalPy.Structures import build_world

SURFACE_RADIUS = 1.0e6
CORE_RADIUS = 0.4 * SURFACE_RADIUS
FREQUENCY = 2.0e-5


def _solve_uneven_two_layer(core_is_liquid):
    """Solve a core and mantle on uneven radii with unequal point counts; the interface radius appears twice."""
    core = CORE_RADIUS * np.sort(np.concatenate(([0.0, 1.0], np.random.default_rng(3).uniform(0.05, 0.95, 38))))
    mantle = CORE_RADIUS + (SURFACE_RADIUS - CORE_RADIUS) * np.linspace(0.0, 1.0, 31) ** 1.7
    radius = np.concatenate((core, mantle))
    num_core, num_mantle = core.size, mantle.size
    density = np.concatenate((np.full(num_core, 7000.0), np.full(num_mantle, 3300.0)))
    bulk = np.full(radius.size, 1.0e11, dtype=np.complex128)
    shear_core = np.full(num_core, 0.0 if core_is_liquid else 8.0e10, dtype=np.complex128)
    shear = np.concatenate((shear_core, Maxwell().calc_complex_modulus_vectorize_modulus(
        np.full(num_mantle, 6.0e10), np.full(num_mantle, 1.0e18), FREQUENCY)))
    solution = radial_solver(
        radius,
        density,
        bulk,
        shear,
        FREQUENCY,
        4000.0,
        ("liquid" if core_is_liquid else "solid", "solid"),
        (True, False),
        (False, False),
        np.array([CORE_RADIUS, SURFACE_RADIUS]),
        degree_l=2,
        solve_for=("tidal",))
    return radius, solution


@pytest.mark.parametrize("core_is_liquid", [False, True])
def test_result_is_sampled_at_the_callers_radii(core_is_liquid):
    """`result` is filled at the caller's radii; a repeated interface takes the lower then the upper layer."""
    radius, solution = _solve_uneven_two_layer(core_is_liquid)
    assert solution.success, solution.message
    result = solution.result
    assert result.shape == (6, radius.size)
    np.testing.assert_array_equal(solution.sample_radii(), radius)
    interface = int(np.flatnonzero(radius == CORE_RADIUS)[0])
    for index, radius_value in enumerate(radius):
        if index in (interface, interface + 1) or radius_value == 0.0:
            continue
        expected = solution.get_radial_solution(radius_value)
        np.testing.assert_allclose(result[:, index], expected, rtol=1e-12, equal_nan=True)
    # The second copy of the interface is the mantle's, so y1..y4 are defined even over a static liquid core.
    assert np.all(np.isfinite(result[:4, interface + 1]))
    if core_is_liquid:
        assert np.all(np.isnan(result[:4, interface]))


def test_a_radius_above_the_surface_has_no_solution():
    """Radii above the surface return NaN; the surface itself is finite."""
    _, solution = _solve_uneven_two_layer(False)
    for factor in (1.0000001, 1.5, 3.0):
        assert np.all(np.isnan(solution.get_radial_solution(factor * SURFACE_RADIUS)))
    assert np.all(np.isfinite(solution.get_radial_solution(SURFACE_RADIUS)))


def test_a_released_world_solution_has_a_filled_result():
    """A world's released solution fills `result` on the solve's grid."""
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
