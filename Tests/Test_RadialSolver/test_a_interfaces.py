"""Interface solver: what it fills, and that it keeps independent solutions independent.

The upper-layer values themselves are checked against the classic solver's in test_comparison/test_compare_interfaces.
"""
import numpy as np
import pytest

from TidalPy.RadialSolver.interfaces.interfaces import solve_upper_y_at_interface

STATIC_LIQUID_DENSITY = 7600.
INTERFACE_GRAVITY = 2.7
G_TO_USE = 6.67430e-11
SOLID, LIQUID = 0, 1
MAX_NUM_Y = 6

# Each layer kind's (solutions, stored ys): a solid carries 3 of y1 to y6, a dynamic liquid 2 of y1, y2, y5, y6, and a
# static liquid 1 of y5 and y7.
LAYER_KIND_SIZES = {
    (SOLID, True): (3, 6),
    (SOLID, False): (3, 6),
    (LIQUID, False): (2, 4),
    (LIQUID, True): (1, 2),
}
# Lower-layer solutions drawn at random are linearly independent, unlike hand-picked rows such as a row and its
# negative, which made some pairings collapse to exactly zero and test nothing.
RANDOM_SEED = 2024
PAIRINGS = [(lower_type, lower_static, upper_type, upper_static)
            for lower_type in (SOLID, LIQUID) for lower_static in (True, False)
            for upper_type in (SOLID, LIQUID) for upper_static in (True, False)]
PAIRING_IDS = [f"{'solid' if lt == SOLID else 'liquid'}_{'static' if ls else 'dynamic'}"
               f"__to__{'solid' if ut == SOLID else 'liquid'}_{'static' if us else 'dynamic'}"
               for lt, ls, ut, us in PAIRINGS]


def _independent_lower_y(layer_type, is_static):
    num_solutions, num_ys = LAYER_KIND_SIZES[(layer_type, is_static)]
    rng = np.random.default_rng(RANDOM_SEED)
    lower_y = np.full((num_solutions, MAX_NUM_Y), np.nan, dtype=np.complex128)
    lower_y[:, :num_ys] = rng.normal(size=(num_solutions, num_ys)) + 1j * rng.normal(size=(num_solutions, num_ys))
    return lower_y


def _solve(lower_type, lower_static, upper_type, upper_static):
    upper_y = np.full((3, MAX_NUM_Y), np.nan, dtype=np.complex128)
    solve_upper_y_at_interface(
        _independent_lower_y(lower_type, lower_static), upper_y, lower_type, lower_static, upper_type, upper_static,
        INTERFACE_GRAVITY, STATIC_LIQUID_DENSITY, G_TO_USE)
    return upper_y


@pytest.mark.parametrize("lower_type, lower_static, upper_type, upper_static", PAIRINGS, ids=PAIRING_IDS)
def test_the_interface_fills_every_solution_of_the_upper_layer(lower_type, lower_static, upper_type, upper_static):
    """Exactly the upper layer's solutions and stored ys are filled, all finite; the rest stay NaN."""
    num_solutions, num_ys = LAYER_KIND_SIZES[(upper_type, upper_static)]
    upper_y = _solve(lower_type, lower_static, upper_type, upper_static)
    assert np.all(np.isfinite(upper_y[:num_solutions, :num_ys]))
    assert np.sum(~np.isnan(upper_y)) == num_solutions * num_ys


@pytest.mark.parametrize("lower_type, lower_static, upper_type, upper_static", PAIRINGS, ids=PAIRING_IDS)
def test_independent_lower_solutions_give_independent_upper_ones(lower_type, lower_static, upper_type, upper_static):
    """The upper layer starts with as many independent solutions as it carries, so no pairing loses one."""
    num_solutions, num_ys = LAYER_KIND_SIZES[(upper_type, upper_static)]
    upper_y = _solve(lower_type, lower_static, upper_type, upper_static)[:num_solutions, :num_ys]
    assert np.linalg.matrix_rank(upper_y) == num_solutions
