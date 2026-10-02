"""The low-level radial-solver wrappers refuse arrays too small for what the C++ writes or reads.

They run with bounds checks off and hand raw pointers to C++ that writes a fixed number of solutions and ys for the
layer kind, so a short array used to be overrun silently.
"""
import numpy as np
import pytest

from TidalPy.RadialSolver.boundaries.boundaries import apply_surface_bc
from TidalPy.RadialSolver.interfaces.interfaces import solve_upper_y_at_interface
from TidalPy.RadialSolver.interfaces.reversed import top_to_bottom_interface_bc
from TidalPy.RadialSolver.love import find_love
from TidalPy.RadialSolver.starting.driver import find_starting_conditions
from TidalPy.RadialSolver.starting.kamata import kamata_solid_dynamic_compressible

SOLID, LIQUID = 0, 1
G = 6.674e-11
NUM_YS = 6


def _y(rows, columns=NUM_YS):
    return np.ones((rows, columns), dtype=np.complex128)


def test_an_upper_buffer_without_a_row_per_solution_is_refused():
    with pytest.raises(ValueError, match="rows"):
        solve_upper_y_at_interface(_y(3), _y(2), SOLID, False, SOLID, False, 1.0, 3000.0, G)


def test_interface_buffers_of_different_widths_are_refused():
    with pytest.raises(ValueError, match="column count"):
        solve_upper_y_at_interface(_y(3), _y(3, 4), SOLID, False, SOLID, False, 1.0, 3000.0, G)


def test_a_max_num_y_that_disagrees_with_the_arrays_is_refused():
    with pytest.raises(ValueError, match="max_num_y"):
        solve_upper_y_at_interface(
            _y(3), _y(3), SOLID, False, SOLID, False, 1.0, 3000.0, G, max_num_y=4)


def test_an_interface_with_full_buffers_still_works():
    upper = _y(3)
    solve_upper_y_at_interface(_y(3), upper, SOLID, False, SOLID, False, 1.0, 3000.0, G)
    assert np.all(np.isfinite(upper))


def test_a_starting_buffer_with_too_few_rows_is_refused():
    with pytest.raises(ValueError, match="rows"):
        kamata_solid_dynamic_compressible(1.0e-5, 1.0e5, 3000.0, 1.0e11, 5.0e10, 2, G, _y(1))
    with pytest.raises(ValueError, match="rows"):
        find_starting_conditions(SOLID, False, False, True, 1.0e-5, 1.0e5, 3000.0, 1.0e11, 5.0e10, 2, G, _y(2))


def test_a_starting_buffer_with_too_few_columns_is_refused():
    with pytest.raises(ValueError, match="columns"):
        kamata_solid_dynamic_compressible(1.0e-5, 1.0e5, 3000.0, 1.0e11, 5.0e10, 2, G, _y(3, 4))


def test_a_short_constant_vector_is_refused():
    with pytest.raises(ValueError, match="constant_vector_view"):
        apply_surface_bc(np.zeros(1, dtype=np.complex128), np.zeros(3), _y(1, 2), 1.0, G, 0, LIQUID, True, False)
    with pytest.raises(ValueError, match="bc_view"):
        apply_surface_bc(np.zeros(3, dtype=np.complex128), np.zeros(3), _y(3), 1.0, G, 1, SOLID, False, False)


def test_a_short_layer_above_constant_vector_is_refused():
    with pytest.raises(ValueError, match="layer_above_constant_vector_view"):
        top_to_bottom_interface_bc(
            np.zeros(3, dtype=np.complex128), np.zeros(1, dtype=np.complex128), _y(3), 1.0, 1.0, 3000.0, 3000.0,
            SOLID, SOLID, False, False)


def test_find_love_needs_all_six_surface_values():
    with pytest.raises(ValueError, match="y1 to y6"):
        find_love(np.ones(3, dtype=np.complex128), 1.0)
