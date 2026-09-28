"""Input validation that protects `radial_solver` from out-of-bounds memory access.

The solver writes one boundary-condition model per `solve_for` entry into a fixed 5-slot buffer, and reads and writes
every radial array over the radius array's length. Both are checked up front, even when `perform_checks` is False.
"""
import numpy as np
import pytest

from TidalPy.exceptions import ArgumentException
from TidalPy.RadialSolver import radial_solver

NUM_SLICES = 10
FREQUENCY = 2.0 * np.pi / (86400.0 * 2.0)


def _inputs(num_density=NUM_SLICES, num_bulk=NUM_SLICES, num_shear=NUM_SLICES):
    radius_array = np.linspace(0.0, 1.0e6, NUM_SLICES)
    density_array = np.full(num_density, 3000.0)
    bulk_array = np.full(num_bulk, 1.0e11, dtype=np.complex128)
    shear_array = np.full(num_shear, 5.0e10 + 1.0e6j, dtype=np.complex128)
    return radius_array, density_array, bulk_array, shear_array


@pytest.mark.parametrize("perform_checks", (True, False))
@pytest.mark.parametrize("short_array", ("density", "bulk", "shear"))
def test_array_length_mismatch_raises(short_array, perform_checks):
    kwargs = {f"num_{short_array}": NUM_SLICES - 3}
    radius_array, density_array, bulk_array, shear_array = _inputs(**kwargs)
    with pytest.raises(ArgumentException, match="must all match"):
        radial_solver(
            radius_array, density_array, bulk_array, shear_array, FREQUENCY, 3000.0,
            ("solid",), (False,), (False,), np.asarray([1.0e6]),
            perform_checks=perform_checks)


@pytest.mark.parametrize("perform_checks", (True, False))
def test_layer_tuple_length_mismatch_raises(perform_checks):
    radius_array, density_array, bulk_array, shear_array = _inputs()
    with pytest.raises(ArgumentException, match="equal lengths"):
        radial_solver(
            radius_array, density_array, bulk_array, shear_array, FREQUENCY, 3000.0,
            ("solid",), (False, False), (False,), np.asarray([1.0e6]),
            perform_checks=perform_checks)


@pytest.mark.parametrize("perform_checks", (True, False))
def test_upper_radius_length_mismatch_raises(perform_checks):
    radius_array, density_array, bulk_array, shear_array = _inputs()
    with pytest.raises(ArgumentException, match="number of layers"):
        radial_solver(
            radius_array, density_array, bulk_array, shear_array, FREQUENCY, 3000.0,
            ("solid",), (False,), (False,), np.asarray([5.0e5, 1.0e6]),
            perform_checks=perform_checks)


def test_empty_radius_array_raises():
    empty = np.zeros(0)
    empty_complex = np.zeros(0, dtype=np.complex128)
    with pytest.raises(ArgumentException, match="must not be empty"):
        radial_solver(
            empty, empty, empty_complex, empty_complex, FREQUENCY, 3000.0,
            ("solid",), (False,), (False,), np.asarray([1.0e6]))


def test_more_than_five_solve_for_entries_raises():
    radius_array, density_array, bulk_array, shear_array = _inputs()
    with pytest.raises(ArgumentException, match="at most 5"):
        radial_solver(
            radius_array, density_array, bulk_array, shear_array, FREQUENCY, 3000.0,
            ("solid",), (False,), (False,), np.asarray([1.0e6]),
            solve_for=("tidal",) * 6)


def test_too_few_slices_fails_without_python_checks():
    """With `perform_checks=False` the solver's own check rejects a layer with fewer than 5 slices.

    The check once compared the error code (`==`) instead of setting it, so the solve carried on regardless.
    """
    from TidalPy.exceptions import SolutionFailedError

    radius_array = np.linspace(0.0, 1.0e6, 4)
    with pytest.raises(SolutionFailedError, match="five layer slices"):
        radial_solver(
            radius_array, np.full(4, 3000.0), np.full(4, 1.0e11, dtype=np.complex128),
            np.full(4, 5.0e10 + 1.0e6j, dtype=np.complex128), FREQUENCY, 3000.0,
            ("solid",), (False,), (False,), np.asarray([1.0e6]),
            perform_checks=False, raise_on_fail=True)
