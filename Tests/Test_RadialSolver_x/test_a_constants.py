import pytest

from TidalPy.RadialSolver_x.rs_constants import get_constants


def test_constants_values():
    """Verify the RadialSolver_x constant values."""
    constants = get_constants()
    assert constants['MAX_NUM_Y'] == 6
    assert constants['MAX_NUM_Y_REAL'] == 12
    assert constants['MAX_NUM_SOL'] == 3
