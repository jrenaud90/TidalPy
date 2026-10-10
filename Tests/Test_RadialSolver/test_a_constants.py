"""RadialSolver compile-time constants."""
from TidalPy.RadialSolver.rs_constants import get_constants


def test_constants_values():
    """The RadialSolver size constants have their expected values."""
    constants = get_constants()
    assert constants['MAX_NUM_Y'] == 6
    assert constants['MAX_NUM_Y_REAL'] == 12
    assert constants['MAX_NUM_SOL'] == 3
