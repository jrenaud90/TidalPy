from math import isclose
from scipy.special import factorial2

import pytest

from TidalPy.utilities.math.special import double_factorial

@pytest.mark.parametrize('l', (0, 2, 3, 4, 5, 10, 20, 100))
def test_double_factorial(l):

    tpy_value = double_factorial(l)
    ref_value = factorial2(l)

    assert isclose(tpy_value, ref_value)
