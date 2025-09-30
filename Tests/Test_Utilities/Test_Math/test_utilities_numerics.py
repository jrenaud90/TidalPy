import pytest

from math import isclose as py_isclose
from math import nan

from TidalPy.utilities.math.numerics import isclose


@pytest.mark.parametrize('a', (0.0, -55.222226, 10., 1.0e200, -1.0e200, 1.0e-9, 1.0e-5, -1.0e-9, -1.0e-5, 1.4562221e34, 1.0e-14, 0.98e-14))
@pytest.mark.parametrize('b', (0.0, -55.222222, 10., 1.0e200, -1.0e200, 1.0e-9, 1.0e-5, -1.0e-9, -1.0e-5, 1.4582221e34, 1.0e-14, 0.98e-14))
@pytest.mark.parametrize('rtol', (3.0, 1.0e-9, 1.0e-3))
@pytest.mark.parametrize('atol', (3.0, 1.0e-3, 0.0))
def test_math_isclose(a, b, rtol, atol):

    pyresult  = py_isclose(a, b, rel_tol=rtol, abs_tol=atol)
    tpyresult = isclose(a, b, rtol, atol)

    assert pyresult == tpyresult

def test_math_isclose_nans():

    assert isclose(nan, 10.0) == False
    assert isclose(10.0, nan) == False

