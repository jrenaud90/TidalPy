import math

import pytest

from TidalPy.RadialSolver_x.love import LoveNumbers


def test_love_numbers_initialization():
    """Test storing and retrieving complex Love numbers."""
    k = 1 + 2j
    h = 3 + 4j
    l = 5 + 6j
    ln = LoveNumbers(k, h, l)
    assert ln.k == k
    assert ln.h == h
    assert ln.l == l


def test_love_numbers_Q_normal():
    """Test Q factor calculation: Q = -|x| / Im(x)."""
    val = 3 + 4j
    ln = LoveNumbers(val, val, val)
    expected_Q = -abs(val) / val.imag  # -5/4 = -1.25
    assert math.isclose(ln.Q_k, expected_Q)
    assert math.isclose(ln.Q_h, expected_Q)
    assert math.isclose(ln.Q_l, expected_Q)


def test_love_numbers_Q_zero_imaginary():
    """Test Q when imaginary part is zero -> infinity."""
    val = 5 + 0j
    ln = LoveNumbers(val, val, val)
    assert math.isinf(ln.Q_k)
    assert math.isinf(ln.Q_h)
    assert math.isinf(ln.Q_l)


def test_love_numbers_lag_normal():
    """Test phase lag: lag = atan(-Im(x) / Re(x))."""
    val = 3 + 4j
    ln = LoveNumbers(val, val, val)
    expected_lag = math.atan(-val.imag / val.real)
    assert math.isclose(ln.lag_k, expected_lag)
    assert math.isclose(ln.lag_h, expected_lag)
    assert math.isclose(ln.lag_l, expected_lag)


def test_love_numbers_lag_zero_imaginary():
    """Test lag when imaginary part is zero -> 0.0."""
    val = 5 + 0j
    ln = LoveNumbers(val, val, val)
    assert ln.lag_k == 0.0
    assert ln.lag_h == 0.0
    assert ln.lag_l == 0.0


def test_love_numbers_lag_zero_real():
    """Test lag when real part is zero -> pi/2."""
    val = 0 + 4j
    ln = LoveNumbers(val, val, val)
    assert math.isclose(ln.lag_k, math.pi / 2.0)
    assert math.isclose(ln.lag_h, math.pi / 2.0)
    assert math.isclose(ln.lag_l, math.pi / 2.0)
