import math
import pytest
from TidalPy.RadialSolver_x.love import LoveNumbers


def test_love_numbers_initialization():
    """Test that complex numbers are correctly passed and retrieved."""
    k_val = 1.0 + 2.0j
    h_val = 3.0 + 4.0j
    l_val = 5.0 + 6.0j
    
    ln = LoveNumbers(k_val, h_val, l_val)
    
    assert ln.k == k_val
    assert ln.h == h_val
    assert ln.l == l_val

def test_love_numbers_Q_normal():
    """Test Q value calculations for normal complex numbers."""
    # Using 3 + 4j because abs(3 + 4j) = 5.0
    k_val = 3.0 + 4.0j
    h_val = 3.0 + 4.0j
    l_val = 3.0 + 4.0j
    
    ln = LoveNumbers(k_val, h_val, l_val)
    
    # Q = -abs(x) / imag(x)
    # expected_Q = -5.0 / 4.0 = -1.25
    expected_Q = -1.25
    
    assert ln.Q_k == pytest.approx(expected_Q)
    assert ln.Q_h == pytest.approx(expected_Q)
    assert ln.Q_l == pytest.approx(expected_Q)

def test_love_numbers_Q_zero_imaginary():
    """Test Q value when imaginary part is zero (should return infinity)."""
    k_val = 5.0 + 0.0j
    h_val = 5.0 + 0.0j
    l_val = 5.0 + 0.0j
    
    ln = LoveNumbers(k_val, h_val, l_val)
    
    # Assuming TidalPyConstants::d_INF maps to float('inf') in Python
    assert math.isinf(ln.Q_k) and ln.Q_k > 0
    assert math.isinf(ln.Q_h) and ln.Q_h > 0
    assert math.isinf(ln.Q_l) and ln.Q_l > 0

def test_love_numbers_lag_normal():
    """Test phase lag calculations for normal complex numbers."""
    k_val = 3.0 + 4.0j
    h_val = 3.0 + 4.0j
    l_val = 3.0 + 4.0j
    
    ln = LoveNumbers(k_val, h_val, l_val)
    
    # lag = atan(-imag(x) / real(x))
    expected_lag = math.atan(-4.0 / 3.0)
    
    assert ln.lag_k == pytest.approx(expected_lag)
    assert ln.lag_h == pytest.approx(expected_lag)
    assert ln.lag_l == pytest.approx(expected_lag)

def test_love_numbers_lag_zero_imaginary():
    """Test phase lag when imaginary part is zero (should return 0.0)."""
    k_val = 5.0 + 0.0j
    h_val = 5.0 + 0.0j
    l_val = 5.0 + 0.0j
    
    ln = LoveNumbers(k_val, h_val, l_val)
    
    assert ln.lag_k == 0.0
    assert ln.lag_h == 0.0
    assert ln.lag_l == 0.0

def test_love_numbers_lag_zero_real():
    """Test phase lag when real part is zero (limit of arctan(inf) -> pi/2)."""
    k_val = 0.0 + 5.0j
    h_val = 0.0 + 5.0j
    l_val = 0.0 + 5.0j
    
    ln = LoveNumbers(k_val, h_val, l_val)
    
    expected_lag = math.pi / 2.0
    
    assert ln.lag_k == pytest.approx(expected_lag)
    assert ln.lag_h == pytest.approx(expected_lag)
    assert ln.lag_l == pytest.approx(expected_lag)