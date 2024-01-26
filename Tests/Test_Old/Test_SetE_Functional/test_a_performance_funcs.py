import numpy as np
from scipy.special import gamma

import TidalPy



def test_find_factorial():
    # Test Load
    from TidalPy.utilities.performance import find_factorial
    assert type(find_factorial(1.)) in [float, np.float64]

    # Test Floats
    assert find_factorial(0.) == 1.
    np.testing.assert_approx_equal(find_factorial(0.1), gamma(0.1 + 1))
    np.testing.assert_approx_equal(find_factorial(0.2), gamma(0.2 + 1))
    np.testing.assert_approx_equal(find_factorial(0.3), gamma(0.3 + 1))
    np.testing.assert_approx_equal(find_factorial(0.45), gamma(0.45 + 1))
    np.testing.assert_approx_equal(find_factorial(0.5), gamma(0.5 + 1))
    np.testing.assert_approx_equal(find_factorial(0.75), gamma(0.75 + 1))
    np.testing.assert_approx_equal(find_factorial(0.9), gamma(0.9 + 1))
    np.testing.assert_approx_equal(find_factorial(3.0), gamma(3.0 + 1))
    np.testing.assert_approx_equal(find_factorial(-1.0), gamma(-1.0 + 1))
    np.testing.assert_approx_equal(find_factorial(-2.0), gamma(-2.0 + 1))
    assert find_factorial(1.) == 1.
