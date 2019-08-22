import numba
import numpy as np

import TidalPy


TidalPy.verbose_level = 0
TidalPy.logging_level = 0
TidalPy.use_disk = False

from TidalPy.performance import njit, find_factorial, use_numba


def test_numba():
    if use_numba:
        assert njit is numba.njit
    else:

        def dummy_func(x):
            return x

        assert njit(dummy_func) is dummy_func


def test_find_factorial():
    # Check on common values to be fed into this function
    values = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    results = [1., 0.9513508, 0.9181687, 0.8974707, 0.88726382, 0.886227, 0.8935153, 0.9086387, 0.9313838, 0.961766, 1.]
    results = np.asarray(results)
    calc_results = np.asarray([find_factorial(v) for v in values])
    np.testing.assert_allclose(results, calc_results, rtol=1.e-6)
