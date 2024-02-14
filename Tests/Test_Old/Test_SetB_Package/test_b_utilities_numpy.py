import numpy as np

import TidalPy
from TidalPy.utilities.numpy_helper.array_other import find_nearest, neg_array_for_log_plot, normalize_dict, value_np_cleanup


def test_normalize_dict():
    """ Tests the normalize_dict function """

    # Test with floats
    test_dict = {
        'a': 102.0,
        'b': 121.2,
        'c': -12.0
        }

    result = normalize_dict(test_dict, pass_negatives=False, new_max=1.0, new_min=0.0)

    assert type(result) == dict
    assert type(result['a']) in [float, np.float64]
    assert type(result['b']) in [float, np.float64]
    assert type(result['c']) in [float, np.float64]

    np.testing.assert_allclose(result['b'], 1.0)
    np.testing.assert_allclose(result['c'] + 1.0, 1.0)
    assert result['b'] > result['a'] > result['c']

    result = normalize_dict(test_dict, pass_negatives=False, new_max=3.2, new_min=-10.0)

    np.testing.assert_allclose(result['b'], 3.2)
    np.testing.assert_allclose(result['c'], -10.0)
    assert result['b'] > result['a'] > result['c']


def test_neg_array_for_plot():
    """ Test function neg_array_for_log_plot """

    x = np.linspace(-10, 10, 10)
    y_p, y_n = neg_array_for_log_plot(x)

    assert np.all(y_p[~np.isnan(y_p)] >= 0.0)
    assert np.all(y_n[~np.isnan(y_n)] >= 0.0)
    assert np.all(np.isnan(y_p[x < 0.0]))
    assert np.all(np.isnan(y_n[x >= 0.0]))