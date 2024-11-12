""" Tests for TidalPy's special math functions """

import numpy as np

import TidalPy


from TidalPy.utilities.math.numba_special import sqrt_neg


def test_sqrt_neg():
    """ Tests the special square root function that works with numba. """

    # Real inputs
    z_float_pos = 22.
    z_float_neg = -22.
    z_array_pos = np.linspace(1., 10., 5)
    z_array_neg = -np.linspace(1., 10., 5)

    # Test that the positive function works
    assert sqrt_neg(z_float_pos, is_real=True) == np.sqrt(z_float_pos)
    assert sqrt_neg(z_float_pos, is_real=False) == np.sqrt(z_float_pos)
    np.testing.assert_allclose(sqrt_neg(z_array_pos, is_real=True), np.sqrt(z_array_pos))
    np.testing.assert_allclose(sqrt_neg(z_array_pos, is_real=False), np.sqrt(z_array_pos))

    # Test that the expanded functionality works as expected
    assert sqrt_neg(z_float_neg, is_real=True) == np.lib.scimath.sqrt(z_float_neg)
    assert sqrt_neg(z_float_neg, is_real=False) == np.lib.scimath.sqrt(z_float_neg)
    np.testing.assert_allclose(sqrt_neg(z_array_neg, is_real=True), np.lib.scimath.sqrt(z_array_neg))
    np.testing.assert_allclose(sqrt_neg(z_array_neg, is_real=False), np.lib.scimath.sqrt(z_array_neg))

    # Now test how imaginary inputs work.
    z_float_pos = 22. + 2.j
    z_float_neg = -22. - 2.j
    z_array_pos = np.linspace(1., 10., 5) + 1.0j * np.linspace(1., 10., 5)
    z_array_neg = -np.linspace(1., 10., 5) - 1.0j * np.linspace(1., 10., 5)

    # Test that the positive function works
    np.testing.assert_allclose(sqrt_neg(z_float_pos, is_real=False), np.sqrt(z_float_pos))
    np.testing.assert_allclose(sqrt_neg(z_array_pos, is_real=False), np.sqrt(z_array_pos))

    # Test that the expanded functionality works as expected
    np.testing.assert_allclose(sqrt_neg(z_float_neg, is_real=False), np.lib.scimath.sqrt(z_float_neg))
    np.testing.assert_allclose(sqrt_neg(z_array_neg, is_real=False), np.lib.scimath.sqrt(z_array_neg))
