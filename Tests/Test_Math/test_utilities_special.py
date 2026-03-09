import cmath
from math import isclose

import numpy as np
from scipy.special import factorial2
from scipy.special import spherical_jn, spherical_yn, spherical_in, spherical_kn

import pytest

from TidalPy.utilities.math.special import double_factorial
from TidalPy.utilities.math.special import (
    _spherical_jn_scalar, _spherical_jn_scalar_complex, _spherical_jn_d_scalar,
    _spherical_yn_scalar, _spherical_yn_scalar_complex, _spherical_yn_d_scalar,
    _spherical_in_scalar, _spherical_in_scalar_complex, _spherical_in_d_scalar,
    _spherical_kn_scalar, _spherical_kn_scalar_complex, _spherical_kn_d_scalar,
    spherical_jn_array, spherical_yn_array, spherical_in_array, spherical_kn_array,
    spherical_jn_complex_array, spherical_yn_complex_array, 
    spherical_in_complex_array, spherical_kn_complex_array
)


@pytest.mark.parametrize('l', (0, 2, 3, 4, 5, 10, 20, 100))
def test_double_factorial(l):

    tpy_value = double_factorial(l)
    ref_value = factorial2(l)

    assert isclose(tpy_value, ref_value)


# Shared test parameters
ORDERS = (0, 1, 2, 5)
REAL_VALUES = (0.5, 1.0, 2.5, 10.0)
COMPLEX_VALUES = (0.5 + 0.5j, 1.0 - 2.0j, 3.0 + 4.0j)
TOL = 1e-10

@pytest.mark.parametrize('n', ORDERS)
@pytest.mark.parametrize('x', REAL_VALUES)
def test_spherical_scalars_real(n, x):
    """Tests real scalar wrappers against SciPy."""
    
    # j_n
    expected_jn = spherical_jn(n, x)
    actual_jn = _spherical_jn_scalar(n, x)
    assert isinstance(actual_jn, float)
    assert isclose(actual_jn, expected_jn, rel_tol=TOL)

    # y_n
    expected_yn = spherical_yn(n, x)
    actual_yn = _spherical_yn_scalar(n, x)
    assert isinstance(actual_yn, float)
    assert isclose(actual_yn, expected_yn, rel_tol=TOL)

    # i_n
    expected_in = spherical_in(n, x)
    actual_in = _spherical_in_scalar(n, x)
    assert isinstance(actual_in, float)
    assert isclose(actual_in, expected_in, rel_tol=TOL)

    # k_n
    expected_kn = spherical_kn(n, x)
    actual_kn = _spherical_kn_scalar(n, x)
    assert isinstance(actual_kn, float)
    assert isclose(actual_kn, expected_kn, rel_tol=TOL)


@pytest.mark.parametrize('n', ORDERS)
@pytest.mark.parametrize('z', COMPLEX_VALUES)
def test_spherical_scalars_complex(n, z):
    """Tests complex scalar wrappers against SciPy."""
    
    # j_n complex
    expected_jn = spherical_jn(n, z)
    actual_jn = _spherical_jn_scalar_complex(n, z)
    assert isinstance(actual_jn, complex)
    assert cmath.isclose(actual_jn, expected_jn, rel_tol=TOL)

    # y_n complex
    expected_yn = spherical_yn(n, z)
    actual_yn = _spherical_yn_scalar_complex(n, z)
    assert isinstance(actual_yn, complex)
    assert cmath.isclose(actual_yn, expected_yn, rel_tol=TOL)

    # i_n complex
    expected_in = spherical_in(n, z)
    actual_in = _spherical_in_scalar_complex(n, z)
    assert isinstance(actual_in, complex)
    assert cmath.isclose(actual_in, expected_in, rel_tol=TOL)

    # k_n complex
    expected_kn = spherical_kn(n, z)
    actual_kn = _spherical_kn_scalar_complex(n, z)
    assert isinstance(actual_kn, complex)
    assert cmath.isclose(actual_kn, expected_kn, rel_tol=TOL)


@pytest.mark.parametrize('n', ORDERS)
@pytest.mark.parametrize('x', REAL_VALUES)
def test_spherical_derivatives_scalar(n, x):
    """Tests scalar derivative wrappers against SciPy."""
    
    # j_n'
    expected_jn_d = spherical_jn(n, x, derivative=True)
    actual_jn_d = _spherical_jn_d_scalar(n, x)
    assert isinstance(actual_jn_d, float)
    assert isclose(actual_jn_d, expected_jn_d, rel_tol=TOL)

    # y_n'
    expected_yn_d = spherical_yn(n, x, derivative=True)
    actual_yn_d = _spherical_yn_d_scalar(n, x)
    assert isinstance(actual_yn_d, float)
    assert isclose(actual_yn_d, expected_yn_d, rel_tol=TOL)

    # i_n'
    expected_in_d = spherical_in(n, x, derivative=True)
    actual_in_d = _spherical_in_d_scalar(n, x)
    assert isinstance(actual_in_d, float)
    assert isclose(actual_in_d, expected_in_d, rel_tol=TOL)

    # k_n'
    expected_kn_d = spherical_kn(n, x, derivative=True)
    actual_kn_d = _spherical_kn_d_scalar(n, x)
    assert isinstance(actual_kn_d, float)
    assert isclose(actual_kn_d, expected_kn_d, rel_tol=TOL)


@pytest.mark.parametrize('n', ORDERS)
def test_spherical_arrays_real(n):
    """Tests 1D real array wrappers (standard and derivative) against SciPy."""
    x_arr = np.array([0.5, 1.0, 2.5, 10.0], dtype=np.float64)

    # Standard evaluations
    np.testing.assert_allclose(spherical_jn_array(n, x_arr), spherical_jn(n, x_arr), rtol=TOL)
    np.testing.assert_allclose(spherical_yn_array(n, x_arr), spherical_yn(n, x_arr), rtol=TOL)
    np.testing.assert_allclose(spherical_in_array(n, x_arr), spherical_in(n, x_arr), rtol=TOL)
    np.testing.assert_allclose(spherical_kn_array(n, x_arr), spherical_kn(n, x_arr), rtol=TOL)

    # Derivative evaluations
    np.testing.assert_allclose(spherical_jn_array(n, x_arr, derivative=True), spherical_jn(n, x_arr, derivative=True), rtol=TOL)
    np.testing.assert_allclose(spherical_yn_array(n, x_arr, derivative=True), spherical_yn(n, x_arr, derivative=True), rtol=TOL)
    np.testing.assert_allclose(spherical_in_array(n, x_arr, derivative=True), spherical_in(n, x_arr, derivative=True), rtol=TOL)
    np.testing.assert_allclose(spherical_kn_array(n, x_arr, derivative=True), spherical_kn(n, x_arr, derivative=True), rtol=TOL)


@pytest.mark.parametrize('n', ORDERS)
def test_spherical_arrays_complex(n):
    """Tests 1D complex array wrappers against SciPy."""
    z_arr = np.array([0.5 + 0.5j, 1.0 - 2.0j, 3.0 + 4.0j], dtype=np.complex128)

    np.testing.assert_allclose(spherical_jn_complex_array(n, z_arr), spherical_jn(n, z_arr), rtol=TOL)
    np.testing.assert_allclose(spherical_yn_complex_array(n, z_arr), spherical_yn(n, z_arr), rtol=TOL)
    np.testing.assert_allclose(spherical_in_complex_array(n, z_arr), spherical_in(n, z_arr), rtol=TOL)
    np.testing.assert_allclose(spherical_kn_complex_array(n, z_arr), spherical_kn(n, z_arr), rtol=TOL)
