""" Tests for `TidalPy.utilities.math.complex` """
import warnings
from math import isclose, isinf, isnan, nan, inf

import pytest
import numpy as np


from TidalPy.utilities.math.complex import hypot, csqrt, cexp, clog, cpow, cipow, cinv, cabs, cabs2


def compare_values(value_1, value_2):

    # Check for nans
    if isnan(value_1):
        assert isnan(value_2)
    elif isinf(value_1):
        assert isinf(value_2)
        assert np.sign(value_1) == np.sign(value_2)
    else:
        # Result should be finite; compare values up to a scale
        assert isclose(value_1, value_2)

def compare_values_complex(value_1, value_2):

    # Check real values
    compare_values(np.real(value_1), np.real(value_2))
    compare_values(np.imag(value_1), np.imag(value_2))


standard_list_complex = (
    0. + 0.j,
    2. + 0.j,
    5. + 0.j,
    # 0. + 5.0j,
    # 10.0 + 10.0j,
    # 10.0 - 10.0j,
    # -10.0 + 10.0j,
    # -10.0 - 10.0j,
    # 1.0e34 + 1.0j * 2.0e20,
    # 1.0e-10 + 1.0j * 2.e-10,
    # nan + 0.j,
    # 0. + nan * 1.j,
    # -inf + 0.j,
    # inf + 0.j,
    # 0 + -inf * 1.j,
    # 0 + inf * 1.j
    )

standard_list_float = (
    0.,
    2.,
    5.,
    -2.,
    1.0e34,
    1.0e-10,
    nan,
    -inf,
    inf
    )

@pytest.mark.parametrize('z', standard_list_complex)
def test_cabs(z):

    np_result = np.abs(z)
    tpy_result = cabs(z)

    compare_values_complex(np_result, tpy_result)

@pytest.mark.parametrize('z', standard_list_complex)
def test_cabs2(z):

    np_result = np.abs(z)**2
    tpy_result = cabs2(z)

    compare_values_complex(np_result, tpy_result)

@pytest.mark.parametrize('z', standard_list_complex)
def test_cinv(z):

    np_result = 1. / np.asarray(z, dtype=np.complex128)
    np_result = np_result
    tpy_result = cinv(z)

    try:
        compare_values_complex(np_result, tpy_result)
    except:
        breakpoint()

@pytest.mark.parametrize('a', standard_list_float)
@pytest.mark.parametrize('b', standard_list_float)
def test_hypot(a, b):

    np_result = np.hypot(a, b)
    tpy_result = hypot(a, b)

    compare_values(np_result, tpy_result)


@pytest.mark.parametrize('z', standard_list_complex)
def test_csqrt(z):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        np_result = np.sqrt(z)
    tpy_result = csqrt(z)

    compare_values_complex(np_result, tpy_result)


@pytest.mark.parametrize('z', standard_list_complex)
def test_cexp(z):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        np_result = np.exp(z)
    tpy_result = cexp(z)

    compare_values_complex(np_result, tpy_result)


@pytest.mark.parametrize('z', standard_list_complex)
def test_clog(z):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        np_result = np.log(z)
    tpy_result = clog(z)

    compare_values_complex(np_result, tpy_result)


@pytest.mark.parametrize('a', standard_list_complex)
@pytest.mark.parametrize('b', standard_list_complex)
def test_cpow(a, b):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        np_result = np.power(a, b)
    tpy_result = cpow(a, b)

    compare_values_complex(np_result, tpy_result)


@pytest.mark.parametrize('a', standard_list_complex)
@pytest.mark.parametrize('bi', (0, 1, 2, 3, -5, -200, 5, 50, 200))
def test_cipow(a, bi):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        np_result = np.power(a, bi)
    tpy_result = cipow(a, bi)

    compare_values_complex(np_result, tpy_result)
