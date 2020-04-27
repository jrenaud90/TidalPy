import pytest

import numpy as np
from scipy.constants import G

import TidalPy
TidalPy.use_disk = False
from TidalPy.exceptions import BadValueError

def test_m2Au():

    # Test Load
    from TidalPy.tools.conversions import m2Au, Au2m
    assert type(m2Au(2.)) is float
    assert type(Au2m(2.)) is float

    # Test Floats
    assert m2Au(1.) == 1.496e11**(-1)
    assert m2Au(1.496e11) == 1.
    assert m2Au(-1.) == -1.496e11**(-1)
    assert m2Au(-1.496e11) == -1.
    assert m2Au(0.) == 0.

    assert Au2m(1.) == 1.496e11
    assert Au2m(1.496e11**(-1)) == 1.
    assert Au2m(-1.) == -1.496e11
    assert Au2m(-1.496e11**(-1)) == -1.
    assert Au2m(0.) == 0.

    # Test Arrays
    test_array_au = np.linspace(0.9, 1.1, 4)
    test_array_m = 1.496e11 * test_array_au

    np.testing.assert_allclose(m2Au(test_array_m), test_array_au)
    np.testing.assert_allclose(m2Au(-test_array_m), -test_array_au)
    np.testing.assert_allclose(m2Au(np.zeros(4)), np.zeros(4))

    np.testing.assert_allclose(Au2m(test_array_au), test_array_m)
    np.testing.assert_allclose(Au2m(-test_array_au), -test_array_m)
    np.testing.assert_allclose(Au2m(np.zeros(4)), np.zeros(4))


def test_rads2days():

    # Test Load
    from TidalPy.tools.conversions import rads2days, days2rads
    assert type(rads2days(.2)) is float
    assert type(days2rads(.2)) is float

    # Test Floats
    assert days2rads(1.) == 2. * np.pi / 86400.
    assert days2rads(2. * np.pi / 86400.) == 1.
    assert days2rads(-1.) == -2. * np.pi / 86400.
    assert days2rads(-2. * np.pi / 86400.) == -1.
    with pytest.raises(ZeroDivisionError):
        assert days2rads(0.)

    assert rads2days(1.) == 86400. / (2. * np.pi)
    assert rads2days(86400. / (2. * np.pi)) == 1.
    assert rads2days(-1.) == -86400. / (2. * np.pi)
    assert rads2days(-86400. / (2. * np.pi)) == -1.
    with pytest.raises(ZeroDivisionError):
        assert rads2days(0.)

    # Test Arrays
    day_array = np.linspace(0.9, 1.1, 5)
    rads_array = 2. * np.pi / (day_array * 86400.)

    np.testing.assert_allclose(days2rads(day_array), rads_array)
    np.testing.assert_allclose(days2rads(-day_array), -rads_array)
    np.testing.assert_allclose(rads2days(rads_array), day_array)
    np.testing.assert_allclose(rads2days(-rads_array), -day_array)
    with pytest.raises(ZeroDivisionError):
        assert rads2days(np.zeros(5))
        assert days2rads(np.zeros(5))

def test_sec2myr():

    # Test Load
    from TidalPy.tools.conversions import sec2myr, myr2sec
    assert type(sec2myr(2.)) is float
    assert type(myr2sec(2.)) is float

    # Test Floats
    assert sec2myr(1.) == 1. / 3.154e13
    assert sec2myr(3.154e13) == 1.
    assert sec2myr(-1.) == -1. / 3.154e13
    assert sec2myr(-3.154e13) == -3.154e13
    assert sec2myr(0.) == 0.

    assert myr2sec(1.) == 3.154e13
    assert myr2sec(1. / 3.154e13) == 1.
    assert myr2sec(-1.) == -3.154e13
    assert myr2sec(-1. / 3.154e13) == -1.
    assert myr2sec(0.) == 0.

    # Test Arrays
    myr_array = np.linspace(0.9, 1.1, 5)
    sec_array = 3.154e13 * myr_array

    np.testing.assert_allclose(sec2myr(sec_array), myr_array)
    np.testing.assert_allclose(sec2myr(-sec_array), -myr_array)
    np.testing.assert_allclose(sec2myr(np.zeros(5)), np.zeros(5))
    np.testing.assert_allclose(myr2sec(myr_array), sec_array)
    np.testing.assert_allclose(myr2sec(-myr_array), -sec_array)
    np.testing.assert_allclose(myr2sec(np.zeros(5)), np.zeros(5))

def test_orbital_motion2semi_a():

    # Test Import
    from TidalPy.tools.conversions import orbital_motion2semi_a, semi_a2orbital_motion
    assert type(orbital_motion2semi_a(1., 2., 3.)) is float
    assert type(semi_a2orbital_motion(1., 2., 3.)) is float

    # Test Floats
    assert orbital_motion2semi_a(1., 1. / G, 0.) == 1.
    assert orbital_motion2semi_a(-1., 1. / G, 0.) == 1.
    assert orbital_motion2semi_a(0.2, 10., 5.) == (G * 15. / 0.2**2)**(1 / 3)
    assert orbital_motion2semi_a(-0.2, 10., 5.) == (G * 15. / 0.2**2)**(1 / 3)

    with pytest.raises(ZeroDivisionError):
        assert orbital_motion2semi_a(0., 1., 1.)
    with pytest.raises(BadValueError):
        assert orbital_motion2semi_a(1., -1., 1.)
    with pytest.raises(BadValueError):
        assert orbital_motion2semi_a(1., 0., 1.)
    with pytest.raises(BadValueError):
        assert orbital_motion2semi_a(1., 1., -1.)

    assert semi_a2orbital_motion(1., 1./G, 0.) == 1.
    assert semi_a2orbital_motion(-1., 1./G, 0.) == 1.
    assert semi_a2orbital_motion(100., 10., 5.) == (G * 15. / 100.**3)**(1 / 2)
    assert semi_a2orbital_motion(-100., 10., 5.) == (G * 15. / 100.**3)**(1 / 2)

    with pytest.raises(ZeroDivisionError):
        assert semi_a2orbital_motion(0., 1., 1.)
    with pytest.raises(BadValueError):
        assert semi_a2orbital_motion(1., -1., 1.)
    with pytest.raises(BadValueError):
        assert semi_a2orbital_motion(1., 0., 1.)
    with pytest.raises(BadValueError):
        assert semi_a2orbital_motion(1., 1., -1.)
