import math

import numpy as np
import pytest

import TidalPy
from TidalPy.constants import G
from TidalPy.exceptions import BadValueError
from TidalPy.utilities.conversions import (Au2m, days2rads, m2Au, myr2sec, orbital_motion2semi_a, rads2days, sec2myr,
                                               semi_a2orbital_motion)

def test_sec2myr():
    # Build arrays for testing
    zero_array = np.zeros(10, dtype=np.float64)
    one_array = np.ones(10, dtype=np.float64)
    sec_array = 3.154e13 * np.ones(10, dtype=np.float64)

    # Test sec2myr - Floats
    assert sec2myr(0.) == 0.
    assert sec2myr(3.154e13) == 1.
    assert sec2myr(-3.154e13) == -1.

    # Test sec2myr - Arrays
    np.testing.assert_array_almost_equal(sec2myr(zero_array), zero_array)
    np.testing.assert_array_almost_equal(sec2myr(sec_array), one_array)
    np.testing.assert_array_almost_equal(sec2myr(-sec_array), -one_array)

    # Test myr2sec - Floats
    assert myr2sec(0.) == 0.
    assert myr2sec(1.) == 3.154e13
    assert myr2sec(-1.) == -3.154e13

    # Test myr2sec - Arrays
    np.testing.assert_array_almost_equal(myr2sec(zero_array), zero_array)
    np.testing.assert_array_almost_equal(myr2sec(one_array), sec_array)
    np.testing.assert_array_almost_equal(myr2sec(-one_array), -sec_array)


def test_Au2m():
    # Build arrays for testing
    zero_array = np.zeros(10, dtype=np.float64)
    one_array = np.ones(10, dtype=np.float64)
    meter_array = 1.496e11 * np.ones(10, dtype=np.float64)

    # Test m2Au - Floats
    assert m2Au(0.) == 0.
    assert m2Au(1.496e11) == 1.
    assert m2Au(-1.496e11) == -1.

    # Test m2Au - Arrays
    np.testing.assert_array_almost_equal(m2Au(zero_array), zero_array)
    np.testing.assert_array_almost_equal(m2Au(meter_array), one_array)
    np.testing.assert_array_almost_equal(m2Au(-meter_array), -one_array)

    # Test Au2m - Floats
    assert Au2m(0.) == 0.
    assert Au2m(1.) == 1.496e11
    assert Au2m(-1.) == -1.496e11

    # Test myr2sec - Arrays
    np.testing.assert_array_almost_equal(Au2m(zero_array), zero_array)
    np.testing.assert_array_almost_equal(Au2m(one_array), meter_array)
    np.testing.assert_array_almost_equal(Au2m(-one_array), -meter_array)


def test_rads2days():
    # Build arrays for testing
    day_to_radians = 2. * np.pi / 86400.
    zero_array = np.zeros(10, dtype=np.float64)
    inf_array = np.inf * np.ones(10, dtype=np.float64)
    one_array = np.ones(10, dtype=np.float64)
    radian_array = day_to_radians * np.ones(10, dtype=np.float64)

    # Test rads2days - Floats
    with pytest.raises(ZeroDivisionError) as e_info:
        _ = rads2days(0.)
    assert rads2days(day_to_radians) == 1.
    assert rads2days(-day_to_radians) == -1.

    # Test rads2days - Arrays
    np.testing.assert_array_almost_equal(rads2days(zero_array), inf_array)
    np.testing.assert_array_almost_equal(rads2days(radian_array), one_array)
    np.testing.assert_array_almost_equal(rads2days(-radian_array), -one_array)

    # Test days2rads - Floats
    with pytest.raises(ZeroDivisionError) as e_info:
        _ = days2rads(0.)
    assert days2rads(1.) == day_to_radians
    assert days2rads(-1.) == -day_to_radians

    # Test days2rads - Arrays
    np.testing.assert_array_almost_equal(days2rads(zero_array), inf_array)
    np.testing.assert_array_almost_equal(days2rads(one_array), radian_array)
    np.testing.assert_array_almost_equal(days2rads(-one_array), -radian_array)


def test_semi_a2orbital_motion():
    # Build Arrays for testing
    earth_distance = 1.49597887e11  # meters
    earth_mass = 5.972e24  # kg
    sun_mass = 1.988435e30  # kg
    earth_orb_motion = (G * (earth_mass + sun_mass) / earth_distance**3)**(1 / 2)  # radians second-1

    zero_array = np.zeros(10, dtype=np.float64)
    distance_array = earth_distance * np.ones(10, dtype=np.float64)
    frequency_array = earth_orb_motion * np.ones(10, dtype=np.float64)

    # Test semi_a2orbital_motion - Floats
    with pytest.raises(ZeroDivisionError) as e_info:
        _ = semi_a2orbital_motion(0., sun_mass, earth_mass)
    np.testing.assert_almost_equal(semi_a2orbital_motion(earth_distance, sun_mass, earth_mass), earth_orb_motion)
    assert math.isnan(semi_a2orbital_motion(-earth_distance, sun_mass, earth_mass))

    # Test semi_a2orbital_motion - Arrays
    assert np.all(np.isinf(semi_a2orbital_motion(zero_array, sun_mass, earth_mass)))
    np.testing.assert_array_almost_equal(
        semi_a2orbital_motion(distance_array, sun_mass, earth_mass), frequency_array, decimal=9
        )
    assert np.all(np.isnan(semi_a2orbital_motion(-1. * distance_array, sun_mass, earth_mass)))

    # Test orbital_motion2semi_a - Floats
    with pytest.raises(ZeroDivisionError) as e_info:
        _ = orbital_motion2semi_a(0., sun_mass, earth_mass)
    np.testing.assert_almost_equal(
        round(orbital_motion2semi_a(earth_orb_motion, sun_mass, earth_mass) * 1.e3) / 1.e3, earth_distance
        )
    np.testing.assert_almost_equal(
        round(orbital_motion2semi_a(-earth_orb_motion, sun_mass, earth_mass) * 1.e3) / 1.e3, earth_distance
        )

    # Test orbital_motion2semi_a - Arrays
    assert np.all(np.isinf(orbital_motion2semi_a(zero_array, sun_mass, earth_mass)))
    np.testing.assert_array_almost_equal(
        orbital_motion2semi_a(frequency_array, sun_mass, earth_mass), distance_array, decimal=3
        )
    np.testing.assert_array_almost_equal(
        orbital_motion2semi_a(-1. * frequency_array, sun_mass, earth_mass), distance_array, decimal=3
        )

    # Test negative masses
    with pytest.raises(BadValueError) as e_info:
        _ = orbital_motion2semi_a(1., -1. * sun_mass, earth_mass)
        _ = orbital_motion2semi_a(1., 0. * sun_mass, earth_mass)
        _ = orbital_motion2semi_a(1., sun_mass, -1. * earth_mass)

    with pytest.raises(BadValueError) as e_info:
        _ = semi_a2orbital_motion(1., -1. * sun_mass, earth_mass)
        _ = semi_a2orbital_motion(1., 0. * sun_mass, earth_mass)
        _ = semi_a2orbital_motion(1., sun_mass, -1. * earth_mass)

    # Make sure that zero target masses are passed okay
    assert orbital_motion2semi_a(G**(1 / 2), 1., 0.) == 1.
    assert orbital_motion2semi_a(G**(1 / 2), .5, .5) == 1.

    #  The below actually fails normally. Thus I added the rounding.
    assert round(1.e12 * semi_a2orbital_motion(G**(1 / 3), 1., 0.)) / 1.e12 == 1.
    assert round(1.e12 * semi_a2orbital_motion(G**(1 / 3), .5, .5)) / 1.e12 == 1.
