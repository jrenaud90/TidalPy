import numpy as np
import pytest

import TidalPy
from TidalPy.constants import G


earth_radius = 6.37101e6
sun_mass = 1.988435e30
earth_semi_a = 1.4959789e11
expected_susceptibility_reduced = 1.5 * G * earth_radius**5 * sun_mass**2
expected_susceptibility = expected_susceptibility_reduced / earth_semi_a**6


def test_tidal_susceptibility_load_and_types():
    from TidalPy.tides.dissipation import calc_tidal_susceptibility, calc_tidal_susceptibility_reduced

    # Test floats
    reduced_float_result = calc_tidal_susceptibility_reduced(1., 2.)
    full_float_result = calc_tidal_susceptibility(1., 2., 3.)
    assert type(reduced_float_result) in [float, np.float64]
    assert type(full_float_result) in [float, np.float64]

    # Test Arrays
    full_float_result = calc_tidal_susceptibility(1., 2., np.linspace(1., 5., 5))
    assert type(full_float_result) is np.ndarray


def test_tidal_susceptibility():
    from TidalPy.tides.dissipation import calc_tidal_susceptibility, calc_tidal_susceptibility_reduced

    # Test float calculations - tidal susceptibility & reduced
    np.testing.assert_approx_equal(
        calc_tidal_susceptibility(1.5, 2., 3.) / G,
        (3. / 2.) * 1.5**2 * 2.**5 / 3.**6
        )
    np.testing.assert_approx_equal(
        calc_tidal_susceptibility_reduced(sun_mass, earth_radius),
        expected_susceptibility_reduced
        )
    np.testing.assert_approx_equal(
        calc_tidal_susceptibility_reduced(-sun_mass, earth_radius),
        expected_susceptibility_reduced
        )
    np.testing.assert_approx_equal(
        calc_tidal_susceptibility_reduced(sun_mass, -earth_radius),
        -expected_susceptibility_reduced
        )
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, 0.), 0.)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(0., earth_radius), 0.)

    np.testing.assert_approx_equal(
        calc_tidal_susceptibility_reduced(1.5, 2.) / G,
        (3. / 2.) * 1.5**2 * 2.**5
        )
    np.testing.assert_approx_equal(
        calc_tidal_susceptibility(sun_mass, earth_radius, earth_semi_a),
        expected_susceptibility
        )
    np.testing.assert_approx_equal(
        calc_tidal_susceptibility(-sun_mass, earth_radius, earth_semi_a),
        expected_susceptibility
        )
    np.testing.assert_approx_equal(
        calc_tidal_susceptibility(sun_mass, earth_radius, -earth_semi_a),
        expected_susceptibility
        )
    np.testing.assert_approx_equal(
        calc_tidal_susceptibility(sun_mass, -earth_radius, earth_semi_a),
        -expected_susceptibility
        )
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, 0., earth_semi_a), 0.)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(0., earth_radius, earth_semi_a), 0.)
    with pytest.raises(ZeroDivisionError):
        assert calc_tidal_susceptibility(sun_mass, earth_radius, 0)

    # Test array calculations - tidal susceptibility & reduced
    earth_semi_a_array = earth_semi_a * np.linspace(0.9, 1.5, 4)
    expected_susceptibility_array = \
        (3. / 2.) * G * sun_mass**2 * earth_radius**5 / earth_semi_a_array**6

    np.testing.assert_allclose(
        calc_tidal_susceptibility(1.5, 2., np.linspace(1., 3., 4)) / G,
        (3. / 2.) * 1.5**2 * 2.**5 / np.linspace(1., 3., 4)**6
        )
    np.testing.assert_allclose(
        calc_tidal_susceptibility(sun_mass, earth_radius, earth_semi_a_array),
        expected_susceptibility_array
        )
    np.testing.assert_allclose(
        calc_tidal_susceptibility(-sun_mass, earth_radius, earth_semi_a_array),
        expected_susceptibility_array
        )
    np.testing.assert_allclose(
        calc_tidal_susceptibility(sun_mass, earth_radius, -earth_semi_a_array),
        expected_susceptibility_array
        )
    np.testing.assert_allclose(
        calc_tidal_susceptibility(sun_mass, -earth_radius, earth_semi_a_array),
        -expected_susceptibility_array
        )
