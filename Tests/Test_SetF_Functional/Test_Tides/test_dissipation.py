
import pytest

import numpy as np
from scipy.constants import G

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
    assert type(reduced_float_result) is float
    assert type(full_float_result) is float

    # Test Arrays
    full_float_result = calc_tidal_susceptibility(1., 2., np.linspace(1., 5., 5))
    assert type(full_float_result) is np.ndarray


def test_tidal_susceptibility():

    from TidalPy.tides.dissipation import calc_tidal_susceptibility, calc_tidal_susceptibility_reduced

    # Test float calculations - tidal susceptibility & reduced
    np.testing.assert_approx_equal(calc_tidal_susceptibility(1.5, 2., 3.) / G,
                                   (3. / 2.) * 1.5**2 * 2.**5 / 3.**6)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, earth_radius),
                                   expected_susceptibility_reduced)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(-sun_mass, earth_radius),
                                   expected_susceptibility_reduced)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, -earth_radius),
                                   -expected_susceptibility_reduced)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, 0.), 0.)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(0., earth_radius), 0.)

    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(1.5, 2.)/G,
                                   (3. / 2.) * 1.5**2 * 2.**5)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, earth_radius, earth_semi_a),
                                   expected_susceptibility)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(-sun_mass, earth_radius, earth_semi_a),
                                   expected_susceptibility)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, earth_radius, -earth_semi_a),
                                   expected_susceptibility)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, -earth_radius, earth_semi_a),
                                   -expected_susceptibility)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, 0., earth_semi_a), 0.)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(0., earth_radius, earth_semi_a), 0.)
    with pytest.raises(ZeroDivisionError):
        assert calc_tidal_susceptibility(sun_mass, earth_radius, 0)

    # Test array calculations - tidal susceptibility & reduced
    earth_semi_a_array = earth_semi_a * np.linspace(0.9, 1.5, 4)
    expected_susceptibility_array = \
        (3. / 2.) * G * sun_mass**2 * earth_radius**5 / earth_semi_a_array**6

    np.testing.assert_allclose(calc_tidal_susceptibility(1.5, 2., np.linspace(1., 3., 4)) / G,
                                   (3. / 2.) * 1.5**2 * 2.**5 / np.linspace(1., 3., 4)**6)
    np.testing.assert_allclose(calc_tidal_susceptibility(sun_mass, earth_radius, earth_semi_a_array),
                               expected_susceptibility_array)
    np.testing.assert_allclose(calc_tidal_susceptibility(-sun_mass, earth_radius, earth_semi_a_array),
                               expected_susceptibility_array)
    np.testing.assert_allclose(calc_tidal_susceptibility(sun_mass, earth_radius, -earth_semi_a_array),
                               expected_susceptibility_array)
    np.testing.assert_allclose(calc_tidal_susceptibility(sun_mass, -earth_radius, earth_semi_a_array),
                               -expected_susceptibility_array)


def _test_calculate_modes(order_l, eccen_trunc_lvl):

    from TidalPy.tides.dissipation import calculate_terms

    # Float testing
    spin_frequency, orbital_frequency = 0.001, 0.002
    semi_major_axis = 100.
    eccentricity, obliquity = 0.1, 0.2
    radius = 10.

    # Array testing
    spin_frequency_array = np.linspace(0.001, 0.005, 5)
    orbital_frequency_array = 0.5 * spin_frequency_array
    semi_major_axis_array = 100. * spin_frequency_array
    eccentricity_array = np.linspace(0.1, 0.5, 6)
    obliquity_array = np.linspace(0.1, 0.5, 6)

    # Matrix Testing
    spin_frequency_mtx, eccentricity_mtx = np.meshgrid(spin_frequency_array, eccentricity_array)
    orbital_frequency_mtx = 0.5 * spin_frequency_mtx
    semi_major_axis_mtx = 100. * spin_frequency_mtx

    # Test pure float version
    unique_freqs_float, results_by_uniquefreq_float = \
        calculate_terms(orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, radius,
                        use_obliquity=True, eccentricity_truncation_lvl=eccen_trunc_lvl,
                        max_order_l=order_l)

    for freq_signature, freq in unique_freqs_float.items():
        assert type(freq_signature) is tuple
        assert len(freq_signature) == 2
        assert type(freq_signature[0]) is int
        assert type(freq_signature[1]) is int
        assert type(freq) is float

    assert type(results_by_uniquefreq_float) is dict
    for _, results_by_orderl in results_by_uniquefreq_float.items():
        assert type(results_by_orderl) is dict

        for order_l, (tidal_heating_terms, dUdM_terms, dUdw_terms, dUdO_terms) in results_by_orderl.items():
            assert type(order_l) in [int, np.int, np.int32]
            assert type(tidal_heating_terms) in [float, np.float, np.float64]
            assert type(dUdM_terms) in [float, np.float, np.float64]
            assert type(dUdw_terms) in [float, np.float, np.float64]
            assert type(dUdO_terms) in [float, np.float, np.float64]


    # Test array version using array'd frequencies
    unique_freqs_freq_array, results_by_uniquefreq_freq_array = \
        calculate_terms(orbital_frequency_array, spin_frequency_array, eccentricity, obliquity, semi_major_axis_array,
                        radius, use_obliquity=True, eccentricity_truncation_lvl=eccen_trunc_lvl, max_order_l=order_l)

    for freq_signature, freq in unique_freqs_freq_array.items():
        assert type(freq_signature) is tuple
        assert len(freq_signature) == 2
        assert type(freq_signature[0]) is int
        assert type(freq_signature[1]) is int
        assert type(freq) is np.ndarray

    assert type(results_by_uniquefreq_freq_array) is dict
    for _, results_by_orderl in results_by_uniquefreq_freq_array.items():
        assert type(results_by_orderl) is dict

        for order_l, (tidal_heating_terms, dUdM_terms, dUdw_terms, dUdO_terms) in results_by_orderl.items():
            assert type(order_l) in [int, np.int, np.int32]
            assert type(tidal_heating_terms) is np.ndarray
            assert type(dUdM_terms) is np.ndarray
            assert type(dUdw_terms) is np.ndarray
            assert type(dUdO_terms) is np.ndarray


    # Test array version using array'd eccentricity
    unique_freqs_eccen_array, results_by_uniquefreq_eccen_array = \
        calculate_terms(orbital_frequency, spin_frequency, eccentricity_array, obliquity, semi_major_axis,
                        radius, use_obliquity=True, eccentricity_truncation_lvl=eccen_trunc_lvl, max_order_l=order_l)

    for freq_signature, freq in unique_freqs_eccen_array.items():
        assert type(freq_signature) is tuple
        assert len(freq_signature) == 2
        assert type(freq_signature[0]) is int
        assert type(freq_signature[1]) is int
        assert type(freq) is float

    assert type(results_by_uniquefreq_eccen_array) is dict
    for _, results_by_orderl in results_by_uniquefreq_eccen_array.items():
        assert type(results_by_orderl) is dict

        for order_l, (tidal_heating_terms, dUdM_terms, dUdw_terms, dUdO_terms) in results_by_orderl.items():
            assert type(order_l) in [int, np.int, np.int32]
            assert type(tidal_heating_terms) is np.ndarray
            assert type(dUdM_terms) is np.ndarray
            assert type(dUdw_terms) is np.ndarray
            assert type(dUdO_terms) is np.ndarray


    # Test array version using array'd obliquity
    unique_freqs_obliq_array, results_by_uniquefreq_obliq_array = \
        calculate_terms(orbital_frequency, spin_frequency, eccentricity, obliquity_array, semi_major_axis,
                        radius, use_obliquity=True, eccentricity_truncation_lvl=eccen_trunc_lvl, max_order_l=order_l)

    for freq_signature, freq in unique_freqs_obliq_array.items():
        assert type(freq_signature) is tuple
        assert len(freq_signature) == 2
        assert type(freq_signature[0]) is int
        assert type(freq_signature[1]) is int
        assert type(freq) is float

    assert type(results_by_uniquefreq_obliq_array) is dict
    for _, results_by_orderl in results_by_uniquefreq_obliq_array.items():
        assert type(results_by_orderl) is dict

        for order_l, (tidal_heating_terms, dUdM_terms, dUdw_terms, dUdO_terms) in results_by_orderl.items():
            assert type(order_l) in [int, np.int, np.int32]
            assert type(tidal_heating_terms) is np.ndarray
            assert type(dUdM_terms) is np.ndarray
            assert type(dUdw_terms) is np.ndarray
            assert type(dUdO_terms) is np.ndarray


    # Test matrix version using meshgrid of frequency and eccentricity
    unique_freqs_mtx, results_by_uniquefreq_mtx = \
        calculate_terms(orbital_frequency_mtx, spin_frequency_mtx, eccentricity_mtx, obliquity, semi_major_axis_mtx,
                        radius, use_obliquity=True, eccentricity_truncation_lvl=eccen_trunc_lvl, max_order_l=order_l)

    for freq_signature, freq in unique_freqs_mtx.items():
        assert type(freq_signature) is tuple
        assert len(freq_signature) == 2
        assert type(freq_signature[0]) is int
        assert type(freq_signature[1]) is int
        assert type(freq) is np.ndarray
        assert len(freq.shape) == 2

    assert type(results_by_uniquefreq_mtx) is dict
    for _, results_by_orderl in results_by_uniquefreq_mtx.items():
        assert type(results_by_orderl) is dict

        for order_l, (tidal_heating_terms, dUdM_terms, dUdw_terms, dUdO_terms) in results_by_orderl.items():
            assert type(order_l) in [int, np.int, np.int32]
            assert type(tidal_heating_terms) is np.ndarray
            assert type(dUdM_terms) is np.ndarray
            assert type(dUdw_terms) is np.ndarray
            assert type(dUdO_terms) is np.ndarray

            assert len(tidal_heating_terms.shape) == 2
            assert len(dUdM_terms.shape) == 2
            assert len(dUdw_terms.shape) == 2
            assert len(dUdO_terms.shape) == 2

    return True


def test_l2e2():
    assert _test_calculate_modes(2, 2)

def test_l2e8():
    assert _test_calculate_modes(2, 8)

def test_l2e10():
    assert _test_calculate_modes(2, 10)

def test_l2e20():
    assert _test_calculate_modes(2, 20)

def test_l4e2():
    assert _test_calculate_modes(4, 2)

def test_l4e10():
    assert _test_calculate_modes(4, 10)

def test_l6e2():
    assert _test_calculate_modes(6, 2)

def test_l6e8():
    assert _test_calculate_modes(6, 8)

def test_l6e10():
    assert _test_calculate_modes(6, 10)

# FIXME: once l6 and l7 e20 are done
# def test_l6e20():
#     assert _test_calculate_modes(6, 20)