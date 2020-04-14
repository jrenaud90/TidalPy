import numpy as np

from scipy.constants import G

earth_radius = 6.37101e6
sun_mass = 1.988435e30
earth_semi_a = 1.4959789e11

def test_tidal_susceptibility():

    from TidalPy.tides.dissipation import calc_tidal_susceptibility, calc_tidal_susceptibility_reduced

    # Test float calculations - tidal suscept
    np.testing.assert_approx_equal(calc_tidal_susceptibility(1.5, 2., 3.),
                                   (3. / 2.) * 1.5**2 * 2.**5 / 3.**6)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, earth_radius, earth_semi_a),
                                   3.707e17)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(-sun_mass, earth_radius, earth_semi_a),
                                   3.707e17)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, earth_radius, -earth_semi_a),
                                   3.707e17)
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, -earth_radius, earth_semi_a),
                                   -3.707e17)

    # Test array calculations - tidal suscept
    np.testing.assert_approx_equal(calc_tidal_susceptibility(1.5, 2., np.linspace(1., 3., 5)),
                                   (3./2.)*1.5**2*2.**5/(np.linspace(1., 3., 5)**6))
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, earth_radius, earth_semi_a*np.linspace(0.9, 1.5, 4)),
                                   (3./2.)*sun_mass**2*earth_radius**5/((earth_semi_a*np.linspace(0.9, 1.5, 4))**6))
    np.testing.assert_approx_equal(calc_tidal_susceptibility(-sun_mass, earth_radius, earth_semi_a*np.linspace(0.9, 1.5, 4)),
                                   (3./2.)*sun_mass**2*earth_radius**5/((earth_semi_a*np.linspace(0.9, 1.5, 4))**6))
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, earth_radius, -earth_semi_a*np.linspace(0.9, 1.5, 4)),
                                   (3./2.)*sun_mass**2*earth_radius**5/((earth_semi_a*np.linspace(0.9, 1.5, 4))**6))
    np.testing.assert_approx_equal(calc_tidal_susceptibility(sun_mass, -earth_radius, earth_semi_a*np.linspace(0.9, 1.5, 4)),
                                   (3./2.)*sun_mass**2*(-earth_radius)**5/((earth_semi_a*np.linspace(0.9, 1.5, 4))**6))

    # Test float calculations - tidal suscept reduced
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(1.5, 2.) / G, (3. / 2.) * 1.5**2 * 2.**5)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, earth_radius), 4.178e84)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(-sun_mass, earth_radius), 4.178e84)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, -earth_radius), -4.178e84)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, 0.), 0.)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(0., earth_radius), 0.)
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, earth_radius),
                                   calc_tidal_susceptibility(sun_mass, earth_radius, 1.))
    np.testing.assert_approx_equal(calc_tidal_susceptibility_reduced(sun_mass, earth_radius),
                                   calc_tidal_susceptibility(sun_mass, earth_radius, earth_semi_a) * earth_semi_a**6)

def test_mode_collapse():

    from TidalPy.tides.dissipation import mode_collapse

    # Float testing
    spin_frequency, orbital_frequency, eccentricity, obliquity = \
        0.001, 0.002, 0.1, 0.2

    for max_order_l in [2, 5]:
        results_by_uniquefreq = mode_collapse(spin_frequency, orbital_frequency, eccentricity, obliquity,
                                use_obliquity=True, max_order_l=max_order_l)

        assert type(results_by_uniquefreq) is dict
        for freq_signature, (freq, result_by_orderl) in results_by_uniquefreq.items():
            assert type(freq_signature) is str
            assert type(freq) is float
            assert type(result_by_orderl) is dict

            # For each max_order_l result_by_orderl should have a length of max_order_l-1
            assert len(result_by_orderl) == max_order_l-1
            for order_l, (tidal_heating_terms, dUdM_terms, dUdw_terms, dUdO_terms) in result_by_orderl.items():
                assert type(order_l) is int
                assert type(tidal_heating_terms) is float
                assert type(dUdM_terms) is float
                assert type(dUdw_terms) is float
                assert type(dUdO_terms) is float

    # Array testing
    spin_frequency, orbital_frequency, eccentricity, obliquity = \
        np.linspace(0.001, 0.003, 4), np.linspace(0.002, 0.003, 4), 0.1, 0.2

    for max_order_l in [2, 5]:
        results_by_uniquefreq = mode_collapse(spin_frequency, orbital_frequency, eccentricity, obliquity,
                                              use_obliquity=True, max_order_l=max_order_l)

        assert type(results_by_uniquefreq) is dict
        for freq_signature, (freq, result_by_orderl) in results_by_uniquefreq.items():
            assert type(freq_signature) is str
            assert type(freq) is np.ndarray
            assert type(result_by_orderl) is dict

            # For each max_order_l result_by_orderl should have a length of max_order_l-1
            assert len(result_by_orderl) == max_order_l - 1
            for order_l, (tidal_heating_terms, dUdM_terms, dUdw_terms, dUdO_terms) in result_by_orderl.items():
                assert type(order_l) is int
                assert type(tidal_heating_terms) is np.ndarray
                assert type(dUdM_terms) is np.ndarray
                assert type(dUdw_terms) is np.ndarray
                assert type(dUdO_terms) is np.ndarray
            