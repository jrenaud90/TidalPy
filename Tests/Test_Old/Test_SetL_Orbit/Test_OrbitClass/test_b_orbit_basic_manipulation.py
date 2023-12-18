import numpy as np

import TidalPy


from TidalPy.structures import build_world, scale_from_world
from TidalPy.structures.orbit import OrbitBase
from TidalPy.utilities.conversions import days2rads, orbital_motion2semi_a

io_config = {
    "name"           : "Io",
    "type"           : "layered",
    "radius"         : 1821.49e3,
    "orbital_period" : 1.769,
    "eccentricity"   : 0.0041,
    "spin_period"    : 1.769,
    "albedo"         : 0.63,
    "force_spin_sync": True,
    "layers"         : {
        "Core"  : {
            "type"    : "iron",
            "is_tidal": False,
            "radius"  : 810.0e3,
            "density" : 5200.
            },
        "Mantle": {
            "type"               : "rock",
            "is_tidal"           : True,
            "radius"             : 1821.49e3,
            "surface_temperature": 100.0,
            "density"            : 3200.
            }
        }
    }


def test_get_orbital_parameters():
    """ Make sure that orbital parameters set in a world_types configuration are accessible as expected. """

    # Build system
    star = build_world('55Cnc')
    io = build_world('io', world_config=io_config)
    small_io = scale_from_world(io, new_name='small_io', radius_scale=0.5)
    orbit = OrbitBase(star, io, small_io)

    # Check that the eccentricity and orbital period were loaded in correctly
    np.testing.assert_approx_equal(orbit.get_eccentricity(small_io), io_config['eccentricity'])
    np.testing.assert_approx_equal(small_io.eccentricity, io_config['eccentricity'])
    np.testing.assert_approx_equal(orbit.get_orbital_period(small_io), io_config['orbital_period'])
    np.testing.assert_approx_equal(small_io.orbital_period, io_config['orbital_period'])
    expected_orb_freq = days2rads(io_config['orbital_period'])
    np.testing.assert_approx_equal(orbit.get_orbital_frequency(small_io), expected_orb_freq)
    np.testing.assert_approx_equal(small_io.orbital_frequency, expected_orb_freq)
    expected_semi_a = orbital_motion2semi_a(expected_orb_freq, io.mass, small_io.mass)
    np.testing.assert_approx_equal(orbit.get_semi_major_axis(small_io), expected_semi_a)
    np.testing.assert_approx_equal(small_io.semi_major_axis, expected_semi_a)


def test_set_orbital_parameters():
    """ Make sure that orbital parameters can be set as expected. """

    # Build system
    star = build_world('55Cnc')
    io = build_world('io', world_config=io_config)
    small_io = scale_from_world(io, new_name='small_io', radius_scale=0.5)
    orbit = OrbitBase(star, io, small_io)

    # Test setting specific parameter as float
    orbit.set_eccentricity(small_io, 0.5)
    assert orbit.get_eccentricity(small_io) == 0.5
    assert small_io.eccentricity == 0.5

    # Test setting specific parameter as array
    eccentricities = np.linspace(0.1, 0.5, 10)
    orbit.set_eccentricity(small_io, eccentricities)
    assert orbit.get_eccentricity(small_io) is eccentricities
    assert small_io.eccentricity is eccentricities

    # Test setting multiple items as floats
    eccentricity = 0.01
    orbital_period = 10.
    expected_orbtial_freq = days2rads(orbital_period)
    expected_semi_a = orbital_motion2semi_a(expected_orbtial_freq, io.mass, small_io.mass)
    orbit.set_state(small_io, eccentricity=eccentricity, orbital_period=orbital_period)
    assert orbit.get_eccentricity(small_io) == eccentricity
    assert small_io.eccentricity == eccentricity
    assert orbit.get_orbital_period(small_io) == orbital_period
    assert small_io.orbital_period == orbital_period
    np.testing.assert_approx_equal(orbit.get_orbital_frequency(small_io), expected_orbtial_freq)
    np.testing.assert_approx_equal(small_io.orbital_frequency, expected_orbtial_freq)
    np.testing.assert_approx_equal(orbit.get_semi_major_axis(small_io), expected_semi_a)
    np.testing.assert_approx_equal(small_io.semi_major_axis, expected_semi_a)

    # Test setting multiple items as arrays
    eccentricity = np.linspace(0.1, 0.5, 10)
    orbital_period = np.linspace(5., 10., 10)
    expected_orbtial_freq = days2rads(orbital_period)
    expected_semi_a = orbital_motion2semi_a(expected_orbtial_freq, io.mass, small_io.mass)
    orbit.set_state(small_io, eccentricity=eccentricity, orbital_period=orbital_period)
    assert orbit.get_eccentricity(small_io) is eccentricity
    assert small_io.eccentricity is eccentricity
    assert orbit.get_orbital_period(small_io) is orbital_period
    assert small_io.orbital_period is orbital_period
    np.testing.assert_allclose(orbit.get_orbital_frequency(small_io), expected_orbtial_freq)
    np.testing.assert_allclose(small_io.orbital_frequency, expected_orbtial_freq)
    np.testing.assert_allclose(orbit.get_semi_major_axis(small_io), expected_semi_a)
    np.testing.assert_allclose(small_io.semi_major_axis, expected_semi_a)


def test_get_tidal_host_orbital_parameters_host_is_star():
    """ This will test how an orbit class handles getting orbital parameters for a tidal host.

    A tidal host's orbital parameters are stored in one of the tidal bodies, the one raising the tides on the host.

    For this test, the tidal host is the central star.
    """

    star = build_world('55Cnc')
    io = build_world('io', world_config=io_config)
    # Build some additional world_types to use as tidal bodies
    small_io = scale_from_world(io, new_name='small_io', radius_scale=0.5)
    tiny_io = scale_from_world(io, new_name='tiny_io', radius_scale=0.1)

    # Construct orbit using all of the tidal bodies
    orbit = OrbitBase(star=star, tidal_host=star, tidal_bodies=[io, small_io, tiny_io], host_tide_raiser=small_io)

    # Make sure tidal host was stored correctly.
    assert orbit.tidal_host is star

    # Check that the constructor set the correct tide raiser.
    assert orbit.host_tide_raiser is small_io

    # Check the host tide raiser setter is working correctly.
    orbit.set_host_tide_raiser(tiny_io)
    assert orbit.host_tide_raiser is tiny_io

    # Set an orbital parameter for the tide raiser and see if the tidal host picks up the change - float
    orbit.set_eccentricity(tiny_io, 0.5)
    assert orbit.get_eccentricity(star) == 0.5
    assert star.eccentricity == 0.5

    # Set an orbital parameter for the tide raiser and see if the tidal host picks up the change - array
    eccentricities = np.linspace(0.1, 0.5, 10)
    orbit.set_eccentricity(tiny_io, eccentricities)
    assert orbit.get_eccentricity(star) is eccentricities
    assert star.eccentricity is eccentricities


def test_get_tidal_host_orbital_parameters_host_is_world():
    """ This will test how an orbit class handles getting orbital parameters for a tidal host.

    A tidal host's orbital parameters are stored in one of the tidal bodies, the one raising the tides on the host.

    For this test, the tidal host is a LayeredWorld.
    """

    star = build_world('55Cnc')
    io = build_world('io', world_config=io_config)
    # Build some additional world_types to use as tidal bodies
    small_io = scale_from_world(io, new_name='small_io', radius_scale=0.5)
    tiny_io = scale_from_world(io, new_name='tiny_io', radius_scale=0.1)
    micro_io = scale_from_world(io, new_name='micro_io', radius_scale=0.1)

    # Construct orbit using all of the tidal bodies
    orbit = OrbitBase(star=star, tidal_host=io, tidal_bodies=[small_io, tiny_io, micro_io], host_tide_raiser=small_io)

    # Make sure tidal host was stored correctly.
    assert orbit.tidal_host is io

    # Check that the constructor set the correct tide raiser.
    assert orbit.host_tide_raiser is small_io

    # Check the host tide raiser setter is working correctly.
    orbit.set_host_tide_raiser(micro_io)
    assert orbit.host_tide_raiser is micro_io

    # Set an orbital parameter for the tide raiser and see if the tidal host picks up the change - float
    orbit.set_eccentricity(micro_io, 0.5)
    assert orbit.get_eccentricity(io) == 0.5
    assert io.eccentricity == 0.5

    # Set an orbital parameter for the tide raiser and see if the tidal host picks up the change - array
    eccentricities = np.linspace(0.1, 0.5, 10)
    orbit.set_eccentricity(micro_io, eccentricities)
    assert orbit.get_eccentricity(io) is eccentricities
    assert io.eccentricity is eccentricities

    # Check that the stellar distances make sense
    orbit.set_orbital_period(io, 100., set_stellar_orbit=True)
    orbit.set_orbital_period(small_io, 1., set_stellar_orbit=True)
    orbit.set_orbital_period(tiny_io, 2., set_stellar_orbit=True)
    orbit.set_orbital_period(micro_io, 3., set_stellar_orbit=True)
    #    Stellar distances should match the tidal host's
    assert orbit.get_stellar_distance(small_io) == orbit.get_stellar_distance(io)
    assert orbit.get_stellar_distance(tiny_io) == orbit.get_stellar_distance(io)
    assert orbit.get_stellar_distance(micro_io) == orbit.get_stellar_distance(io)
    #    They should not match the tidal world's semi-major axis
    assert orbit.get_stellar_distance(small_io) != orbit.get_semi_major_axis(small_io)
    assert orbit.get_stellar_distance(tiny_io) != orbit.get_semi_major_axis(tiny_io)
    assert orbit.get_stellar_distance(micro_io) != orbit.get_semi_major_axis(micro_io)

    # Check that the stellar eccentricities make sense
    orbit.set_eccentricity(io, .1, set_stellar_orbit=True)
    orbit.set_eccentricity(small_io, 0.01)
    orbit.set_eccentricity(tiny_io, 0.02)
    orbit.set_eccentricity(micro_io, 0.03)
    #    Stellar eccentricities should match the tidal host's
    assert orbit.get_stellar_eccentricity(small_io) == orbit.get_stellar_eccentricity(io)
    assert orbit.get_stellar_eccentricity(tiny_io) == orbit.get_stellar_eccentricity(io)
    assert orbit.get_stellar_eccentricity(micro_io) == orbit.get_stellar_eccentricity(io)
    #    They should not match the tidal world's eccentricity
    assert orbit.get_stellar_eccentricity(small_io) != orbit.get_eccentricity(small_io)
    assert orbit.get_stellar_eccentricity(tiny_io) != orbit.get_eccentricity(tiny_io)
    assert orbit.get_stellar_eccentricity(micro_io) != orbit.get_eccentricity(micro_io)
