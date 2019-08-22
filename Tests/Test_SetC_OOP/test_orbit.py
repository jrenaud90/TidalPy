import numpy as np

import TidalPy


TidalPy.verbose_level = 0
TidalPy.logging_level = 0
TidalPy.use_disk = False

from TidalPy.orbit import Orbit
from TidalPy.planets import build_planet
from TidalPy.utilities.conversions import days2rads, orbital_motion2semi_a, semi_a2orbital_motion


io = build_planet('io')
jupiter = build_planet('jupiter')
sol = build_planet('sol')
orbit = Orbit(sol, jupiter, io)


def test_build_orbit():
    assert isinstance(orbit, Orbit)


def test_orbit_object_getters():
    # Check that orbit was loaded into its objects correctly
    assert io.orbit is orbit
    assert jupiter.orbit is orbit
    assert sol.orbit is orbit

    # Check that the orbit can access all of the objects
    assert orbit.all_objects[0] is sol
    assert orbit.all_objects[1] is jupiter
    assert orbit.all_objects[2] is io
    assert orbit.all_objects_byname['sol'] is sol
    assert orbit.all_objects_byname['jupiter'] is jupiter
    assert orbit.all_objects_byname['io'] is io

    # Check planet pointer method
    assert orbit.find_planet_pointer(0) is sol
    assert orbit.find_planet_pointer(1) is jupiter
    assert orbit.find_planet_pointer(2) is io
    assert orbit.find_planet_pointer(sol) is sol
    assert orbit.find_planet_pointer(jupiter) is jupiter
    assert orbit.find_planet_pointer(io) is io
    assert orbit.find_planet_pointer('sol') is sol
    assert orbit.find_planet_pointer('jupiter') is jupiter
    assert orbit.find_planet_pointer('io') is io
    assert orbit.find_planet_pointer('Sol') is sol
    assert orbit.find_planet_pointer('Jupiter') is jupiter
    assert orbit.find_planet_pointer('Io') is io


def test_orbit_parameter_getters():
    eccen = 0.0041
    orb_freq = days2rads(1.769)

    # See if values were passed on to the objects
    np.testing.assert_approx_equal(io.eccentricity, eccen)
    np.testing.assert_approx_equal(io.orbital_freq, orb_freq)

    # See if the orbit can find the same values
    np.testing.assert_approx_equal(orbit.get_eccentricity(2), eccen)
    np.testing.assert_approx_equal(orbit.get_orbital_freq(2), orb_freq)
    np.testing.assert_approx_equal(orbit.get_semi_major_axis(2), orbital_motion2semi_a(orb_freq, jupiter.mass, io.mass))


def test_orbit_parameter_setters():
    new_eccen = 0.01
    new_semi_a = orbital_motion2semi_a(days2rads(1.), jupiter.mass, io.mass)

    # Set from orbit
    orbit.set_eccentricity(2, new_eccen)
    np.testing.assert_approx_equal(io.eccentricity, new_eccen)
    orbit.set_semi_major_axis(2, new_semi_a)
    np.testing.assert_approx_equal(io.semi_major_axis, new_semi_a)
    np.testing.assert_approx_equal(io.orbital_freq, semi_a2orbital_motion(new_semi_a, jupiter.mass, io.mass))

    # Set from planet
    newer_eccen = 0.1
    newer_semi_a = orbital_motion2semi_a(days2rads(1.5), jupiter.mass, io.mass)
    io.eccentricity = newer_eccen
    np.testing.assert_approx_equal(orbit.get_eccentricity(2), newer_eccen)
    io.semi_major_axis = newer_semi_a
    np.testing.assert_approx_equal(orbit.get_semi_major_axis(2), newer_semi_a)
    np.testing.assert_approx_equal(orbit.get_orbital_freq(2),
                                   semi_a2orbital_motion(newer_semi_a, jupiter.mass, io.mass))


def test_star_host():

    sh_orbit = Orbit(sol, None, io)
    assert sh_orbit.star_host
    sh_orbit = Orbit(sol, sol, io)
    assert sh_orbit.star_host
