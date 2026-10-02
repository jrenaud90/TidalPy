"""System container: membership, tidal host and star roles, orbital elements, Kepler helpers, and insolation."""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Structures.worlds.base import BaseWorld

MASS_SOLAR = 1.988435e30
RADIUS_SOLAR = 6.957e8
MASS_EARTH = 5.9721986e24
AU = 1.495978707e11


def _sun():
    return StarWorld("sun", RADIUS_SOLAR, MASS_SOLAR)


def _planet(name="earth", mass=MASS_EARTH, radius=6.371e6):
    return BaseWorld(name, radius, mass)


def _moon():
    return _planet("moon", mass=7.342e22, radius=1.7374e6)


def _star_and_earth(stellar_semi_major_axis=None, stellar_eccentricity=None):
    """The sun as star and tidal host of the earth at 1 AU."""
    system = System()
    sun = _sun()
    earth = _planet()
    system.add_world(sun, is_star=True)
    system.add_world(earth, tidal_host=0, semi_major_axis=AU)
    if stellar_semi_major_axis is not None:
        system.set_stellar_semi_major_axis("earth", stellar_semi_major_axis)
    if stellar_eccentricity is not None:
        system.set_stellar_eccentricity("earth", stellar_eccentricity)
    return system, sun, earth


def test_empty_system():
    system = System("sol")
    assert system.name == "sol"
    assert system.num_worlds == 0
    assert len(system) == 0
    assert system.has_star is False


def test_add_worlds_and_host():
    """add_world returns indices and records the tidal host; an unhosted world reports none."""
    system = System()
    star = _sun()
    planet = _planet()
    assert system.add_world(star) == 0
    assert system.add_world(planet, tidal_host=star, semi_major_axis=AU, eccentricity=0.0167) == 1
    assert system.num_worlds == 2
    assert system.has_tidal_host(planet) is True
    assert system.get_tidal_host_index(planet) == 0
    assert system.get_tidal_host(planet) is star
    assert system.has_tidal_host(star) is False
    assert system.get_tidal_host_index(star) == -1
    assert system.get_tidal_host(star) is None


def test_set_tidal_host_by_index_name_object():
    """set_tidal_host accepts indices, names, and objects, and None removes the host."""
    system = System()
    star = _sun()
    planet = _planet()
    moon = _moon()
    system.add_world(star)
    system.add_world(planet)
    system.add_world(moon)
    system.set_tidal_host(2, 1)
    assert system.get_tidal_host(moon) is planet
    system.set_tidal_host("moon", "sun")
    assert system.get_tidal_host(moon) is star
    system.set_tidal_host(moon, planet)
    assert system.get_tidal_host("moon") is planet
    system.set_tidal_host(moon, None)
    assert system.has_tidal_host(moon) is False


def test_a_world_cannot_host_itself():
    system = System()
    system.add_world(_sun())
    with pytest.raises(ValueError, match="own tidal host"):
        system.set_tidal_host("sun", "sun")


def test_add_world_with_an_unknown_host_adds_nothing():
    system = System()
    system.add_world(_sun())
    with pytest.raises(KeyError):
        system.add_world(_planet(), tidal_host="nobody", semi_major_axis=AU)
    assert system.num_worlds == 1


def test_mutual_pair_shares_one_orbit():
    """Two worlds hosting each other share one orbit: either may carry it, and both may if they agree."""
    system = System()
    earth = _planet("earth")
    moon = _moon()
    system.add_world(earth)
    system.add_world(moon, tidal_host=earth, semi_major_axis=3.844e8, eccentricity=0.0549)
    assert system.is_mutual_pair(moon) is False
    system.set_tidal_host(earth, moon)
    assert system.is_mutual_pair(earth) is True and system.is_mutual_pair(moon) is True

    # The earth carries no elements of its own, so it reports the moon's.
    assert system.get_semi_major_axis(earth) == 3.844e8
    assert system.get_eccentricity(earth) == 0.0549
    assert system.calc_orbital_frequency(earth) == system.calc_orbital_frequency(moon)
    assert system.calc_gravitational_parameter(earth) == system.calc_gravitational_parameter(moon)

    # Setting an element on either member sets the shared orbit.
    system.set_semi_major_axis(earth, 4.0e8)
    system.set_eccentricity(earth, 0.06)
    assert system.get_semi_major_axis(moon) == 4.0e8
    assert system.get_eccentricity(moon) == 0.06


def test_a_pair_formed_from_two_different_orbits_is_refused():
    """Two worlds that carried different elements before they hosted each other give no one orbit."""
    system = System()
    earth = _planet("earth")
    moon = _moon()
    system.add_world(earth)
    system.add_world(moon, tidal_host=earth, semi_major_axis=3.844e8, eccentricity=0.0549)
    system.set_semi_major_axis(earth, 4.0e8)
    system.set_tidal_host(earth, moon)
    with pytest.raises(ValueError, match="share one orbit"):
        system.get_semi_major_axis(moon)
    with pytest.raises(ValueError, match="share one orbit"):
        system.calc_orbital_frequency(earth)


def test_add_none_world_raises():
    with pytest.raises(TypeError):
        System().add_world(None)


def test_orbital_element_getters_setters():
    system = System()
    system.add_world(_sun())
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU, eccentricity=0.0167)
    assert math.isclose(system.get_semi_major_axis("earth"), AU, rel_tol=1e-15)
    assert math.isclose(system.get_eccentricity("earth"), 0.0167, rel_tol=1e-15)
    system.set_semi_major_axis("earth", 1.5 * AU)
    system.set_eccentricity("earth", 0.1)
    assert math.isclose(system.get_semi_major_axis(1), 1.5 * AU, rel_tol=1e-15)
    assert math.isclose(system.get_eccentricity(1), 0.1, rel_tol=1e-15)


def test_gravitational_parameter_and_mean_motion():
    """The tidal orbit's mu is G (M + m) and its mean motion follows Kepler's third law."""
    system = System()
    system.add_world(_sun())
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    mu = system.calc_gravitational_parameter("earth")
    assert math.isclose(mu, G * (MASS_SOLAR + MASS_EARTH), rel_tol=1e-14)
    mean_motion = system.calc_orbital_frequency("earth")
    assert math.isclose(mean_motion, math.sqrt(mu / AU ** 3), rel_tol=1e-14)
    period_days = 2.0 * math.pi / mean_motion / 86400.0
    assert math.isclose(period_days, 365.25, rel_tol=2e-3)


def test_semi_major_axis_from_frequency_round_trip():
    system = System()
    system.add_world(_sun())
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    mean_motion = system.calc_orbital_frequency("earth")
    assert math.isclose(system.calc_semi_major_axis_from_frequency("earth", mean_motion), AU, rel_tol=1e-12)


def test_world_without_a_tidal_host_has_no_orbit():
    """A world with no tidal host (with or without elements) has a NaN mu and mean motion."""
    system = System()
    system.add_world(_sun())
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    assert np.isnan(system.calc_gravitational_parameter(0))
    assert np.isnan(system.calc_orbital_frequency(0))

    no_host = System()
    no_host.add_world(_planet(), semi_major_axis=AU)
    assert np.isnan(no_host.calc_gravitational_parameter(0))
    assert np.isnan(no_host.calc_orbital_frequency(0))


def test_unset_semi_major_axis_is_nan():
    system = System()
    system.add_world(_sun())
    system.add_world(_planet(), tidal_host=0)
    assert np.isnan(system.get_semi_major_axis("earth"))
    assert np.isnan(system.calc_orbital_frequency("earth"))


def test_sequence_protocol():
    """A system iterates, indexes, slices, and resolves names and attributes to the added worlds."""
    system = System()
    star = _sun()
    planet_b = _planet("planet_b")
    planet_c = _planet("planet_c")
    system.add_world(star)
    system.add_world(planet_b, tidal_host=0, semi_major_axis=AU)
    system.add_world(planet_c, tidal_host=0, semi_major_axis=2 * AU)
    assert len(system) == 3
    assert list(system) == [star, planet_b, planet_c]
    assert system[0] is star
    assert system[-1] is planet_c
    assert system[1:] == [planet_b, planet_c]
    assert system.planet_b is planet_b
    assert system["planet_c"] is planet_c
    assert system.worlds == [star, planet_b, planet_c]


@pytest.mark.parametrize("lookup, error", [
    pytest.param(lambda system: system.not_a_world, AttributeError, id="unknown_attribute"),
    pytest.param(lambda system: system[5], IndexError, id="index_out_of_range"),
    pytest.param(lambda system: system.get_semi_major_axis("nope"), KeyError, id="unknown_name"),
    pytest.param(lambda system: system.get_semi_major_axis(_planet("stray")), ValueError, id="world_not_in_system"),
])
def test_bad_world_identifiers_raise(lookup, error):
    system = System()
    system.add_world(_sun())
    with pytest.raises(error):
        lookup(system)


def test_world_stays_usable_and_alive_after_add():
    """The system co-owns an added world, which survives dropping the local reference."""
    system = System()
    star = _sun()
    system.add_world(star)
    assert math.isclose(star.mass, MASS_SOLAR, rel_tol=1e-15)
    del star
    held = system["sun"]
    assert math.isclose(held.mass, MASS_SOLAR, rel_tol=1e-15)
    assert held.name == "sun"


def test_star_membership():
    """A star can be both the tidal host and the insolation source."""
    system, sun, _ = _star_and_earth()
    assert system.has_star is True
    assert system.star_index == 0
    assert system.star is sun
    assert system.get_tidal_host("earth") is sun
    assert math.isclose(system.get_star_luminosity(), sun.luminosity, rel_tol=1e-15)


def test_tidal_host_and_star_are_independent():
    """In an earth-moon-sun system the moon's tidal host (earth) differs from the star (sun)."""
    system = System()
    sun = _sun()
    earth = _planet("earth")
    moon = _moon()
    system.add_world(sun, is_star=True)
    system.add_world(earth)
    system.add_world(moon, tidal_host=earth, semi_major_axis=3.844e8, eccentricity=0.0549)
    for name in ("earth", "moon"):
        system.set_stellar_semi_major_axis(name, AU)
        system.set_stellar_eccentricity(name, 0.0167)

    assert system.get_tidal_host(moon) is earth
    assert system.star is sun
    assert system.get_tidal_host_index(moon) != system.star_index
    assert math.isclose(system.get_semi_major_axis("moon"), 3.844e8, rel_tol=1e-12)
    assert math.isclose(system.get_stellar_semi_major_axis("moon"), AU, rel_tol=1e-12)


def test_stellar_orbit_getters_setters():
    system, _, _ = _star_and_earth(stellar_semi_major_axis=1.2 * AU, stellar_eccentricity=0.05)
    assert math.isclose(system.get_stellar_semi_major_axis("earth"), 1.2 * AU, rel_tol=1e-15)
    assert math.isclose(system.get_stellar_eccentricity("earth"), 0.05, rel_tol=1e-15)


def test_stellar_gravitational_parameter_and_frequency():
    system, _, _ = _star_and_earth(stellar_semi_major_axis=AU)
    mu = system.calc_stellar_gravitational_parameter("earth")
    assert math.isclose(mu, G * (MASS_SOLAR + MASS_EARTH), rel_tol=1e-14)
    assert math.isclose(system.calc_stellar_orbital_frequency("earth"), math.sqrt(mu / AU ** 3), rel_tol=1e-14)


@pytest.mark.parametrize("eccentricity", [0.0, 0.4])
def test_insolation_flux(eccentricity):
    """The orbit-averaged flux is L / (4 pi a^2 sqrt(1 - e^2))."""
    system, sun, _ = _star_and_earth(stellar_semi_major_axis=AU, stellar_eccentricity=eccentricity)
    flux = system.calc_insolation_flux("earth")
    expected = sun.luminosity / (4.0 * math.pi * AU ** 2 * math.sqrt(1.0 - eccentricity ** 2))
    assert math.isclose(flux, expected, rel_tol=1e-14)
    if eccentricity == 0.0:
        # The solar constant at 1 AU is about 1361 W/m^2.
        assert 1300.0 < flux < 1420.0


def test_equilibrium_temperature():
    """The system value is the world's own radiative balance for its flux (about 255 K for the earth)."""
    system, _, earth = _star_and_earth(stellar_semi_major_axis=AU)
    flux = system.calc_insolation_flux("earth")
    temperature = system.calc_equilibrium_temperature("earth")
    assert math.isclose(temperature, earth.calc_equilibrium_temperature(flux), rel_tol=1e-14)
    assert 245.0 < temperature < 265.0


def test_insolation_guards():
    """Insolation raises with no star and is NaN for the star itself or an unset stellar orbit."""
    no_star = System()
    no_star.add_world(_planet(), semi_major_axis=AU)
    with pytest.raises(RuntimeError):
        no_star.calc_insolation_flux("earth")

    system = System()
    system.add_world(_sun(), is_star=True)
    # Not hosted by the star, so its stellar orbit is separate and left unset.
    system.add_world(_planet(), semi_major_axis=AU)
    assert np.isnan(system.calc_insolation_flux(0))
    assert np.isnan(system.calc_insolation_flux("earth"))
    assert np.isnan(system.calc_equilibrium_temperature("earth"))
    # With the star as its tidal host the two orbits are one, so the insolation follows the tidal orbit.
    system.set_tidal_host("earth", 0)
    assert np.isfinite(system.calc_insolation_flux("earth"))


@pytest.mark.parametrize("identifier", ["index", "name", "object"])
def test_set_star_by_index_name_object(identifier):
    system = System()
    sun = _sun()
    system.add_world(_planet("earth"))
    system.add_world(sun)
    system.set_star({"index": 1, "name": "sun", "object": sun}[identifier])
    assert system.star is sun
