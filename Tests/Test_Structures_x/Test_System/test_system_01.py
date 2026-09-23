"""System container (``TidalPy.structures_x.system.System``).

Covers world membership + co-ownership, host designation, per-world orbital elements, the Kepler
mean-motion / gravitational-parameter helpers, world identification by index/name/object, and the
sequence protocol (iterate / len / index / slice / attribute access).
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.structures_x.system import System
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld

MASS_SOLAR = 1.988435e30
RADIUS_SOLAR = 6.957e8
AU = 1.495978707e11


def _sun():
    return StarWorld("sun", RADIUS_SOLAR, MASS_SOLAR)


def _planet(name="earth", mass=5.9721986e24, radius=6.371e6):
    return LayeredWorld(name, radius, mass)


# =====================================================================================================================
# Membership + host
# =====================================================================================================================
def test_empty_system():
    system = System("sol")
    assert system.name == "sol"
    assert system.num_worlds == 0
    assert len(system) == 0
    assert system.has_star is False


def test_add_worlds_and_host():
    system = System()
    star = _sun()
    planet = _planet()
    idx_star = system.add_world(star)
    idx_planet = system.add_world(planet, tidal_host=star, semi_major_axis=AU, eccentricity=0.0167)
    assert idx_star == 0
    assert idx_planet == 1
    assert system.num_worlds == 2
    assert system.has_tidal_host(planet) is True
    assert system.get_tidal_host_index(planet) == 0
    assert system.get_tidal_host(planet) is star
    # The star is forced by nothing here.
    assert system.has_tidal_host(star) is False
    assert system.get_tidal_host_index(star) == -1
    assert system.get_tidal_host(star) is None


def test_set_tidal_host_by_index_name_object():
    system = System()
    star = _sun()
    planet = _planet()
    moon = _planet("moon", mass=7.342e22, radius=1.7374e6)
    system.add_world(star)
    system.add_world(planet)
    system.add_world(moon)
    system.set_tidal_host(2, 1)
    assert system.get_tidal_host(moon) is planet
    system.set_tidal_host("moon", "sun")
    assert system.get_tidal_host(moon) is star
    system.set_tidal_host(moon, planet)
    assert system.get_tidal_host("moon") is planet
    # None takes the host away again.
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
    """Two worlds that host each other describe one orbit: either may carry it, and both may if they agree."""
    system = System()
    earth = _planet("earth")
    moon = _planet("moon", mass=7.342e22, radius=1.7374e6)
    system.add_world(earth)
    system.add_world(moon, tidal_host=earth, semi_major_axis=3.844e8, eccentricity=0.0549)
    assert system.is_mutual_pair(moon) is False
    system.set_tidal_host(earth, moon)
    assert system.is_mutual_pair(earth) is True and system.is_mutual_pair(moon) is True

    # The Earth carries no elements of its own, so it reports the Moon's, and both see one mean motion.
    assert system.get_semi_major_axis(earth) == 3.844e8
    assert system.get_eccentricity(earth) == 0.0549
    assert system.calc_orbital_frequency(earth) == system.calc_orbital_frequency(moon)
    assert system.calc_gravitational_parameter(earth) == system.calc_gravitational_parameter(moon)

    # The same elements on both is fine; different ones contradict each other.
    system.set_semi_major_axis(earth, 3.844e8)
    system.set_eccentricity(earth, 0.0549)
    assert system.get_semi_major_axis(moon) == 3.844e8
    system.set_semi_major_axis(earth, 4.0e8)
    with pytest.raises(ValueError, match="share one orbit"):
        system.get_semi_major_axis(moon)
    with pytest.raises(ValueError, match="share one orbit"):
        system.calc_orbital_frequency(earth)


def test_add_none_world_raises():
    system = System()
    with pytest.raises(TypeError):
        system.add_world(None)


# =====================================================================================================================
# Orbital elements
# =====================================================================================================================
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
    system = System()
    star = _sun()
    earth_mass = 5.9721986e24
    system.add_world(star)
    system.add_world(_planet(mass=earth_mass), tidal_host=0, semi_major_axis=AU)

    mu = system.calc_gravitational_parameter("earth")
    assert math.isclose(mu, G * (MASS_SOLAR + earth_mass), rel_tol=1e-14)

    n = system.calc_orbital_frequency("earth")
    assert math.isclose(n, math.sqrt(mu / AU ** 3), rel_tol=1e-14)
    # Earth's mean motion is 2*pi / ~365.25 days.
    period_days = 2.0 * math.pi / n / 86400.0
    assert math.isclose(period_days, 365.25, rel_tol=2e-3)


def test_semi_major_axis_from_frequency_round_trip():
    system = System()
    system.add_world(_sun())
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    n = system.calc_orbital_frequency("earth")
    a = system.calc_semi_major_axis_from_frequency("earth", n)
    assert math.isclose(a, AU, rel_tol=1e-12)


def test_world_without_a_tidal_host_has_no_orbit():
    system = System()
    system.add_world(_sun())
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    # The star is hosted by nothing, so it has no orbit within the system.
    assert np.isnan(system.calc_gravitational_parameter(0))
    assert np.isnan(system.calc_orbital_frequency(0))

    # Elements with no host to be about give no orbit either.
    no_host = System()
    no_host.add_world(_planet(), semi_major_axis=AU)
    assert np.isnan(no_host.calc_gravitational_parameter(0))
    assert np.isnan(no_host.calc_orbital_frequency(0))


def test_unset_semi_major_axis_is_nan():
    system = System()
    system.add_world(_sun())
    system.add_world(_planet(), tidal_host=0)   # no semi-major axis supplied
    assert np.isnan(system.get_semi_major_axis("earth"))
    assert np.isnan(system.calc_orbital_frequency("earth"))


# =====================================================================================================================
# Sequence protocol + identification
# =====================================================================================================================
def test_sequence_protocol():
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
    # Attribute + name + object access all resolve to the same wrapper.
    assert system.planet_b is planet_b
    assert system["planet_c"] is planet_c
    assert system.worlds == [star, planet_b, planet_c]


def test_attribute_access_unknown_raises():
    system = System()
    system.add_world(_sun())
    with pytest.raises(AttributeError):
        system.not_a_world


def test_resolve_index_errors():
    system = System()
    system.add_world(_sun())
    with pytest.raises(IndexError):
        system[5]
    with pytest.raises(KeyError):
        system.get_semi_major_axis("nope")
    stray = _planet("stray")
    with pytest.raises(ValueError):
        system.get_semi_major_axis(stray)


def test_world_stays_usable_and_alive_after_add():
    """The added wrapper is co-owned; it stays fully functional and survives dropping the local ref."""
    system = System()
    star = _sun()
    system.add_world(star)
    # Star wrapper still works (co-owned C++ world).
    assert math.isclose(star.mass, MASS_SOLAR, rel_tol=1e-15)
    del star
    # The system still holds the world; its API is intact.
    held = system["sun"]
    assert math.isclose(held.mass, MASS_SOLAR, rel_tol=1e-15)
    assert held.name == "sun"


# =====================================================================================================================
# Star membership + stellar orbits + insolation
# =====================================================================================================================
def test_star_membership():
    """Exoplanet-like: the star is both the tidal host and the insolation source."""
    system = System()
    sun = _sun()
    system.add_world(sun, is_star=True)
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    assert system.has_star is True
    assert system.star_index == 0
    assert system.star is sun
    assert system.get_tidal_host("earth") is sun
    assert math.isclose(system.get_star_luminosity(), sun.luminosity, rel_tol=1e-15)


def test_tidal_host_and_star_are_independent():
    """Earth-Moon-Sun: the Moon's tidal host (Earth) and the star (Sun) are different bodies."""
    system = System()
    sun = _sun()
    earth = _planet("earth")
    moon = _planet("moon", mass=7.342e22, radius=1.7374e6)
    system.add_world(sun, is_star=True)                                   # star, nobody's tidal host
    system.add_world(earth)                                               # the Moon's tidal host, not the star
    system.add_world(moon, tidal_host=earth, semi_major_axis=3.844e8, eccentricity=0.0549)  # moon about Earth
    # Stellar orbits: Earth and the Moon both ~1 AU about the Sun.
    for name in ("earth", "moon"):
        system.set_stellar_semi_major_axis(name, AU)
        system.set_stellar_eccentricity(name, 0.0167)

    assert system.get_tidal_host(moon) is earth
    assert system.star is sun
    assert system.get_tidal_host_index(moon) != system.star_index
    # The tidal orbit (Moon about Earth) is distinct from the stellar orbit (Moon about Sun).
    assert math.isclose(system.get_semi_major_axis("moon"), 3.844e8, rel_tol=1e-12)
    assert math.isclose(system.get_stellar_semi_major_axis("moon"), AU, rel_tol=1e-12)


def test_stellar_orbit_getters_setters():
    system = System()
    system.add_world(_sun(), is_star=True)
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    system.set_stellar_semi_major_axis("earth", 1.2 * AU)
    system.set_stellar_eccentricity("earth", 0.05)
    assert math.isclose(system.get_stellar_semi_major_axis("earth"), 1.2 * AU, rel_tol=1e-15)
    assert math.isclose(system.get_stellar_eccentricity("earth"), 0.05, rel_tol=1e-15)


def test_stellar_gravitational_parameter_and_frequency():
    system = System()
    sun = _sun()
    earth_mass = 5.9721986e24
    system.add_world(sun, is_star=True)
    system.add_world(_planet(mass=earth_mass), tidal_host=0, semi_major_axis=AU)
    system.set_stellar_semi_major_axis("earth", AU)
    mu = system.calc_stellar_gravitational_parameter("earth")
    assert math.isclose(mu, G * (MASS_SOLAR + earth_mass), rel_tol=1e-14)
    n = system.calc_stellar_orbital_frequency("earth")
    assert math.isclose(n, math.sqrt(mu / AU ** 3), rel_tol=1e-14)


def test_insolation_flux_circular():
    system = System()
    sun = _sun()
    system.add_world(sun, is_star=True)
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    system.set_stellar_semi_major_axis("earth", AU)
    system.set_stellar_eccentricity("earth", 0.0)
    flux = system.calc_insolation_flux("earth")
    assert math.isclose(flux, sun.luminosity / (4.0 * math.pi * AU ** 2), rel_tol=1e-14)
    # The solar constant at 1 AU is ~1361 W/m^2.
    assert 1300.0 < flux < 1420.0


def test_insolation_flux_eccentric():
    system = System()
    sun = _sun()
    system.add_world(sun, is_star=True)
    system.add_world(_planet(), tidal_host=0, semi_major_axis=AU)
    system.set_stellar_semi_major_axis("earth", AU)
    system.set_stellar_eccentricity("earth", 0.4)
    flux = system.calc_insolation_flux("earth")
    expected = sun.luminosity / (4.0 * math.pi * AU ** 2 * math.sqrt(1.0 - 0.4 ** 2))
    assert math.isclose(flux, expected, rel_tol=1e-14)


def test_equilibrium_temperature():
    system = System()
    sun = _sun()
    earth = _planet()   # albedo 0.3, emissivity 1.0
    system.add_world(sun, is_star=True)
    system.add_world(earth, tidal_host=0, semi_major_axis=AU)
    system.set_stellar_semi_major_axis("earth", AU)
    flux = system.calc_insolation_flux("earth")
    t_system = system.calc_equilibrium_temperature("earth")
    # The system value equals the world's own radiative balance for that flux.
    assert math.isclose(t_system, earth.calc_equilibrium_temperature(flux), rel_tol=1e-14)
    # Earth's equilibrium temperature is ~255 K.
    assert 245.0 < t_system < 265.0


def test_insolation_guards():
    # No star set -> raises.
    no_star = System()
    no_star.add_world(_planet(), semi_major_axis=AU)
    with pytest.raises(RuntimeError):
        no_star.calc_insolation_flux("earth")

    system = System()
    system.add_world(_sun(), is_star=True)
    # Not tidally hosted by the star, so its orbit about the star is set on its own.
    system.add_world(_planet(), semi_major_axis=AU)
    # The star's own entry has no insolation.
    assert np.isnan(system.calc_insolation_flux(0))
    # Unset stellar semi-major axis -> NaN.
    assert np.isnan(system.calc_insolation_flux("earth"))
    assert np.isnan(system.calc_equilibrium_temperature("earth"))
    # Once the star is its tidal host the two orbits are one, and the insolation follows the tidal orbit.
    system.set_tidal_host("earth", 0)
    assert np.isfinite(system.calc_insolation_flux("earth"))


def test_set_star_by_index_name_object():
    system = System()
    sun = _sun()
    system.add_world(_planet("earth"))
    system.add_world(sun)
    system.set_star(1)
    assert system.star is sun
    system.set_star("sun")
    assert system.star is sun
    system.set_star(sun)
    assert system.star is sun
