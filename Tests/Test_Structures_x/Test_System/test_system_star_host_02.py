"""A world whose tidal host is the star has one orbit, whichever order its elements and roles are set in."""
import pytest

from TidalPy.structures_x.system import System
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld

AU = 1.495978707e11  # [m]


def _system_with_planet(semi_major_axis=None):
    system = System("s")
    system.add_world(StarWorld("sun", 6.957e8, 1.988e30), is_star=True)
    planet_index = system.add_world(LayeredWorld("earth", 6.371e6, 5.97e24), semi_major_axis=semi_major_axis)
    return system, planet_index


def test_stellar_orbit_set_before_the_host_becomes_the_tidal_orbit():
    system, planet = _system_with_planet()
    system.set_stellar_semi_major_axis(planet, AU)
    system.set_tidal_host(planet, "sun")
    assert system.get_semi_major_axis(planet) == pytest.approx(AU)
    assert system.get_stellar_semi_major_axis(planet) == pytest.approx(AU)


def test_two_different_orbits_are_refused():
    system, planet = _system_with_planet(semi_major_axis=AU)
    system.set_stellar_semi_major_axis(planet, 2.0 * AU)
    with pytest.raises(ValueError, match="differ"):
        system.set_tidal_host(planet, "sun")
    assert not system.has_tidal_host(planet)
