"""Eccentricity rate precision at small e, and the checks a system makes on its orbits, names, and indices."""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Tides.classes.collapse import collapse_global_tides
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Dynamics import OrbitSolver
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.stellar import StarWorld

_N = 2.0 * np.pi / 86400.0
_HOST = 1.989e30
_MASS = 6.0e24
_SMA = (G * (_HOST + _MASS) / _N ** 2) ** (1.0 / 3.0)
_CPL = {"fixed_k": [0.3], "fixed_q": [100.0]}


def _collapse(eccentricity, spin_factor=3.0):
    return collapse_global_tides(
        planet_radius=6.4e6,
        orbital_frequency=_N,
        spin_frequency=spin_factor * _N,
        eccentricity=eccentricity,
        obliquity=0.0,
        semi_major_axis=_SMA,
        host_mass=_HOST,
        G_to_use=G,
        tide_model="cpl",
        tide_config=_CPL,
        eccentricity_truncation=4,
        obliquity_truncation=0)


def _de_dt(eccentricity, **kwargs):
    tide = _collapse(eccentricity)
    return OrbitSolver().calc_de_dt(_N, _SMA, eccentricity, _MASS, _HOST, tide["dUdM"], tide["dUdw"], **kwargs), tide


def _star(name="star", luminosity=3.828e26):
    return StarWorld(name, 7.0e8, _HOST, luminosity=luminosity)


def _planet(name="planet"):
    planet = StarWorld(name, 6.4e6, _MASS)  # a point mass with a fixed-Q tide is all these checks need
    planet.set_tide_model(make_tide("cpl", _CPL))
    return planet


def _star_system(**planet_orbit):
    system = System("checks")
    system.add_world(_star(), is_star=True)
    if planet_orbit:
        system.add_world(_planet(), **planet_orbit)
    return system


def test_the_per_mode_difference_is_the_difference_of_the_sums():
    tide = _collapse(0.05)
    assert tide["dUdM_minus_dw"] == pytest.approx(tide["dUdM"] - tide["dUdw"], rel=1.0e-10)


def test_de_dt_over_e_is_constant_as_e_vanishes():
    """With the per-mode difference, de/dt / e tends to a nonzero constant instead of losing precision."""
    rates = []
    for eccentricity in (1.0e-5, 1.0e-7, 1.0e-9, 1.0e-11):
        _, tide = _de_dt(eccentricity)
        exact = OrbitSolver().calc_de_dt(
            _N, _SMA, eccentricity, _MASS, _HOST, tide["dUdM"], tide["dUdw"], tide["dUdM_minus_dw"])
        rates.append(exact / eccentricity)
    assert rates[0] != 0.0
    for rate in rates[1:]:
        assert rate == pytest.approx(rates[0], rel=1.0e-8)


def test_both_forms_agree_at_a_moderate_eccentricity():
    de_dt, tide = _de_dt(0.05)
    exact, _ = _de_dt(0.05, dU_dM_minus_dw=tide["dUdM_minus_dw"])
    assert exact == pytest.approx(de_dt, rel=1.0e-9)


@pytest.mark.parametrize("orbit", [dict(eccentricity=1.0), dict(eccentricity=-0.1), dict(semi_major_axis=-1.0),
                                   dict(semi_major_axis=0.0), dict(semi_major_axis=math.inf)])
def test_an_unbound_orbit_is_refused(orbit):
    system = _star_system()
    with pytest.raises(ValueError):
        system.add_world(_planet(), tidal_host=0, **orbit)
    assert system.num_worlds == 1


@pytest.mark.parametrize("setter, value", [
    ("set_eccentricity", 1.2),
    ("set_semi_major_axis", -_SMA),
    ("set_stellar_eccentricity", -0.5),
])
def test_orbit_setters_refuse_an_unbound_orbit(setter, value):
    system = _star_system(tidal_host=0, semi_major_axis=_SMA)
    with pytest.raises(ValueError):
        getattr(system, setter)("planet", value)
    assert system.get_semi_major_axis("planet") == _SMA


@pytest.mark.parametrize("same_object", [True, False], ids=["same_world", "same_name"])
def test_names_are_unique_and_a_world_joins_once(same_object):
    system = System("checks")
    star = _star()
    system.add_world(star, is_star=True)
    with pytest.raises(ValueError):
        system.add_world(star if same_object else _star())
    assert system.num_worlds == 1


def test_world_indices_take_numpy_integers_but_not_bools():
    system = _star_system(tidal_host=0, semi_major_axis=_SMA)
    assert system[np.int64(1)] is system["planet"]
    assert system.get_semi_major_axis(np.int32(1)) == _SMA
    with pytest.raises(TypeError):
        system.get_semi_major_axis(True)


def test_a_star_hosted_orbit_moves_its_insolation():
    """A world whose tidal host is the star has one orbit, so moving it moves the insolation."""
    system = _star_system(tidal_host="star", semi_major_axis=_SMA)
    flux = system.calc_insolation_flux("planet")
    assert system.get_stellar_semi_major_axis("planet") == _SMA
    system.set_semi_major_axis("planet", 2.0 * _SMA)
    assert system.calc_insolation_flux("planet") == pytest.approx(flux / 4.0, rel=1.0e-12)
    system.set_stellar_semi_major_axis("planet", _SMA)
    assert system.get_semi_major_axis("planet") == _SMA
