"""calc_dissipation, and calc_world_evolution and calc_pair_evolution built on it: one per-body tide calculation."""
import math

import pytest

from TidalPy.Material import Material, Phase
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.classes.tide import make_tide

AU = 1.495978707e11
TIDE_KEYS = ("tidal_heating", "dU_dM", "dU_dw", "dU_dO", "dU_dM_minus_dw", "moment_of_inertia")


def _body(name, radius, density, fixed_q=None):
    """A uniform layered body, with an analytic fixed-Q tide model when fixed_q is given (rigid otherwise)."""
    mass = 4.0 / 3.0 * math.pi * radius**3 * density
    world = BaseWorld(name, radius, mass)
    material = Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": density},
        shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10}))
    world.add_layer(Layer("mantle", 0, 0.0, radius, mass, material))
    if fixed_q is not None:
        world.set_tide_model(make_tide("fixed_q", {"fixed_q": [fixed_q]}))
        world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    world.solve_eos()
    return world


def _star_planet_moon(mutual):
    """The Sun lighting a planet and its moon. With `mutual` the planet and the moon host each other (the usual
    setup); otherwise the planet's tidal host is the Sun."""
    system = System("three")
    system.add_world(StarWorld("sun", 6.957e8, 1.988e30), is_star=True)
    planet = _body("planet", radius=6.0e6, density=5000.0, fixed_q=100.0)
    moon = _body("moon", radius=1.5e6, density=3300.0, fixed_q=50.0)
    planet.set_spin_frequency(7.0e-5)
    moon.set_spin_frequency(4.0e-6)
    if mutual:
        system.add_world(planet, stellar_semi_major_axis=AU, stellar_eccentricity=0.0167)
    else:
        system.add_world(planet, tidal_host="sun", semi_major_axis=AU, eccentricity=0.0167)
    system.add_world(
        moon,
        tidal_host="planet",
        semi_major_axis=3.844e8,
        eccentricity=0.0549,
        stellar_semi_major_axis=AU,
        stellar_eccentricity=0.0167)
    if mutual:
        system.set_tidal_host("planet", "moon")
    return system


# =====================================================================================================================
# calc_dissipation
# =====================================================================================================================
def test_dissipation_is_the_tide_part_of_the_world_evolution():
    system = _star_planet_moon(mutual=True)
    dissipation = system.calc_dissipation("moon")
    evolution = system.calc_world_evolution("moon")
    assert dissipation["solved"] and dissipation["has_tide_model"]
    assert (dissipation["world_name"], dissipation["companion_name"]) == ("moon", "planet")
    for key in TIDE_KEYS:
        assert dissipation[key] == evolution[key], key
    assert dissipation["companion_mass"] == evolution["host_mass"]
    assert dissipation["tidal_heating"] > 0.0
    assert "da_dt" not in dissipation


def test_a_world_without_a_host_is_not_solved():
    dissipation = _star_planet_moon(mutual=True).calc_dissipation("sun")
    assert dissipation["solved"] is False
    assert dissipation["companion_name"] is None


def test_a_rigid_world_is_solved_with_no_tide():
    system = _star_planet_moon(mutual=True)
    system.add_world(_body("rock", radius=1.0e5, density=2000.0), tidal_host="planet", semi_major_axis=1.0e8)
    dissipation = system.calc_dissipation("rock")
    assert dissipation["solved"] is True
    assert dissipation["has_tide_model"] is False
    assert dissipation["tidal_heating"] == 0.0


# =====================================================================================================================
# calc_pair_evolution on top of it
# =====================================================================================================================
@pytest.mark.parametrize("first, second", [("moon", "planet"), ("planet", "moon")])
def test_a_mutual_pair_is_each_world_s_own_dissipation(first, second):
    """Each member's part of the pair is its own dissipation, in either order, and the totals do not depend on it."""
    system = _star_planet_moon(mutual=True)
    pair = system.calc_pair_evolution(first, second)
    assert pair["world_names"] == (first, second)
    assert list(pair["worlds"]) == [first, second]
    for name, part in pair["worlds"].items():
        assert part["world_name"] == name
        assert part == system.calc_world_evolution(name)
    swapped = system.calc_pair_evolution(second, first)
    for key in ("da_dt", "de_dt", "dn_dt", "tidal_heating_total", "energy_residual"):
        assert swapped[key] == pytest.approx(pair[key], rel=1e-14, abs=0.0), key


def test_the_partner_defaults_to_the_tidal_host():
    system = _star_planet_moon(mutual=True)
    assert system.calc_pair_evolution("moon") == system.calc_pair_evolution("moon", "planet")


def test_a_one_sided_pair_takes_the_hosted_world_s_orbit():
    """The planet's own tidal host is the Sun, so in its pair with the moon it dissipates on the moon's orbit, not as
    calc_dissipation reports it (about the Sun)."""
    system = _star_planet_moon(mutual=False)
    pair = system.calc_pair_evolution("planet", "moon")
    planet_part = pair["worlds"]["planet"]
    assert planet_part["semi_major_axis"] == pytest.approx(3.844e8, rel=1e-15)
    assert planet_part["host_mass"] == system["moon"].mass
    assert planet_part["tidal_heating"] != system.calc_dissipation("planet")["tidal_heating"]
    reverse = system.calc_pair_evolution("moon")
    assert reverse["worlds"]["planet"]["tidal_heating"] == planet_part["tidal_heating"]
    assert reverse["da_dt"] == pytest.approx(pair["da_dt"], rel=1e-14, abs=0.0)


@pytest.mark.parametrize("first, second, match", [
    ("moon", "moon", "two different worlds"),
    ("sun", "moon", "neither"),
])
def test_a_pair_needs_two_worlds_one_hosting_the_other(first, second, match):
    system = _star_planet_moon(mutual=True)
    with pytest.raises(ValueError, match=match):
        system.calc_pair_evolution(first, second)


def test_a_world_without_a_host_has_no_pair():
    pair = _star_planet_moon(mutual=True).calc_pair_evolution("sun")
    assert pair["evolved"] is False
    assert pair["world_names"] == ("sun", None)
    assert pair["worlds"] == {}
