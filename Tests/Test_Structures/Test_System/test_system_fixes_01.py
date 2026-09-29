"""Two element sets that describe one orbit merge element by element, and concurrent evolution calls that share a
world each read the result of their own tidal solve.
"""
import math
import threading

import pytest

from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Structures.configs.system_builder import build_system_from_dict
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.classes.tide import make_tide

AU = 1.495978707e11  # [m]


def _body(name, radius=1.0e6, density=3000.0, fixed_q=None):
    """A uniform layered body, with an analytic fixed-Q tide model when fixed_q is given."""
    mass = 4.0 / 3.0 * math.pi * radius**3 * density
    world = LayeredWorld(name, radius, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, radius, mass)
    layer.set_eos(ConstantDensityEOS(reference_density=density, shear_modulus_static=5.0e10))
    world.add_layer(layer)
    if fixed_q is not None:
        world.set_tide_model(make_tide("fixed_q", {"fixed_q": [fixed_q]}))
        world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    return world


def _star_and_planet(**planet_orbit):
    system = System("s")
    system.add_world(StarWorld("sun", 6.957e8, 1.988e30), is_star=True)
    planet = system.add_world(_body("planet"), **planet_orbit)
    return system, planet


# =====================================================================================================================
# A world hosted by the star
# =====================================================================================================================
def test_star_hosted_tidal_eccentricity_with_a_stellar_semi_major_axis():
    """The tidal set gives e and the stellar set gives a: the one orbit keeps both."""
    system, planet = _star_and_planet(eccentricity=0.1)
    system.set_stellar_semi_major_axis(planet, AU)
    system.set_tidal_host(planet, "sun")
    assert system.get_eccentricity(planet) == 0.1
    assert system.get_stellar_eccentricity(planet) == 0.1
    assert system.get_semi_major_axis(planet) == AU
    assert system.get_stellar_semi_major_axis(planet) == AU


def test_star_hosted_stellar_eccentricity_with_a_tidal_semi_major_axis():
    """The tidal set gives a and the stellar set gives e: the one orbit keeps both."""
    system, planet = _star_and_planet(semi_major_axis=AU)
    system.set_stellar_eccentricity(planet, 0.3)
    system.set_tidal_host(planet, "sun")
    assert system.get_eccentricity(planet) == 0.3
    assert system.get_stellar_eccentricity(planet) == 0.3
    assert system.get_semi_major_axis(planet) == AU


def test_star_hosted_elements_given_on_both_sets_must_agree():
    """An eccentricity both sets give with different values is refused, and the world keeps no host."""
    system, planet = _star_and_planet(semi_major_axis=AU, eccentricity=0.1)
    system.set_stellar_eccentricity(planet, 0.2)
    with pytest.raises(ValueError, match="different values"):
        system.set_tidal_host(planet, "sun")
    assert not system.has_tidal_host(planet)
    assert system.get_eccentricity(planet) == 0.1
    assert system.get_stellar_eccentricity(planet) == 0.2


def test_an_unset_eccentricity_reads_as_circular():
    system, planet = _star_and_planet(semi_major_axis=AU)
    assert system.get_eccentricity(planet) == 0.0
    assert system.get_stellar_eccentricity(planet) == 0.0
    system.set_tidal_host(planet, "sun")
    assert system.get_eccentricity(planet) == 0.0
    assert math.isfinite(system.calc_insolation_flux(planet))


@pytest.mark.parametrize("eccentricity", [-0.1, 1.0, math.inf])
def test_an_unbound_eccentricity_is_still_refused(eccentricity):
    system, planet = _star_and_planet(semi_major_axis=AU)
    with pytest.raises(ValueError, match="bound orbit"):
        system.set_eccentricity(planet, eccentricity)


# =====================================================================================================================
# A mutual pair
# =====================================================================================================================
def _mutual_pair(earth_orbit, moon_orbit):
    system = System("pair")
    earth = system.add_world(_body("earth", radius=6.371e6, density=5500.0), **earth_orbit)
    moon = system.add_world(_body("moon", radius=1.737e6), tidal_host="earth", **moon_orbit)
    system.set_tidal_host(earth, moon)
    return system


def test_mutual_pair_takes_each_element_from_the_member_that_gives_it():
    """The audit case: the eccentricity on the Earth and the semi-major axis on the Moon make one orbit."""
    system = _mutual_pair({"eccentricity": 0.055}, {"semi_major_axis": 3.844e8})
    for name in ("earth", "moon"):
        assert system.get_eccentricity(name) == 0.055
        assert system.get_semi_major_axis(name) == 3.844e8
    assert system.calc_orbital_frequency("earth") == system.calc_orbital_frequency("moon")


def test_mutual_pair_elements_given_by_both_must_agree():
    system = _mutual_pair({"eccentricity": 0.05}, {"semi_major_axis": 3.844e8, "eccentricity": 0.06})
    with pytest.raises(ValueError, match="share one orbit"):
        system.get_eccentricity("earth")
    system.set_eccentricity("earth", 0.06)
    assert system.get_eccentricity("earth") == 0.06


def test_mutual_pair_merged_orbit_survives_the_config_round_trip():
    """The live configuration writes the merged orbit on both members, so the rebuilt pair agrees."""
    system = _mutual_pair({"eccentricity": 0.055}, {"semi_major_axis": 3.844e8})
    config = system.get_config_dict()
    for name in ("earth", "moon"):
        assert config["worlds"][name]["eccentricity"] == 0.055
        assert config["worlds"][name]["semi_major_axis_m"] == 3.844e8
    rebuilt = build_system_from_dict(config)
    assert rebuilt.get_eccentricity("earth") == 0.055
    assert rebuilt.get_semi_major_axis("moon") == 3.844e8


def test_an_unset_eccentricity_survives_a_binary_round_trip(tmp_path):
    """An element left unset stays unset through a save and load, so the pair still merges after it."""
    system = _mutual_pair({"eccentricity": 0.055}, {"semi_major_axis": 3.844e8})
    path = str(tmp_path / "pair.tpyb")
    system.save_binary(path)
    loaded = System()
    loaded.load_binary(path)
    assert loaded.get_eccentricity("moon") == 0.055
    assert loaded.get_semi_major_axis("earth") == 3.844e8


# =====================================================================================================================
# Concurrent evolution calls
# =====================================================================================================================
def _three_body_system():
    """A planet about the star and a moon about the planet, both dissipating with an analytic tide model."""
    system = System("three")
    system.add_world(StarWorld("sun", 6.957e8, 1.988e30), is_star=True)
    planet = _body("planet", radius=6.0e6, density=5000.0, fixed_q=100.0)
    moon = _body("moon", radius=1.5e6, density=3300.0, fixed_q=50.0)
    for world in (planet, moon):
        world.solve_eos()
        world.set_spin_frequency(1.0e-5)
    system.add_world(planet, tidal_host="sun", semi_major_axis=0.05 * AU, eccentricity=0.05)
    system.add_world(moon, tidal_host="planet", semi_major_axis=4.0e8, eccentricity=0.02)
    return system


def test_concurrent_evolution_calls_each_read_their_own_solve():
    """The planet dissipates under the star in calc_world_evolution and under the moon in calc_pair_evolution; run
    together on two threads, each call reports exactly what it reports alone."""
    system = _three_body_system()
    world_reference = system.calc_world_evolution("planet")
    pair_reference = system.calc_pair_evolution("moon")
    assert world_reference["tidal_heating"] != pair_reference["host"]["tidal_heating"]
    mismatches = []
    errors = []

    def run(call, reference, pick):
        try:
            for _ in range(2000):
                result = pick(call())
                if result != pick(reference):
                    mismatches.append(result)
        except Exception as error:   # reported on the main thread
            errors.append(error)

    def host_heating(result):
        return (result["host"]["tidal_heating"], result["host"]["da_dt"], result["host"]["dspin_dt"])

    def world_heating(result):
        return (result["tidal_heating"], result["da_dt"], result["dspin_dt"])

    threads = [
        threading.Thread(target=run, args=(lambda: system.calc_world_evolution("planet"), world_reference,
                                           world_heating)),
        threading.Thread(target=run, args=(lambda: system.calc_pair_evolution("moon"), pair_reference,
                                           host_heating)),
    ]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert errors == []
    assert mismatches == []
