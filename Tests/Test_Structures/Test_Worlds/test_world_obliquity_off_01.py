"""A tidal solve at nonzero obliquity with the obliquity truncation off ignores the obliquity and warns once per world,
naming the setting that turns the obliquity tides on."""
import math

import pytest

from TidalPy.constants import G
from TidalPy.Rheology.rheology import Elastic, Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.classes.tide import make_tide

from shared_materials import constant_solid


_WARNING_TEXT = "obliquity truncation is off"
_SETTING_TEXT = "set_tide_config(obliquity_truncation=...) or its [tides] obliquity_trunc_lvl"
_ORBIT = dict(orbital_frequency=2.0e-5, spin_frequency=3.0e-5, eccentricity=0.01, semi_major_axis=1.0e10,
              host_mass=2.0e30)


def _world(name, truncation="off"):
    world = StarWorld(name, 7.0e7, 1.9e27)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [1.0e4]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=truncation)
    return world


def test_warns_once_per_world_and_leaves_the_heating_unchanged(spdlog_text):
    world = _world("tilted")
    aligned = world.calc_tides(obliquity=0.0, **_ORBIT)["tidal_heating"]
    assert _WARNING_TEXT not in spdlog_text()

    tilted = world.calc_tides(obliquity=0.2, **_ORBIT)["tidal_heating"]
    world.calc_tides(obliquity=-0.3, **_ORBIT)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 1
    assert "world 'tilted'" in text and _SETTING_TEXT in text
    assert tilted == aligned

    _world("second").calc_tides(obliquity=0.2, **_ORBIT)
    assert spdlog_text().count(_WARNING_TEXT) == 2


@pytest.mark.parametrize("truncation", [2, 4, "gen"])
def test_no_warning_with_the_truncation_on(spdlog_text, truncation):
    world = _world("on", truncation)
    aligned = world.calc_tides(obliquity=0.0, **_ORBIT)["tidal_heating"]
    assert world.calc_tides(obliquity=0.2, **_ORBIT)["tidal_heating"] != aligned
    assert _WARNING_TEXT not in spdlog_text()


def test_a_system_solve_warns(spdlog_text):
    host = _world("host")
    planet = _world("planet")
    planet.set_obliquity(0.2)
    system = System("pair")
    system.add_world(host)
    system.add_world(planet, tidal_host="host", semi_major_axis=1.0e10, eccentricity=0.01)
    system.calc_world_evolution(planet)
    system.calc_world_evolution(planet)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 1 and "world 'planet'" in text


def test_the_3d_heating_warns(spdlog_text):
    radius, density = 1.0e6, 5000.0
    mass = (4.0 / 3.0) * math.pi * radius ** 3 * density
    material = constant_solid(density, bulk_modulus=1.0e11, shear_modulus=5.0e10, shear_viscosity=1.0e19,
                              bulk_viscosity=1.0e19)
    world = BaseWorld("solid", radius, mass)
    world.add_layer(Layer("mantle", 0, 0.0, radius, mass, material, is_static=False, shear_rheology=Maxwell(),
                          bulk_rheology=Elastic()))
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(obliquity_truncation="off")
    world.solve_eos(G_to_use=G)
    world.get_3d_tidal_heating(radius=0.5 * radius, colatitude=0.8, orbital_frequency=2.0e-5, spin_frequency=2.0e-5,
                               eccentricity=0.01, obliquity=0.2, semi_major_axis=1.0e9, host_mass=1.0e27)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 1 and "world 'solid'" in text
