"""Small world values (moment-of-inertia factor, tidal heat flux, a star's insolation flux, vectorized getters), the
keyword-only solve_eos with raise_on_fail, the solved-mass check, and what a saved configuration leaves out."""
import math

import numpy as np
import pytest

from TidalPy.exceptions import SolutionFailedError
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds import BaseWorld

from shared_materials import constant_solid

ORBIT = {"orbital_frequency": 4.11e-5, "spin_frequency": 4.11e-5, "eccentricity": 0.0041, "obliquity": 0.0,
         "semi_major_axis": 4.217e8, "host_mass": 1.898e27}


def uniform_world(name, density, stated_mass):
    """A one-layer world of constant density with a stated mass of its own."""
    world = BaseWorld(name, 1.0e6, stated_mass)
    world.add_layer(Layer("body", radius_outer=1.0e6, material=constant_solid(density, bulk_modulus=1.0e11,
                                                                              shear_modulus=5.0e10)))
    return world


# =====================================================================================================================
# Moment-of-inertia factor, tidal heat flux, insolation flux, vectorized getters (W11)
# =====================================================================================================================
def test_moment_of_inertia_factor():
    world = build_world("io")
    assert math.isnan(world.moment_of_inertia_factor)
    world.solve_eos()
    expected = world.planet_moi_eos / (world.planet_mass_eos * world.radius ** 2)
    assert world.moment_of_inertia_factor == pytest.approx(expected, rel=1.0e-14)
    assert 0.3 < world.moment_of_inertia_factor < 0.4


def test_tidal_heat_flux():
    world = build_world("io")
    assert math.isnan(world.get_tidal_heat_flux())
    world.solve_eos()
    world.calc_tides(**ORBIT)
    expected = world.get_tidal_heating() / (4.0 * math.pi * world.radius ** 2)
    assert world.get_tidal_heat_flux() == pytest.approx(expected, rel=1.0e-14)


def test_star_insolation_flux_matches_the_system():
    star = build_world("sol")
    planet = build_world("earth_simple")
    system = System("pair")
    system.add_world(star, is_star=True)
    system.add_world(planet, tidal_host=star, semi_major_axis=1.496e11, eccentricity=0.2)
    assert star.calc_insolation_flux(1.496e11, eccentricity=0.2) == system.calc_insolation_flux(planet)
    assert star.calc_insolation_flux(1.496e11) == pytest.approx(
        star.luminosity / (4.0 * math.pi * 1.496e11 ** 2), rel=1.0e-14)
    distances = np.array([1.0e11, 2.0e11])
    np.testing.assert_array_equal(
        star.calc_insolation_flux(distances), [star.calc_insolation_flux(value) for value in distances])
    with pytest.raises(ValueError, match="eccentricity in"):
        star.calc_insolation_flux(1.0e11, eccentricity=1.0)
    with pytest.raises(ValueError, match="positive distance"):
        star.calc_insolation_flux(0.0)


def test_equilibrium_temperature_takes_arrays():
    world = build_world("earth_simple")
    fluxes = np.array([[1361.0, 900.0], [0.0, 2000.0]])
    temperatures = world.calc_equilibrium_temperature(fluxes)
    assert temperatures.shape == fluxes.shape
    for index in np.ndindex(fluxes.shape):
        assert temperatures[index] == world.calc_equilibrium_temperature(float(fluxes[index]))
    assert isinstance(world.calc_equilibrium_temperature(1361.0), float)


def test_radial_functions_take_arrays(io):
    io.solve_love_numbers(frequency=4.1e-5)
    radii = np.linspace(0.3, 1.0, 5) * io.radius
    values = io.get_love_radial_y(radii, y_idx=1)
    assert values.dtype == np.complex128
    for radius, value in zip(radii, values):
        assert value == io.get_love_radial_y(float(radius), y_idx=1)


def test_3d_heating_array_takes_a_scalar_radius(io):
    world = io.copy()
    world.solve_eos()
    colatitudes = np.array([0.3, 1.0, 2.0])
    radius = 0.9 * world.radius
    np.testing.assert_array_equal(
        world.get_3d_tidal_heating_array(**ORBIT, radii=radius, colatitudes=colatitudes, num_threads=1),
        world.get_3d_tidal_heating_array(**ORBIT, radii=np.full(3, radius), colatitudes=colatitudes, num_threads=1))


# =====================================================================================================================
# Keyword-only solve_eos and raise_on_fail (W12)
# =====================================================================================================================
def test_solve_eos_is_keyword_only():
    world = build_world("io")
    with pytest.raises(TypeError):
        world.solve_eos(1500.0)


def test_solve_eos_g_defaults_to_the_configured_value():
    from TidalPy.constants import G
    world = build_world("io")
    default = world.solve_eos()
    explicit = world.solve_eos(G_to_use=G)
    assert default["planet_mass"] == explicit["planet_mass"]
    assert default["central_pressure"] == explicit["central_pressure"]


def test_solve_eos_raise_on_fail():
    world = build_world("io")
    failed = world.solve_eos(max_iters=1, pressure_tol=1.0e-14)
    assert not failed["success"]
    with pytest.raises(SolutionFailedError, match="EOS solve of world 'Io' failed"):
        world.solve_eos(max_iters=1, pressure_tol=1.0e-14, raise_on_fail=True)
    with pytest.raises(RuntimeError):
        world.solve_eos(max_iters=1, pressure_tol=1.0e-14, raise_on_fail=True)
    assert world.solve_eos(raise_on_fail=True)["success"]


# =====================================================================================================================
# Solved-mass check (W13)
# =====================================================================================================================
def test_a_solved_mass_off_the_stated_one_warns_once(spdlog_text):
    volume = (4.0 / 3.0) * math.pi * 1.0e18
    world = uniform_world("heavy_claim", 3000.0, 1.1 * 3000.0 * volume)
    assert world.solve_eos()["success"]
    world.solve_eos()
    text = spdlog_text()
    assert text.count("world 'heavy_claim' solved to a mass of") == 1
    assert "percent off its stated mass" in text


def test_a_matching_mass_does_not_warn(spdlog_text):
    volume = (4.0 / 3.0) * math.pi * 1.0e18
    world = uniform_world("honest", 3000.0, 1.005 * 3000.0 * volume)
    assert world.solve_eos()["success"]
    assert "solved to a mass" not in spdlog_text()


@pytest.mark.parametrize("name", ["io", "earth_simple", "earth_prem", "europa", "luna", "mercury", "pluto",
                                  "charon", "triton", "trappist1e", "jupiter_simple"])
def test_bundled_worlds_solve_to_their_stated_mass(spdlog_text, name):
    world = build_world(name)
    assert world.solve_eos()["success"]
    assert "solved to a mass" not in spdlog_text()


# =====================================================================================================================
# What a saved configuration leaves out (W15)
# =====================================================================================================================
def test_an_unsolved_world_writes_no_layer_masses_or_default_moment_of_inertia():
    world = build_world("earth_simple")
    config = world.get_config_dict()
    assert all("mass_kg" not in layer for layer in config["layers"].values())
    assert "moment_of_inertia_factor" not in config
    world.solve_eos()
    assert all(layer["mass_kg"] > 0.0 for layer in world.get_config_dict()["layers"].values())


def test_a_given_moment_of_inertia_factor_is_written():
    config = build_world("earth_simple").get_config_dict()
    config["moment_of_inertia_factor"] = 0.4
    rebuilt = build_world(config)
    assert rebuilt.get_config_dict()["moment_of_inertia_factor"] == 0.4
    config["moment_of_inertia_factor"] = 0.33
    assert build_world(config).get_config_dict()["moment_of_inertia_factor"] == 0.33
