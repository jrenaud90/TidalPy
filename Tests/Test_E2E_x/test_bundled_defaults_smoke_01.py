"""Defaults smoke tests: every bundled world through every flagship entry point, with no argument left to chance.

Each entry point is called the way a first-time user calls it, with only its required arguments, so whatever the
package defaults and the ``_x`` configuration resolve to is what runs. The checks are coarse on purpose (it ran, it
succeeded, the numbers are finite and have the right sign); the point is that a default which breaks a bundled world
fails here, whichever module it lives in.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.structures_x import available_systems, available_worlds, build_system, build_world
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.worlds.stellar import StarWorld


_HOST_MASS = 1.898e27     # [kg]
_ORBITAL_PERIOD_DAYS = 3.0
_ECCENTRICITY = 0.01


def _tidal_state(world):
    """A synchronous, slightly eccentric orbit about a Jupiter-mass host."""
    mean_motion = 2.0 * math.pi / (_ORBITAL_PERIOD_DAYS * 86400.0)
    semi_major_axis = (G * (_HOST_MASS + world.mass) / mean_motion**2) ** (1.0 / 3.0)
    return dict(
        orbital_frequency=mean_motion, spin_frequency=mean_motion, eccentricity=_ECCENTRICITY, obliquity=0.0,
        semi_major_axis=semi_major_axis, host_mass=_HOST_MASS)


_WORLD_NAMES = available_worlds()
_LAYERED_NAMES = [name for name in _WORLD_NAMES if isinstance(build_world(name), LayeredWorld)]


def test_the_pack_holds_the_worlds_the_documentation_names():
    for name in ("earth_prem", "earth_simple", "europa", "io", "jupiter_simple", "luna", "mercury", "sol"):
        assert name in _WORLD_NAMES


@pytest.mark.parametrize("world_name", _WORLD_NAMES)
def test_build_and_bulk_properties(world_name):
    world = build_world(world_name)
    assert world.radius > 0.0 and world.mass > 0.0
    assert math.isfinite(world.calc_surface_gravity()) and world.calc_surface_gravity() > 0.0
    assert world.tide_model_set, "every bundled world is given a tide model by the builder"
    if isinstance(world, StarWorld):
        assert world.luminosity > 0.0 and world.effective_temperature > 0.0


@pytest.mark.parametrize("world_name", _LAYERED_NAMES)
def test_solve_eos_defaults(world_name):
    world = build_world(world_name)
    result = world.solve_eos()
    assert world.eos_solved and result["success"], result["message"]
    assert not result["max_iters_hit"]
    # Coarse on purpose: the worlds fit to their mass are held to it in test_bundled_bodies_01, and the simple
    # worlds are not fit at all.
    assert 0.5 * world.mass < world.planet_mass_eos < 2.0 * world.mass
    factor = world.planet_moi_eos / (world.planet_mass_eos * world.radius**2)
    assert 0.2 < factor <= 0.4 * (1.0 + 1.0e-9)   # 0.4 is a uniform sphere, which jupiter_simple is
    radii = np.linspace(0.0, world.radius, 25)
    for getter in (world.get_density, world.get_pressure, world.get_gravity):
        values = np.asarray(getter(radii))
        assert np.all(np.isfinite(values))
        # The surface pressure is zero to the tolerance of the solve, so it can land a hair below.
        assert np.all(values >= -1.0e-6 * values.max())
    assert result["planet_mass"] == world.planet_mass_eos


@pytest.mark.parametrize("world_name", _LAYERED_NAMES)
def test_solve_love_numbers_defaults(world_name):
    world = build_world(world_name)
    world.solve_eos()
    result = world.solve_love_numbers()
    assert result["success"], result["message"]
    k2, h2 = result["love_number_k"], result["love_number_h"]
    # Between a rigid body and a homogeneous fluid (3/2 and 5/2, which a gas giant reaches to roundoff), and
    # dissipating or lossless but never amplifying.
    assert 0.0 < k2.real < 1.5 * (1.0 + 1.0e-6)
    assert k2.imag <= 0.0
    # A static liquid carries no displacement solution, so a world whose surface is one has k alone.
    surface_layer = world[-1]
    if surface_layer.is_solid or not surface_layer.is_static:
        assert 0.0 < h2.real < 2.5 * (1.0 + 1.0e-6)
    else:
        assert math.isnan(h2.real)


@pytest.mark.parametrize("world_name", _WORLD_NAMES)
def test_calc_tides_defaults(world_name):
    world = build_world(world_name)
    if isinstance(world, LayeredWorld):
        world.solve_eos()
    world.calc_tides(**_tidal_state(world))
    assert world.tides_solved
    heating = world.get_tidal_heating()
    assert math.isfinite(heating) and heating >= 0.0
    assert world.get_num_tidal_modes() > 0
    assert all(math.isfinite(value) for value in world.get_tidal_potential_derivatives())
    if isinstance(world, LayeredWorld):
        per_layer = [world.get_layer_tidal_heating(index) for index in range(world.num_layers)]
        assert all(math.isfinite(value) and value >= 0.0 for value in per_layer)


@pytest.mark.parametrize("world_name", _LAYERED_NAMES)
def test_calc_3d_tides_defaults(world_name):
    world = build_world(world_name)
    if world.get_config_dict()["tides"]["global_tidal_model"] != "rheology":
        pytest.skip("The 3D paths need the rheology tide model.")
    world.solve_eos()
    state = _tidal_state(world)
    world.calc_tides(**state)
    summed = world.calc_3d_tides(**state, radial_summed=True, latitude_summed=True, longitude_summed=True)
    assert math.isfinite(summed["total"]) and summed["total"] >= 0.0
    # The volume integral of the 3D heating is the 1D heating, on the default grids.
    if world.get_tidal_heating() > 0.0:
        assert summed["total"] == pytest.approx(world.get_tidal_heating(), rel=5.0e-2)
    assert len(summed["per_layer"]) == world.num_layers


@pytest.mark.parametrize("system_name", available_systems())
def test_system_evolution_defaults(system_name):
    system = build_system(system_name)
    for world in system:
        if isinstance(world, LayeredWorld):
            world.solve_eos()
    rates = system.calc_system_evolution()
    assert len(rates) == system.num_worlds
    evolved = [entry for entry in rates if entry["evolved"]]
    assert evolved, "a bundled system has at least one tidally forced world"
    for entry in evolved:
        for key, value in entry.items():
            if isinstance(value, float):
                assert math.isfinite(value), key
