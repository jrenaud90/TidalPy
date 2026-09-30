"""A layered world's solved state: what invalidates it, and what a released solution keeps."""
import math
import os

import pytest

from TidalPy.Material.eos import make_material_eos
from TidalPy.PartialMelt import make_partial_melt
from TidalPy.Structures import build_world
from TidalPy.Structures.layers.base import BaseLayer

_IO_FREQUENCY = 4.11e-5   # [rad s-1]
_MANTLE_RADIUS = 1.2e6    # [m], inside Io's mantle


def _solved_io(**eos_kwargs):
    world = build_world("io")
    world.solve_eos(**eos_kwargs)
    return world


def test_a_released_solution_does_not_follow_a_later_solve():
    """A released radial solution keeps its structure and moduli after the world re-solves."""
    world = _solved_io(temperature=1600.0)
    world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)
    released = world.release_radial_solution()
    density = released.get_density(_MANTLE_RADIUS)
    shear = released.get_shear_modulus(_MANTLE_RADIUS)
    viscosity = released.get_shear_viscosity(_MANTLE_RADIUS)

    world.solve_eos(temperature=1750.0)
    world.solve_eos(nondimensionalize=False)
    assert released.get_density(_MANTLE_RADIUS) == density
    assert released.get_shear_modulus(_MANTLE_RADIUS) == shear
    assert released.get_shear_viscosity(_MANTLE_RADIUS) == viscosity


def test_changing_a_layer_material_after_a_solve_unsolves_the_world():
    """A material change clears the solved profile, so no Love solve can use the stale one."""
    world = _solved_io()
    density_before = world.get_density(_MANTLE_RADIUS)
    world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)

    world.mantle.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 2000.0}))
    assert not world.eos_solved
    assert math.isnan(world.get_density(_MANTLE_RADIUS))
    with pytest.raises(ValueError):
        world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)
    world.solve_eos()
    assert world.get_density(_MANTLE_RADIUS) != density_before


def test_a_new_melt_model_takes_effect_at_the_next_solve():
    """A new partial-melt model unsolves the world; the next solve_eos applies it."""
    world = _solved_io()
    shear_before = world.get_shear_modulus(_MANTLE_RADIUS)
    world.mantle.set_partial_melt(make_partial_melt("henning", {"solidus_k": 100.0, "liquidus_k": 200.0}))
    assert not world.eos_solved
    assert math.isnan(world.get_shear_modulus(_MANTLE_RADIUS))
    assert world.molten_regions == []
    world.solve_eos()
    assert world.get_shear_modulus(_MANTLE_RADIUS) < shear_before


def test_loading_a_binary_file_leaves_the_world_unsolved(tmp_path):
    """Loading a binary file leaves the world unsolved; a fresh solve matches the original."""
    world = _solved_io()
    k2_original = world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)["love_number_k"]
    path = os.path.join(str(tmp_path), "io.tpyb")
    world.save_binary(path)

    world.load_binary(path)
    assert not world.eos_solved
    assert not world.love_solved
    assert math.isnan(world.get_density(_MANTLE_RADIUS))
    assert world.mantle.radius > 0.0
    with pytest.raises(ValueError):
        world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)
    world.solve_eos()
    k2_reloaded = world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)["love_number_k"]
    assert k2_reloaded == pytest.approx(k2_original, rel=1e-12)


@pytest.mark.parametrize("name", ["io", "earth_prem", "sol"])
def test_a_load_clears_the_configuration_the_world_was_built_from(name, tmp_path):
    """After a load the build configurations no longer describe the world, so they are cleared."""
    world = build_world(name)
    assert world.config is not None
    path = os.path.join(str(tmp_path), f"{name}.tpyb")
    world.save_binary(path)

    world.load_binary(path)
    assert world.config is None
    assert world.portable_config is None


def test_a_layer_view_cannot_be_loaded_in_place(tmp_path):
    """load_binary on a world-owned layer view raises."""
    world = _solved_io()
    path = os.path.join(str(tmp_path), "mantle.tpyb")
    world.mantle.save_binary(path)
    with pytest.raises(ValueError):
        world.mantle.load_binary(path)


def test_a_layer_past_the_world_radius_is_refused():
    """A layer reaching past the world radius is refused and the world stays solved."""
    world = _solved_io()
    top = world.radius
    shell = BaseLayer(
        "shell",
        world.num_layers,
        top,
        top + 1.0e4,
        0.0,
    )
    shell.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 1000.0}))
    with pytest.raises(ValueError, match="past the radius"):
        world.add_layer(shell)
    assert world.eos_solved
    assert math.isfinite(world.get_density(_MANTLE_RADIUS))


def test_moving_a_layer_leaves_the_world_unsolved():
    """Moving a layer leaves the world unsolved and its profile reads NaN."""
    world = _solved_io()
    mantle = world.mantle
    mantle.set_radii(mantle.radius_inner, mantle.radius_outer)
    assert not world.eos_solved
    assert math.isnan(world.get_density(_MANTLE_RADIUS))


def test_a_failed_solve_leaves_the_world_unsolved():
    """A failed re-solve leaves the world unsolved until a good solve restores it."""
    world = _solved_io()
    result = world.solve_eos(integration_method="LSODA", rtol=1.0e-9, atol=1.0e-12)
    if result["success"]:
        pytest.skip("LSODA solved the whole-planet EOS here, so this check has no failing solve to use.")
    assert not world.eos_solved
    assert math.isnan(world.get_density(_MANTLE_RADIUS))
    assert all(math.isnan(value) for value in result["density"])
    with pytest.raises(ValueError):
        world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)
    world.solve_eos()
    assert world.eos_solved
    assert world.get_density(_MANTLE_RADIUS) > 0.0
