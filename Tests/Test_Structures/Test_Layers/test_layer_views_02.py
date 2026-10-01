"""Layer views stay safe across a world binary load and layer edits; add_layer and the constructor reject bad radii."""
import math
import os

import pytest

from TidalPy.constants import G
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.layers import Layer
from TidalPy.Material import Material, Phase

PLANET_RADIUS = 1.0e6  # [m]


def _material(density=5000.0):
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 1.0e11},
        shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10}))


def _layer(name, index, radius_inner, radius_outer, density=5000.0):
    return Layer(name, index, radius_inner, radius_outer, 0.0, _material(density))


def _world():
    world = BaseWorld("w", PLANET_RADIUS, (4.0 / 3.0) * math.pi * PLANET_RADIUS ** 3 * 5000.0)
    world.add_layer(_layer("core", 0, 0.0, 0.5 * PLANET_RADIUS, 6000.0))
    world.add_layer(_layer("mantle", 1, 0.5 * PLANET_RADIUS, PLANET_RADIUS))
    return world


def test_views_held_across_a_binary_load_raise(tmp_path):
    """A view that outlives a binary load raises instead of reading freed memory."""
    world = _world()
    held = [world.core, world.mantle]
    path = os.path.join(tmp_path, "w.tpyb")
    world.save_binary(path)
    world.load_binary(path)
    for view in held:
        with pytest.raises(RuntimeError, match="loaded a binary file"):
            view.name
    # Fresh views read the loaded layers.
    assert [layer.name for layer in world] == ["core", "mantle"]


def test_layer_passed_to_add_layer_is_detached_by_a_load(tmp_path):
    world = BaseWorld("w", PLANET_RADIUS, 1.0e22)
    added = _layer("only", 0, 0.0, PLANET_RADIUS)
    world.add_layer(added)
    path = os.path.join(tmp_path, "w.tpyb")
    world.save_binary(path)
    world.load_binary(path)
    with pytest.raises(RuntimeError):
        added.radius_outer


def test_moving_a_layer_through_a_view_forgets_the_solve():
    """Changing radii through a view drops the solved profile that no longer lines up."""
    world = _world()
    world.solve_eos(G_to_use=G)
    assert world.eos_solved
    mantle = world.mantle
    mantle.set_radii(0.5 * PLANET_RADIUS, PLANET_RADIUS)
    assert not world.eos_solved
    assert math.isnan(mantle.get_density(0.8 * PLANET_RADIUS))


def test_resolve_after_a_view_change_uses_the_new_material():
    world = _world()
    world.solve_eos(G_to_use=G)
    world.mantle.material = _material(3000.0)
    world.solve_eos(G_to_use=G)
    assert world.mantle.get_density(0.8 * PLANET_RADIUS) == pytest.approx(3000.0)


@pytest.mark.parametrize("layer_args, message", (
    (("core", 1, 0.5 * PLANET_RADIUS, PLANET_RADIUS), "already has a layer named"),
    (("mantle", 1, 0.5 * PLANET_RADIUS, 1.2 * PLANET_RADIUS), "past the radius"),
    (("mantle", 1, 0.4 * PLANET_RADIUS, PLANET_RADIUS), "not continuous"),
))
def test_add_layer_rejections_leave_the_layer_usable(layer_args, message):
    """add_layer refuses a layer that breaks the stack and leaves it usable."""
    world = BaseWorld("w", PLANET_RADIUS, 1.0e22)
    world.add_layer(_layer("core", 0, 0.0, 0.5 * PLANET_RADIUS))
    rejected = _layer(*layer_args)
    with pytest.raises(ValueError, match=message):
        world.add_layer(rejected)
    assert rejected.name == layer_args[0]
    assert world.num_layers == 1


@pytest.mark.parametrize("radii", ((-1.0, 1.0e6), (5.0e5, 1.0e5), (0.0, math.inf)))
def test_inverted_or_negative_layers_are_rejected(radii):
    with pytest.raises(ValueError, match="radius_inner"):
        Layer("bad", 0, radii[0], radii[1], 0.0)
