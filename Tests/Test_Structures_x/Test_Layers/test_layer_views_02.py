"""Layer views of a world stay safe and consistent when the layer or the world changes.

A view that outlives a world binary load raises instead of reading freed memory, moving a layer through a view
forgets the solved profile it no longer lines up with, and add_layer refuses a layer that would break the stack
while leaving the rejected layer usable.
"""
import math
import os

import pytest

from TidalPy.constants import G
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS

PLANET_RADIUS = 1.0e6  # [m]


def _layer(name, index, radius_inner, radius_outer, density=5000.0):
    layer = PhysicsLayer(name, index, radius_inner, radius_outer, 0.0)
    layer.set_eos(ConstantDensityEOS(
        reference_density=density,
        shear_modulus_static=5.0e10,
        bulk_modulus_static=1.0e11))
    return layer


def _world():
    world = LayeredWorld("w", PLANET_RADIUS, (4.0 / 3.0) * math.pi * PLANET_RADIUS ** 3 * 5000.0)
    world.add_layer(_layer("core", 0, 0.0, 0.5 * PLANET_RADIUS, 6000.0))
    world.add_layer(_layer("mantle", 1, 0.5 * PLANET_RADIUS, PLANET_RADIUS))
    return world


def test_views_held_across_a_binary_load_raise(tmp_path):
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
    world = LayeredWorld("w", PLANET_RADIUS, 1.0e22)
    added = _layer("only", 0, 0.0, PLANET_RADIUS)
    world.add_layer(added)
    path = os.path.join(tmp_path, "w.tpyb")
    world.save_binary(path)
    world.load_binary(path)
    with pytest.raises(RuntimeError):
        added.radius_outer


def test_moving_a_layer_through_a_view_forgets_the_solve():
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
    world.mantle.set_eos(ConstantDensityEOS(reference_density=3000.0))
    world.solve_eos(G_to_use=G)
    assert world.mantle.get_density(0.8 * PLANET_RADIUS) == pytest.approx(3000.0)


@pytest.mark.parametrize("layer_args, message", (
    (("core", 1, 0.5 * PLANET_RADIUS, PLANET_RADIUS), "already has a layer named"),
    (("mantle", 1, 0.5 * PLANET_RADIUS, 1.2 * PLANET_RADIUS), "past the radius"),
    (("mantle", 1, 0.4 * PLANET_RADIUS, PLANET_RADIUS), "not continuous"),
))
def test_add_layer_rejections_leave_the_layer_usable(layer_args, message):
    world = LayeredWorld("w", PLANET_RADIUS, 1.0e22)
    world.add_layer(_layer("core", 0, 0.0, 0.5 * PLANET_RADIUS))
    rejected = _layer(*layer_args)
    with pytest.raises(ValueError, match=message):
        world.add_layer(rejected)
    assert rejected.name == layer_args[0]
    assert world.num_layers == 1


@pytest.mark.parametrize("radii", ((-1.0, 1.0e6), (5.0e5, 1.0e5), (0.0, math.inf)))
def test_inverted_or_negative_layers_are_rejected(radii):
    with pytest.raises(ValueError, match="radius_inner"):
        PhysicsLayer("bad", 0, radii[0], radii[1], 0.0)
