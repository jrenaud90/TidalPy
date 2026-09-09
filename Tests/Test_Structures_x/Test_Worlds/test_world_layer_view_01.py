"""Tests for the non-owning layer views returned by a LayeredWorld.

A built world owns its layers as C++ objects. ``world.get_layer(i)``, ``world.layers``, and
attribute access ``world.<layer_name>`` hand back a Cython wrapper around the world-owned layer
(dispatched to the matching ``PhysicsLayer`` / ``SolidLiquidLayer`` / ``GasLayer`` / ``BaseLayer``
subclass) so the layer's full API is reachable. The view is non-owning and keeps the world alive.
"""
import gc
import math

import pytest

from TidalPy.constants import G
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.base import BaseLayer
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.layers.solidliquid import SolidLiquidLayer
from TidalPy.Tides_x.classes.tide import make_tide


_R       = 1.6e6
_R_CORE  = 0.5 * _R


def _two_layer_world():
    """A solidliquid core + physics mantle (distinct subclasses to test dispatch)."""
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0
    world = LayeredWorld("planet", _R, mass)
    core = SolidLiquidLayer("core", 0, 0.0, _R_CORE, 0.0,
                            shear_modulus_static_pa=8.0e10, bulk_modulus_static_pa=2.5e11,
                            tidal_scale=0.3)
    mantle = PhysicsLayer("mantle", 1, _R_CORE, _R, 0.0,
                          shear_modulus_static_pa=6.0e10, bulk_modulus_static_pa=2.0e11,
                          tidal_scale=0.7)
    world.add_layer(core)
    world.add_layer(mantle)
    return world


# =====================================================================================================================
# Index + attribute access, subclass dispatch
# =====================================================================================================================
def test_get_layer_dispatches_to_subclass():
    world = _two_layer_world()
    assert isinstance(world.get_layer(0), SolidLiquidLayer)
    assert isinstance(world.get_layer(1), PhysicsLayer)


def test_attribute_access_by_layer_name():
    world = _two_layer_world()
    assert isinstance(world.core, SolidLiquidLayer)
    assert isinstance(world.mantle, PhysicsLayer)
    assert world.core.name == "core"
    assert world.mantle.name == "mantle"


def test_view_exposes_base_and_subclass_api():
    world = _two_layer_world()
    mantle = world.mantle
    # Base-layer API:
    assert mantle.radius_outer == _R
    assert math.isclose(mantle.tidal_scale, 0.7)
    assert mantle.tidal_scale_method == "user_provided_scale"
    # PhysicsLayer-specific API:
    assert math.isclose(mantle.shear_modulus_static, 6.0e10)
    assert math.isclose(world.core.shear_modulus_static, 8.0e10)


def test_layers_property_lists_views_inner_to_outer():
    world = _two_layer_world()
    layers = world.layers
    assert len(layers) == 2
    assert [layer_view.name for layer_view in layers] == ["core", "mantle"]


def test_negative_index_and_out_of_range():
    world = _two_layer_world()
    assert world.get_layer(-1).name == "mantle"
    with pytest.raises(IndexError):
        world.get_layer(2)
    with pytest.raises(IndexError):
        world.get_layer(-3)


def test_unknown_layer_name_raises_attributeerror():
    world = _two_layer_world()
    with pytest.raises(AttributeError):
        world.nonexistent_layer
    # A defined attribute/method must still win over the layer lookup.
    assert world.num_layers == 2


# =====================================================================================================================
# Sequence protocol: for layer in world / len(world) / world[i]
# =====================================================================================================================
def test_world_is_iterable_over_layers():
    world = _two_layer_world()
    names = [layer.name for layer in world]
    assert names == ["core", "mantle"]


def test_len_and_indexing():
    world = _two_layer_world()
    assert len(world) == 2
    assert world[0].name == "core"
    assert world[-1].name == "mantle"
    with pytest.raises(IndexError):
        world[2]


def test_slice_returns_view_list():
    world = _two_layer_world()
    first = world[:1]
    assert [layer.name for layer in first] == ["core"]
    assert isinstance(first, list)


def test_iteration_yields_cached_views():
    world = _two_layer_world()
    assert list(world)[1] is world.mantle


# =====================================================================================================================
# Views are built once (cached) and the cache is invalidated when layers change
# =====================================================================================================================
def test_views_are_cached_built_once():
    world = _two_layer_world()
    # Repeated access returns the same object (no per-access reallocation).
    assert world.get_layer(1) is world.get_layer(1)
    assert world.mantle is world.mantle
    assert world.mantle is world.get_layer(1)
    assert world.layers[0] is world.get_layer(0)


def test_cache_invalidated_on_add_layer():
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0
    world = LayeredWorld("planet", _R, mass)
    world.add_layer(PhysicsLayer("core", 0, 0.0, _R_CORE, 0.0))
    core_before = world.core               # builds the cache
    world.add_layer(PhysicsLayer("mantle", 1, _R_CORE, _R, 0.0))   # invalidates it
    assert world.core is not core_before   # rebuilt
    assert world.mantle.name == "mantle"
    assert len(world.layers) == world.num_layers == 2


# =====================================================================================================================
# Non-owning semantics: the view keeps the world alive and never double-frees
# =====================================================================================================================
def test_view_keeps_world_alive_after_del():
    world = _two_layer_world()
    mantle = world.mantle
    del world
    gc.collect()
    # The view holds a reference to the world, so the C++ layer is still valid.
    assert mantle.name == "mantle"
    assert mantle.radius_outer == _R


def test_many_views_do_not_double_free():
    world = _two_layer_world()
    views = [world.get_layer(i % 2) for i in range(50)]
    del views
    gc.collect()
    # World still usable after dropping all views (no double-free / corruption).
    assert world.mantle.name == "mantle"


# =====================================================================================================================
# Tidal heating is reachable on the layer view
# =====================================================================================================================
def test_layer_view_reports_tidal_heating():
    world = _two_layer_world()
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=2, obliquity_truncation=0)
    world.calc_tides(orbital_frequency=2.05e-5, spin_frequency=2.05e-5, eccentricity=0.0041,
                     obliquity=0.0, semi_major_axis=4.2e8, host_mass=1.898e27)
    total = world.get_tidal_heating()
    # layer.get_tidal_heating() (on the view) == world heating * the layer's tidal_scale.
    assert math.isclose(world.mantle.get_tidal_heating(), total * 0.7, rel_tol=1.0e-9)
    assert math.isclose(world.core.get_tidal_heating(), total * 0.3, rel_tol=1.0e-9)
    # Agrees with the world-side accessor.
    assert math.isclose(world.get_layer_tidal_heating(1), world.mantle.get_tidal_heating(),
                        rel_tol=1.0e-12)
