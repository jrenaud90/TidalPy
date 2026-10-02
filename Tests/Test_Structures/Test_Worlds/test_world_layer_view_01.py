"""Non-owning layer views returned by a BaseWorld (access, sequence protocol, caching, lifetime)."""
import gc
import math

import pytest

from TidalPy.Material import Material
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.layers import Layer
from TidalPy.Tides.classes.tide import make_tide

from shared_materials import constant_solid


_R       = 1.6e6
_R_CORE  = 0.5 * _R


def _material(shear_modulus, bulk_modulus):
    """A constant-density elastic solid."""
    return constant_solid(3500.0, bulk_modulus=bulk_modulus, shear_modulus=shear_modulus)


def _two_layer_world():
    """A core and a mantle of different materials and tidal scales."""
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0
    world = BaseWorld("planet", _R, mass)
    core = Layer(
        "core",
        0,
        0.0,
        _R_CORE,
        0.0,
        _material(8.0e10, 2.5e11),
        tidal_scale=0.3,
    )
    mantle = Layer(
        "mantle",
        1,
        _R_CORE,
        _R,
        0.0,
        _material(6.0e10, 2.0e11),
        tidal_scale=0.7,
    )
    world.add_layer(core)
    world.add_layer(mantle)
    return world


@pytest.mark.parametrize("index, name", [(0, "core"), (1, "mantle")], ids=["core", "mantle"])
def test_layer_access_returns_the_named_layer(index, name):
    """get_layer and attribute access by name return a view of that layer."""
    world = _two_layer_world()
    assert isinstance(world.get_layer(index), Layer)
    assert isinstance(getattr(world, name), Layer)
    assert getattr(world, name).name == name
    assert world.get_layer(index).layer_index == index


def test_view_exposes_the_layer_api():
    """A view exposes the layer's geometry, settings, and material."""
    world = _two_layer_world()
    mantle = world.mantle
    assert mantle.radius_outer == _R
    assert math.isclose(mantle.tidal_scale, 0.7)
    assert math.isclose(mantle.calc_state(0.0)["shear_modulus"], 6.0e10)
    assert math.isclose(world.core.calc_state(0.0)["shear_modulus"], 8.0e10)


def test_layers_property_lists_views_inner_to_outer():
    """world.layers lists the views from inner to outer."""
    world = _two_layer_world()
    layers = world.layers
    assert len(layers) == 2
    assert [layer_view.name for layer_view in layers] == ["core", "mantle"]


def test_negative_index_and_out_of_range():
    """get_layer supports negative indices and raises IndexError out of range."""
    world = _two_layer_world()
    assert world.get_layer(-1).name == "mantle"
    with pytest.raises(IndexError):
        world.get_layer(2)
    with pytest.raises(IndexError):
        world.get_layer(-3)


def test_unknown_layer_name_raises_attributeerror():
    """An unknown attribute raises AttributeError; defined attributes win over layer lookup."""
    world = _two_layer_world()
    with pytest.raises(AttributeError):
        world.nonexistent_layer
    assert world.num_layers == 2


def test_world_is_iterable_over_layers():
    """Iterating a world yields its layers."""
    world = _two_layer_world()
    names = [layer.name for layer in world]
    assert names == ["core", "mantle"]


def test_len_and_indexing():
    """len(world) and world[i] follow the sequence protocol."""
    world = _two_layer_world()
    assert len(world) == 2
    assert world[0].name == "core"
    assert world[-1].name == "mantle"
    with pytest.raises(IndexError):
        world[2]


def test_slice_returns_view_list():
    """Slicing a world returns a list of views."""
    world = _two_layer_world()
    first = world[:1]
    assert [layer.name for layer in first] == ["core"]
    assert isinstance(first, list)


def test_iteration_yields_cached_views():
    """Iteration yields the same cached view objects as attribute access."""
    world = _two_layer_world()
    assert list(world)[1] is world.mantle


def test_views_are_cached_built_once():
    """Repeated access returns the same view object."""
    world = _two_layer_world()
    assert world.get_layer(1) is world.get_layer(1)
    assert world.mantle is world.mantle
    assert world.mantle is world.get_layer(1)
    assert world.layers[0] is world.get_layer(0)


def test_cache_invalidated_on_add_layer():
    """Adding a layer invalidates the view cache."""
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0
    world = BaseWorld("planet", _R, mass)
    world.add_layer(Layer(
        "core",
        0,
        0.0,
        _R_CORE,
        0.0,
    ))
    core_before = world.core
    world.add_layer(Layer(
        "mantle",
        1,
        _R_CORE,
        _R,
        0.0,
    ))
    assert world.core is not core_before
    assert world.mantle.name == "mantle"
    assert len(world.layers) == world.num_layers == 2


def test_view_keeps_world_alive_after_del():
    """A view keeps its world alive after the world name is deleted."""
    world = _two_layer_world()
    mantle = world.mantle
    del world
    gc.collect()
    assert mantle.name == "mantle"
    assert mantle.radius_outer == _R


def test_many_views_do_not_double_free():
    """Dropping many views leaves the world usable."""
    world = _two_layer_world()
    views = [world.get_layer(i % 2) for i in range(50)]
    del views
    gc.collect()
    assert world.mantle.name == "mantle"


def test_layer_view_reports_tidal_heating():
    """Layer-view tidal heating equals world heating times the layer's tidal_scale."""
    world = _two_layer_world()
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    world.calc_tides(
        orbital_frequency=2.05e-5,
        spin_frequency=2.05e-5,
        eccentricity=0.0041,
        obliquity=0.0,
        semi_major_axis=4.2e8,
        host_mass=1.898e27,
    )
    total = world.get_tidal_heating()
    assert math.isclose(world.mantle.get_tidal_heating(), total * 0.7, rel_tol=1.0e-9)
    assert math.isclose(world.core.get_tidal_heating(), total * 0.3, rel_tol=1.0e-9)
    assert math.isclose(world.get_layer_tidal_heating(1), world.mantle.get_tidal_heating(), rel_tol=1.0e-12)
