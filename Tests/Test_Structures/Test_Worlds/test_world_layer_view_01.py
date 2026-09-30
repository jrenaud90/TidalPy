"""Non-owning layer views returned by a LayeredWorld (dispatch, sequence protocol, caching, lifetime)."""
import gc
import math

import pytest

from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.layers.solidliquid import SolidLiquidLayer
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Tides.classes.tide import make_tide


_R       = 1.6e6
_R_CORE  = 0.5 * _R


def _two_layer_world():
    """A solidliquid core and physics mantle (distinct subclasses to test dispatch)."""
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0
    world = LayeredWorld("planet", _R, mass)
    core = SolidLiquidLayer(
        "core",
        0,
        0.0,
        _R_CORE,
        0.0,
        tidal_scale=0.3,
    )
    core.set_eos(ConstantDensityEOS(shear_modulus_static=8.0e10, bulk_modulus_static=2.5e11))
    mantle = BaseLayer(
        "mantle",
        1,
        _R_CORE,
        _R,
        0.0,
        tidal_scale=0.7,
    )
    mantle.set_eos(ConstantDensityEOS(shear_modulus_static=6.0e10, bulk_modulus_static=2.0e11))
    world.add_layer(core)
    world.add_layer(mantle)
    return world


@pytest.mark.parametrize(
    "index, name, layer_class",
    [(0, "core", SolidLiquidLayer), (1, "mantle", BaseLayer)],
    ids=["core", "mantle"],
)
def test_layer_access_dispatches_to_subclass(index, name, layer_class):
    """get_layer and attribute access by name return the matching layer subclass."""
    world = _two_layer_world()
    assert isinstance(world.get_layer(index), layer_class)
    assert isinstance(getattr(world, name), layer_class)
    assert getattr(world, name).name == name


def test_view_exposes_base_and_subclass_api():
    """A view exposes both the base-layer and the subclass API."""
    world = _two_layer_world()
    mantle = world.mantle
    assert mantle.radius_outer == _R
    assert math.isclose(mantle.tidal_scale, 0.7)
    assert math.isclose(mantle.shear_modulus_static, 6.0e10)
    assert math.isclose(world.core.shear_modulus_static, 8.0e10)


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
    world = LayeredWorld("planet", _R, mass)
    world.add_layer(BaseLayer(
        "core",
        0,
        0.0,
        _R_CORE,
        0.0,
    ))
    core_before = world.core
    world.add_layer(BaseLayer(
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
