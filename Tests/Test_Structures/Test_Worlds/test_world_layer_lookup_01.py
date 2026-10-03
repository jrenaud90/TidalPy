"""Layers given by index or name in every world method that takes one, and layers constructed without an index or an
inner radius."""
import math

import pytest

from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import BaseWorld
from TidalPy.Tides.classes.tide import make_tide

from shared_materials import constant_solid

RADIUS = 1.0e6
MASS = 1.0e22
ORBIT = {"orbital_frequency": 2.0e-5, "spin_frequency": 2.0e-5, "eccentricity": 0.01, "obliquity": 0.0,
         "semi_major_axis": 4.0e8, "host_mass": 1.9e27}


def three_layer_world():
    """A world built from layers given only their outer radii."""
    world = BaseWorld("stack", RADIUS, MASS)
    material = constant_solid(3000.0, bulk_modulus=1.0e11, shear_modulus=5.0e10)
    for name, fraction in (("core", 0.3), ("mantle", 0.8), ("crust", 1.0)):
        world.add_layer(Layer(name, radius_outer=fraction * RADIUS, material=material))
    world.set_tide_model(make_tide("cpl", fixed_k=[0.3], fixed_q=[50.0]))
    return world


# =====================================================================================================================
# Layer constructor and add_layer (W14)
# =====================================================================================================================
def test_add_layer_fills_the_index_and_inner_radius():
    world = three_layer_world()
    assert [layer.layer_index for layer in world] == [0, 1, 2]
    assert [layer.radius_inner for layer in world] == [0.0, 0.3 * RADIUS, 0.8 * RADIUS]


def test_positional_construction_still_works():
    world = BaseWorld("old_style", RADIUS, MASS)
    world.add_layer(Layer("core", 0, 0.0, 0.5 * RADIUS))
    world.add_layer(Layer("mantle", 1, 0.5 * RADIUS, RADIUS, 0.0))
    assert [layer.name for layer in world] == ["core", "mantle"]


def test_a_standalone_layer_without_them_reports_the_defaults():
    layer = Layer("lonely", radius_outer=RADIUS)
    assert layer.layer_index == 0
    assert layer.radius_inner == 0.0
    with pytest.raises(TypeError, match="radius_outer"):
        Layer("no_outer")


def test_a_given_index_or_inner_radius_that_disagrees_is_refused():
    world = BaseWorld("strict", RADIUS, MASS)
    world.add_layer(Layer("core", radius_outer=0.5 * RADIUS))
    wrong_index = Layer("mantle", 5, radius_outer=RADIUS)
    with pytest.raises(ValueError, match="layer_index 5, but it would be layer 1"):
        world.add_layer(wrong_index)
    # A refused layer is left as it was constructed and stays usable.
    assert wrong_index.layer_index == 5
    with pytest.raises(ValueError, match="not continuous.*starts at 400000"):
        world.add_layer(Layer("mantle", radius_inner=0.4 * RADIUS, radius_outer=RADIUS))
    with pytest.raises(ValueError, match="below the top of the stack"):
        world.add_layer(Layer("mantle", radius_outer=0.2 * RADIUS))
    world.add_layer(Layer("mantle", 1, radius_outer=RADIUS))
    assert world.num_layers == 2


# =====================================================================================================================
# Layers by index or name (W9)
# =====================================================================================================================
@pytest.mark.parametrize("key, name", [(0, "core"), (-1, "crust"), (-3, "core"), ("mantle", "mantle")])
def test_get_layer_and_indexing_take_an_index_or_a_name(key, name):
    world = three_layer_world()
    assert world.get_layer(key).name == name
    assert world[key].name == name


@pytest.mark.parametrize("method", ["get_layer", "get_layer_tidal_heating", "get_layer_tidal_scale",
                                    "calc_layer_temperature_rate"])
def test_lookup_errors_are_consistent(method):
    world = three_layer_world()
    with pytest.raises(IndexError, match="out of range"):
        getattr(world, method)(99)
    with pytest.raises(IndexError, match="out of range"):
        getattr(world, method)(-4)
    with pytest.raises(KeyError, match="no layer named 'mantel'.*mantle"):
        getattr(world, method)("mantel")
    with pytest.raises(TypeError, match="integer index"):
        getattr(world, method)(1.5)


def test_layer_tidal_heating_and_scale_by_name():
    world = three_layer_world()
    assert all(math.isnan(value) for value in world.get_layer_tidal_heating().values())
    world.calc_tides(**ORBIT)
    heating = world.get_layer_tidal_heating()
    assert list(heating) == ["core", "mantle", "crust"]
    assert heating["mantle"] == world.get_layer_tidal_heating(1) == world.get_layer_tidal_heating("mantle")
    assert world.get_layer_tidal_heating(-1) == heating["crust"]
    assert world.get_layer_tidal_scale("crust") == world.get_layer_tidal_scale(2)
    assert sum(heating.values()) == pytest.approx(world.get_tidal_heating(), rel=1.0e-12)
