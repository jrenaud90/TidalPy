"""Replacing a layer's material keeps the viscosity and partial-melt models attached through the layer.

The layer setters (`set_shear_viscosity`, `set_bulk_viscosity`, `set_partial_melt`) store their models on the layer's
material, so a later `set_eos` must carry them over to the new material unless it brings its own, whether the layer is
standalone or already part of a solved world.
"""
import math

import numpy as np

from TidalPy.Material_x.eos import make_material_eos
from TidalPy.partial_melt_x import make_partial_melt
from TidalPy.rheology_x import make_rheology
from TidalPy.structures_x.layers import SolidLiquidLayer
from TidalPy.structures_x.worlds import LayeredWorld
from TidalPy.Tides_x.classes import make_tide
from TidalPy.viscosity_x import make_viscosity

RADIUS = 1.8e6
MASS = 8.9e22
DENSITY = MASS / ((4.0 / 3.0) * math.pi * RADIUS ** 3)


def _material(shear_modulus=6.0e10, shear_viscosity=None):
    config = {"reference_density_kg_m3": DENSITY, "shear_modulus_static_pa": shear_modulus}
    if shear_viscosity is not None:
        config["shear_viscosity"] = {"model": "constant", "reference_viscosity_pas": shear_viscosity}
    return make_material_eos("constant", config)


def _layer():
    layer = SolidLiquidLayer("mantle", 0, 0.0, RADIUS, 0.0)
    layer.set_eos(_material())
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e18}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e20}))
    layer.set_partial_melt(make_partial_melt("henning"))
    layer.set_shear_rheology(make_rheology("maxwell"))
    return layer


def test_replacing_the_material_keeps_the_layer_models():
    layer = _layer()
    layer.set_eos(_material(shear_modulus=2.0e10))
    assert layer.shear_viscosity_set
    assert layer.bulk_viscosity_set
    assert layer.partial_melt_set


def test_a_new_material_keeps_its_own_viscosity():
    layer = _layer()
    layer.set_eos(_material(shear_modulus=2.0e10, shear_viscosity=3.0e15))
    assert layer.shear_viscosity_set
    world = LayeredWorld("world", RADIUS, MASS)
    world.add_layer(layer)
    world.solve_eos()
    assert np.isclose(world.get_shear_viscosity(0.5 * RADIUS), 3.0e15, rtol=1.0e-12)


def test_material_swap_in_a_solved_world_keeps_its_tides():
    world = LayeredWorld("world", RADIUS, MASS)
    world.add_layer(_layer())
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0,
                          love_method="homogeneous")
    world.solve_eos()
    orbit = (2.0e-5, 2.0e-5, 0.01, 0.0, 4.2e8, 1.9e27)
    world.calc_tides(*orbit)
    stiff_heating = world.get_tidal_heating()

    world.get_layer(0).set_eos(_material(shear_modulus=2.0e10))
    world.solve_eos()
    world.calc_tides(*orbit)
    soft_heating = world.get_tidal_heating()
    assert np.isfinite(soft_heating) and soft_heating > 0.0
    assert not np.isclose(soft_heating, stiff_heating, rtol=1.0e-3)
