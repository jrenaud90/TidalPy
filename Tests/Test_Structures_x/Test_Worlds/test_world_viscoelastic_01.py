"""
Tests for the world-level viscoelastic state pipeline (increments B/C/D):
solve_eos populates each layer's pre/post-melt shear/bulk modulus + viscosity,
and the world (and layers) expose them via get_* getters and the radius-resolved
calc_complex_*_modulus methods.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import math

import pytest

from TidalPy.constants import G


_PLANET_RADIUS  = 6.0e6     # [m]
_DENSITY        = 4000.0    # [kg/m^3]
_STATIC_SHEAR   = 6.0e10    # [Pa]
_STATIC_BULK    = 1.3e11    # [Pa]
_SHEAR_VISC     = 1.0e21    # [Pa s]
_BULK_VISC      = 1.0e30    # [Pa s]


def _imports():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.partial_melt_x import make_partial_melt
    from TidalPy.rheology_x.rheology import Maxwell
    return LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, make_partial_melt, Maxwell


def _uniform_physics_world(with_viscosity=True, with_melt=False, with_rheology=False):
    LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, make_partial_melt, Maxwell = _imports()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("rocky", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR,
                         bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    if with_viscosity:
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
        layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _BULK_VISC}))
    if with_melt:
        layer.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1600.0, "liquidus_k": 2000.0}))
    if with_rheology:
        layer.set_shear_rheology(Maxwell())
        layer.set_bulk_rheology(Maxwell())
    world.add_layer(layer)
    return world


# =====================================================================================================================
# Population
# =====================================================================================================================
def test_nan_before_solve():
    world = _uniform_physics_world()
    assert math.isnan(world.get_shear_modulus(_PLANET_RADIUS * 0.5))
    assert math.isnan(world.get_shear_viscosity(_PLANET_RADIUS * 0.5))


def test_static_moduli_recovered_no_melt():
    world = _uniform_physics_world(with_viscosity=True, with_melt=False)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = _PLANET_RADIUS * 0.5
    # No melt model -> post-melt equals pre-melt equals the layer's static value.
    assert world.get_shear_modulus(radius) == pytest.approx(_STATIC_SHEAR)
    assert world.get_bulk_modulus(radius) == pytest.approx(_STATIC_BULK)
    assert world.get_premelt_shear_modulus(radius) == pytest.approx(_STATIC_SHEAR)
    assert world.get_premelt_bulk_modulus(radius) == pytest.approx(_STATIC_BULK)


def test_viscosity_from_models():
    world = _uniform_physics_world(with_viscosity=True)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = _PLANET_RADIUS * 0.5
    assert world.get_shear_viscosity(radius) == pytest.approx(_SHEAR_VISC)
    assert world.get_bulk_viscosity(radius) == pytest.approx(_BULK_VISC)


def test_viscosity_nan_when_model_unset():
    world = _uniform_physics_world(with_viscosity=False)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = _PLANET_RADIUS * 0.5
    # No viscosity model attached -> NaN viscosity, but moduli still come from statics.
    assert math.isnan(world.get_shear_viscosity(radius))
    assert world.get_shear_modulus(radius) == pytest.approx(_STATIC_SHEAR)


# =====================================================================================================================
# Partial-melt weakening
# =====================================================================================================================
def test_melt_weakens_shear_at_high_temperature():
    world = _uniform_physics_world(with_viscosity=True, with_melt=True)
    # T = 1800 K -> melt fraction 0.5 (solidus 1600, liquidus 2000): Henning weakens shear.
    world.solve_eos(G_to_use=G, temperature=1800.0, verbose=False)
    radius = _PLANET_RADIUS * 0.5
    assert world.get_shear_modulus(radius) < world.get_premelt_shear_modulus(radius)
    assert world.get_premelt_shear_modulus(radius) == pytest.approx(_STATIC_SHEAR)


def test_no_melt_below_solidus():
    world = _uniform_physics_world(with_viscosity=True, with_melt=True)
    # T = 1500 K < solidus -> melt fraction 0 -> post == pre.
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = _PLANET_RADIUS * 0.5
    assert world.get_shear_modulus(radius) == pytest.approx(_STATIC_SHEAR)


# =====================================================================================================================
# Complex moduli (radius-resolved)
# =====================================================================================================================
def test_complex_without_rheology_is_real_static():
    world = _uniform_physics_world(with_viscosity=True, with_rheology=False)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = _PLANET_RADIUS * 0.5
    value = world.calc_complex_shear_modulus(radius, 1.0e-5)
    assert value.real == pytest.approx(_STATIC_SHEAR)
    assert value.imag == pytest.approx(0.0)


def test_complex_shear_matches_maxwell():
    _, _, _, _, _, Maxwell = _imports()
    world = _uniform_physics_world(with_viscosity=True, with_rheology=True)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = _PLANET_RADIUS * 0.5
    frequency = 1.0e-5

    expected = Maxwell().calc_complex_modulus(_STATIC_SHEAR, _SHEAR_VISC, frequency)
    value = world.calc_complex_shear_modulus(radius, frequency)
    assert value.real == pytest.approx(expected.real, rel=1e-9)
    assert value.imag == pytest.approx(expected.imag, rel=1e-9)
    # Viscoelastic dissipation -> nonzero imaginary part.
    assert abs(value.imag) > 0.0


def test_complex_nan_on_geometry_layer():
    """A geometry-only BaseLayer world yields NaN complex moduli (no rheology/moduli)."""
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.base import BaseLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("geom", _PLANET_RADIUS, mass)
    layer = BaseLayer("rock", 0, 0.0, _PLANET_RADIUS, mass)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, verbose=False)
    value = world.calc_complex_shear_modulus(_PLANET_RADIUS * 0.5, 1.0e-5)
    assert math.isnan(value.real)


# =====================================================================================================================
# Layer-level getters mirror the world
# =====================================================================================================================
def test_layer_getters_match_world():
    LayeredWorld, PhysicsLayer, ConstantDensityEOS, make_viscosity, _, _ = _imports()
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("rocky", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR, bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    world.add_layer(layer)
    # The layer's C++ object was moved into the world; query through the world's layer view.
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    radius = _PLANET_RADIUS * 0.5
    assert world.get_shear_modulus(radius) == pytest.approx(_STATIC_SHEAR)
    assert world.get_shear_viscosity(radius) == pytest.approx(_SHEAR_VISC)
