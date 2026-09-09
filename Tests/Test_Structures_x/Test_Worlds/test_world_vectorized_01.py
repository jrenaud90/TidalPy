"""
Tests for the vectorized world profile getters: every get_*/calc_complex_* accepts
a scalar radius (returns a float/complex) OR a NumPy array of radii (returns an
array of the same shape), and the array result matches the elementwise scalar.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import math

import numpy as np
import pytest

from TidalPy.constants import G


_PLANET_RADIUS = 6.0e6
_DENSITY       = 4000.0
_STATIC_SHEAR  = 6.0e10
_STATIC_BULK   = 1.3e11
_SHEAR_VISC    = 1.0e21


def _solved_world():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell

    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = LayeredWorld("rocky", _PLANET_RADIUS, mass)
    layer = PhysicsLayer("mantle", 0, 0.0, _PLANET_RADIUS, mass,
                         shear_modulus_static_pa=_STATIC_SHEAR, bulk_modulus_static_pa=_STATIC_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _SHEAR_VISC}))
    layer.set_shear_rheology(Maxwell())
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    return world


_REAL_GETTERS = [
    "get_density", "get_gravity", "get_pressure",
    "get_shear_modulus", "get_bulk_modulus", "get_shear_viscosity", "get_bulk_viscosity",
    "get_premelt_shear_modulus", "get_premelt_bulk_modulus",
    "get_premelt_shear_viscosity", "get_premelt_bulk_viscosity",
]


@pytest.mark.parametrize("name", _REAL_GETTERS)
def test_real_scalar_is_float(name):
    world = _solved_world()
    value = getattr(world, name)(_PLANET_RADIUS * 0.5)
    assert isinstance(value, float)


@pytest.mark.parametrize("name", _REAL_GETTERS)
def test_real_array_matches_scalar(name):
    world = _solved_world()
    getter = getattr(world, name)
    radii = np.linspace(0.0, _PLANET_RADIUS, 17)
    out = getter(radii)
    assert isinstance(out, np.ndarray)
    assert out.shape == radii.shape
    assert out.dtype == np.float64
    for i, radius in enumerate(radii):
        scalar = getter(float(radius))
        if math.isnan(scalar):
            assert math.isnan(out[i])
        else:
            assert out[i] == pytest.approx(scalar)


def test_array_shape_preserved_2d():
    world = _solved_world()
    radii = np.linspace(0.0, _PLANET_RADIUS, 12).reshape(3, 4)
    out = world.get_density(radii)
    assert out.shape == (3, 4)
    assert out[1, 2] == pytest.approx(world.get_density(float(radii[1, 2])))


def test_complex_scalar_is_complex():
    world = _solved_world()
    value = world.calc_complex_shear_modulus(_PLANET_RADIUS * 0.5, 1.0e-5)
    assert isinstance(value, complex)


@pytest.mark.parametrize("name", ["calc_complex_shear_modulus", "calc_complex_bulk_modulus"])
def test_complex_array_matches_scalar(name):
    world = _solved_world()
    getter = getattr(world, name)
    radii = np.linspace(0.0, _PLANET_RADIUS, 13)
    out = getter(radii, 1.0e-5)
    assert isinstance(out, np.ndarray)
    assert out.dtype == np.complex128
    assert out.shape == radii.shape
    for i, radius in enumerate(radii):
        scalar = getter(float(radius), 1.0e-5)
        if math.isnan(scalar.real):
            assert math.isnan(out[i].real)
        else:
            assert out[i] == pytest.approx(scalar)
