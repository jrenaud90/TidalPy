"""Vectorized world profile getters: a scalar radius gives a scalar, an array gives a same-shape array matching it."""

import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material import Material
from TidalPy.Rheology.rheology import Maxwell
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld

from shared_materials import constant_solid


_PLANET_RADIUS = 6.0e6
_DENSITY       = 4000.0
_STATIC_SHEAR  = 6.0e10
_STATIC_BULK   = 1.3e11
_SHEAR_VISC    = 1.0e21
_FREQUENCY     = 1.0e-5


def _material():
    return constant_solid(_DENSITY, bulk_modulus=_STATIC_BULK, shear_modulus=_STATIC_SHEAR, shear_viscosity=_SHEAR_VISC)


def _solved_world():
    mass = (4.0 / 3.0) * math.pi * _PLANET_RADIUS ** 3 * _DENSITY
    world = BaseWorld("rocky", _PLANET_RADIUS, mass)
    layer = Layer(
        "mantle",
        0,
        0.0,
        _PLANET_RADIUS,
        mass,
        _material(),
        shear_rheology=Maxwell(),
    )
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, temperature=1500.0, verbose=False)
    return world


_REAL_GETTERS = [
    "get_density", "get_gravity", "get_pressure",
    "get_shear_modulus", "get_bulk_modulus", "get_shear_viscosity", "get_bulk_viscosity", "get_melt_fraction",
]
_COMPLEX_GETTERS = ["calc_complex_shear_modulus", "calc_complex_bulk_modulus"]

# (getter name, extra arguments after the radius, scalar type, array dtype, number of radii in the array test)
_GETTER_CASES = (
    [pytest.param(name, (), float, np.float64, 17, id=name) for name in _REAL_GETTERS]
    + [pytest.param(name, (_FREQUENCY,), complex, np.complex128, 13, id=name) for name in _COMPLEX_GETTERS]
)


@pytest.mark.parametrize("name, extra_args, scalar_type, array_dtype, num_radii", _GETTER_CASES)
def test_scalar_radius_gives_a_scalar(
        name,
        extra_args,
        scalar_type,
        array_dtype,
        num_radii,
):
    """A scalar radius returns a Python float (complex for the complex moduli)."""
    world = _solved_world()
    value = getattr(world, name)(_PLANET_RADIUS * 0.5, *extra_args)
    assert isinstance(value, scalar_type)


@pytest.mark.parametrize("name, extra_args, scalar_type, array_dtype, num_radii", _GETTER_CASES)
def test_array_matches_scalar(
        name,
        extra_args,
        scalar_type,
        array_dtype,
        num_radii,
):
    """An array of radii returns a same-shape array matching the elementwise scalar calls."""
    world = _solved_world()
    getter = getattr(world, name)
    radii = np.linspace(0.0, _PLANET_RADIUS, num_radii)
    out = getter(radii, *extra_args)
    assert isinstance(out, np.ndarray)
    assert out.shape == radii.shape
    assert out.dtype == array_dtype
    for i, radius in enumerate(radii):
        scalar = getter(float(radius), *extra_args)
        if math.isnan(scalar.real):
            assert math.isnan(out[i].real)
        else:
            assert out[i] == pytest.approx(scalar)


def test_array_shape_preserved_2d():
    """A 2D radius array keeps its shape."""
    world = _solved_world()
    radii = np.linspace(0.0, _PLANET_RADIUS, 12).reshape(3, 4)
    out = world.get_density(radii)
    assert out.shape == (3, 4)
    assert out[1, 2] == pytest.approx(world.get_density(float(radii[1, 2])))
