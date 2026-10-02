"""Layer radius getters and radius-resolved complex moduli take a float or an array (same-shape array out)."""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.layers import Layer
from TidalPy.Material import Material
from TidalPy.Rheology.rheology import Maxwell, Elastic

from shared_materials import constant_solid

_RADIUS = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_MASS = (4.0 / 3.0) * math.pi * _RADIUS ** 3 * _DENSITY
_FREQ = 1.0e-5

_PROFILE_GETTERS = ("get_density", "get_gravity", "get_pressure")
_REAL_GETTERS = _PROFILE_GETTERS + (
    "get_shear_modulus", "get_bulk_modulus", "get_shear_viscosity", "get_bulk_viscosity", "get_melt_fraction")


def _material():
    return constant_solid(
        _DENSITY, bulk_modulus=_BULK, shear_modulus=_SHEAR, shear_viscosity=_VISC, bulk_viscosity=_VISC)


def _solved_world():
    """A homogeneous Maxwell world with its EOS solved."""
    world = BaseWorld("world", _RADIUS, _MASS)
    world.add_layer(Layer("mantle", 0, 0.0, _RADIUS, _MASS, _material(), shear_rheology=Maxwell(),
                          bulk_rheology=Elastic()))
    world.solve_eos(G_to_use=G)
    return world


def _standalone_layer():
    """A layer with a directly populated EOS profile and no world."""
    layer = Layer("mantle", 0, 0.0, _RADIUS, _MASS)
    layer.update_eos_data(
        np.linspace(0.0, _RADIUS, 11),
        np.full(11, _DENSITY),
        np.linspace(0.0, 9.0, 11),
        np.linspace(1.0e11, 0.0, 11))
    return layer


@pytest.mark.parametrize("source, getter_name", (
    [pytest.param("standalone", name, id=f"standalone-{name}") for name in _PROFILE_GETTERS]
    + [pytest.param("solved", name, id=f"solved-{name}") for name in _REAL_GETTERS]))
def test_array_matches_scalar(source, getter_name):
    """An array of radii gives the same finite values as scalar calls, in the input's shape."""
    # The world owns a solved layer's view, so it must stay alive for the whole test.
    world = _solved_world() if source == "solved" else None
    layer = world.mantle if world is not None else _standalone_layer()
    getter = getattr(layer, getter_name)
    radii = np.linspace(0.1 * _RADIUS, 0.9 * _RADIUS, 7)
    array_result = getter(radii)
    assert isinstance(array_result, np.ndarray)
    assert array_result.shape == radii.shape
    for i, radius in enumerate(radii):
        scalar_result = getter(float(radius))
        assert isinstance(scalar_result, float)
        assert math.isclose(array_result[i], scalar_result, rel_tol=1e-14)
        assert np.isfinite(array_result[i])


@pytest.mark.parametrize(
    "getter_name",
    ("get_density", "get_gravity", "get_pressure", "get_shear_modulus", "get_shear_viscosity"))
def test_layer_array_matches_world(getter_name):
    """Inside the layer, the layer getters agree with the world's."""
    world = _solved_world()
    radii = np.linspace(0.2 * _RADIUS, 0.8 * _RADIUS, 5)
    np.testing.assert_allclose(
        getattr(world.mantle, getter_name)(radii), getattr(world, getter_name)(radii), rtol=1e-14)


def test_shape_and_noncontiguous_input():
    """A 2D array keeps its shape and a strided view gives the scalar values."""
    layer = _standalone_layer()
    radii_2d = np.linspace(0.1 * _RADIUS, 0.9 * _RADIUS, 12).reshape(3, 4)
    assert layer.get_density(radii_2d).shape == (3, 4)
    strided = np.linspace(0.1 * _RADIUS, 0.9 * _RADIUS, 10)[::2]
    result_strided = layer.get_density(strided)
    assert result_strided.shape == (5,)
    for i, radius in enumerate(strided):
        assert math.isclose(result_strided[i], layer.get_density(float(radius)), rel_tol=1e-14)


def test_get_static_viscoelastics_bundle():
    """The bundle matches the individual getters, and a scalar radius gives floats."""
    world = _solved_world()
    layer = world.mantle
    radii = np.linspace(0.2 * _RADIUS, 0.8 * _RADIUS, 4)
    shear_mod, shear_visc, bulk_mod, bulk_visc = layer.get_static_viscoelastics(radii)
    np.testing.assert_allclose(shear_mod, layer.get_shear_modulus(radii), rtol=1e-14)
    np.testing.assert_allclose(shear_visc, layer.get_shear_viscosity(radii), rtol=1e-14)
    np.testing.assert_allclose(bulk_mod, layer.get_bulk_modulus(radii), rtol=1e-14)
    np.testing.assert_allclose(bulk_visc, layer.get_bulk_viscosity(radii), rtol=1e-14)
    scalar_bundle = layer.get_static_viscoelastics(0.5 * _RADIUS)
    assert all(isinstance(value, float) for value in scalar_bundle)


def test_get_state_bundle():
    """get_state returns every profile quantity, as scalars or radius-shaped arrays."""
    world = _solved_world()
    layer = world.mantle
    state = layer.get_state(0.5 * _RADIUS)
    expected_keys = {"density", "gravity", "pressure", "shear_modulus", "shear_viscosity",
                     "bulk_modulus", "bulk_viscosity", "melt_fraction"}
    assert set(state.keys()) == expected_keys
    assert math.isclose(state["density"], _DENSITY, rel_tol=1e-9)
    radii = np.linspace(0.2 * _RADIUS, 0.8 * _RADIUS, 3)
    state_arrays = layer.get_state(radii)
    for key in expected_keys:
        assert state_arrays[key].shape == radii.shape


def test_layer_constant_complex_form_unchanged():
    """The one-argument calc_complex_* form returns the layer's material at zero pressure, as a complex."""
    world = _solved_world()
    layer = world.mantle
    mu = layer.calc_complex_shear_modulus(_FREQ)
    bulk = layer.calc_complex_bulk_modulus(_FREQ)
    assert isinstance(mu, complex)
    assert mu.imag != 0.0        # Maxwell dissipates
    assert bulk.imag == 0.0      # Elastic does not
    assert math.isclose(bulk.real, _BULK, rel_tol=1e-12)


def test_radius_resolved_complex_scalar_and_array():
    """The two-argument form gives complex arrays matching scalar calls; the elastic bulk stays real."""
    world = _solved_world()
    layer = world.mantle
    radii = np.linspace(0.2 * _RADIUS, 0.8 * _RADIUS, 6)

    mu_array = layer.calc_complex_shear_modulus(radii, _FREQ)
    assert isinstance(mu_array, np.ndarray)
    assert mu_array.dtype == np.complex128
    assert mu_array.shape == radii.shape
    for i, radius in enumerate(radii):
        mu_scalar = layer.calc_complex_shear_modulus(float(radius), _FREQ)
        assert isinstance(mu_scalar, complex)
        assert math.isclose(mu_array[i].real, mu_scalar.real, rel_tol=1e-14)
        assert math.isclose(mu_array[i].imag, mu_scalar.imag, rel_tol=1e-14)

    bulk_array = layer.calc_complex_bulk_modulus(radii, _FREQ)
    assert bulk_array.shape == radii.shape
    np.testing.assert_allclose(bulk_array.imag, 0.0, atol=1e-30)


@pytest.mark.parametrize("method", ("calc_complex_shear_modulus", "calc_complex_bulk_modulus"))
def test_radius_resolved_complex_matches_world(method):
    """The layer's radius-resolved complex moduli equal the world's."""
    world = _solved_world()
    radii = np.linspace(0.2 * _RADIUS, 0.8 * _RADIUS, 5)
    np.testing.assert_allclose(
        getattr(world.mantle, method)(radii, _FREQ), getattr(world, method)(radii, _FREQ), rtol=1e-14)
