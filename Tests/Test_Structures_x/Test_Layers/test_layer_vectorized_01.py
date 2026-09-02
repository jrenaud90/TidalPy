"""Layer-level vectorized getter surface (float-or-ndarray radius dispatch).

Layers expose the same radius-query surface as the world: every real-valued ``get_*`` profile getter,
the ``get_static_viscoelastics`` / ``get_state`` shorthand bundles, and the radius-resolved complex
moduli accept a float or an np.ndarray of radii (array in, same-shape array out). The one-argument
layer-constant ``calc_complex_*(frequency)`` form is unchanged.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.rheology_x.rheology import Maxwell, Elastic

_RADIUS = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_MASS = (4.0 / 3.0) * math.pi * _RADIUS ** 3 * _DENSITY
_FREQ = 1.0e-5

_REAL_GETTERS = (
    "get_density", "get_gravity", "get_pressure",
    "get_shear_modulus", "get_bulk_modulus", "get_shear_viscosity", "get_bulk_viscosity",
    "get_premelt_shear_modulus", "get_premelt_bulk_modulus",
    "get_premelt_shear_viscosity", "get_premelt_bulk_viscosity",
)


def _solved_world():
    """A homogeneous Maxwell world with its EOS solved (populates the layer profiles)."""
    world = LayeredWorld("world", _RADIUS, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _RADIUS, _MASS,
                         shear_modulus_static_pa=_SHEAR, bulk_modulus_static_pa=_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.solve_eos(G_to_use=G)
    return world


def _standalone_layer():
    """A directly-populated layer (no world) for the base EOS-profile getters."""
    layer = PhysicsLayer("mantle", 0, 0.0, _RADIUS, _MASS,
                         shear_modulus_static_pa=_SHEAR, bulk_modulus_static_pa=_BULK)
    radius = np.linspace(0.0, _RADIUS, 11)
    density = np.full(11, _DENSITY)
    gravity = np.linspace(0.0, 9.0, 11)
    pressure = np.linspace(1.0e11, 0.0, 11)
    layer.update_eos_data(radius, density, gravity, pressure)
    return layer


# =====================================================================================================================
# Real-valued getters: scalar vs array agreement, shape preservation
# =====================================================================================================================
@pytest.mark.parametrize("getter_name", ("get_density", "get_gravity", "get_pressure"))
def test_standalone_layer_array_matches_scalar(getter_name):
    layer = _standalone_layer()
    getter = getattr(layer, getter_name)
    radii = np.linspace(0.1 * _RADIUS, 0.9 * _RADIUS, 7)
    array_result = getter(radii)
    assert isinstance(array_result, np.ndarray)
    assert array_result.shape == radii.shape
    for i, radius in enumerate(radii):
        assert math.isclose(array_result[i], getter(float(radius)), rel_tol=1e-14)


@pytest.mark.parametrize("getter_name", _REAL_GETTERS)
def test_solved_layer_array_matches_scalar(getter_name):
    world = _solved_world()
    layer = world.mantle
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


def test_layer_array_matches_world_surface():
    """Inside the layer, the layer getters agree with the world's (which delegate by radius)."""
    world = _solved_world()
    layer = world.mantle
    radii = np.linspace(0.2 * _RADIUS, 0.8 * _RADIUS, 5)
    for getter_name in ("get_density", "get_gravity", "get_pressure",
                        "get_shear_modulus", "get_shear_viscosity"):
        layer_vals = getattr(layer, getter_name)(radii)
        world_vals = getattr(world, getter_name)(radii)
        np.testing.assert_allclose(layer_vals, world_vals, rtol=1e-14)


def test_shape_and_noncontiguous_input():
    layer = _standalone_layer()
    radii_2d = np.linspace(0.1 * _RADIUS, 0.9 * _RADIUS, 12).reshape(3, 4)
    result_2d = layer.get_density(radii_2d)
    assert result_2d.shape == (3, 4)
    strided = np.linspace(0.1 * _RADIUS, 0.9 * _RADIUS, 10)[::2]
    result_strided = layer.get_density(strided)
    assert result_strided.shape == (5,)
    for i, radius in enumerate(strided):
        assert math.isclose(result_strided[i], layer.get_density(float(radius)), rel_tol=1e-14)


# =====================================================================================================================
# Shorthand bundles
# =====================================================================================================================
def test_get_static_viscoelastics_bundle():
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
    world = _solved_world()
    layer = world.mantle
    state = layer.get_state(0.5 * _RADIUS)
    expected_keys = {"density", "gravity", "pressure", "shear_modulus", "shear_viscosity",
                     "bulk_modulus", "bulk_viscosity"}
    assert set(state.keys()) == expected_keys
    assert math.isclose(state["density"], _DENSITY, rel_tol=1e-9)
    radii = np.linspace(0.2 * _RADIUS, 0.8 * _RADIUS, 3)
    state_arrays = layer.get_state(radii)
    for key in expected_keys:
        assert state_arrays[key].shape == radii.shape


# =====================================================================================================================
# Complex moduli: layer-constant (1-arg) and radius-resolved (2-arg) forms
# =====================================================================================================================
def test_layer_constant_complex_form_unchanged():
    world = _solved_world()
    layer = world.mantle
    mu = layer.calc_complex_shear_modulus(_FREQ)
    bulk = layer.calc_complex_bulk_modulus(_FREQ)
    assert isinstance(mu, complex)
    assert mu.imag != 0.0        # Maxwell dissipates
    assert bulk.imag == 0.0      # Elastic does not
    assert math.isclose(bulk.real, _BULK, rel_tol=1e-12)


def test_radius_resolved_complex_scalar_and_array():
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


def test_radius_resolved_complex_matches_world():
    """The layer's radius-resolved complex moduli equal the world's at the same radii."""
    world = _solved_world()
    layer = world.mantle
    radii = np.linspace(0.2 * _RADIUS, 0.8 * _RADIUS, 5)
    np.testing.assert_allclose(
        layer.calc_complex_shear_modulus(radii, _FREQ),
        world.calc_complex_shear_modulus(radii, _FREQ), rtol=1e-14)
    np.testing.assert_allclose(
        layer.calc_complex_bulk_modulus(radii, _FREQ),
        world.calc_complex_bulk_modulus(radii, _FREQ), rtol=1e-14)
