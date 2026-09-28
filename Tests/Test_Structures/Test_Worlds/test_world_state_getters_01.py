"""Bundled profile getters (get_state, get_static_viscoelastics) match the single-quantity getters exactly."""
import numpy as np
import pytest

from TidalPy.Structures import build_world

_KEYS = ("density", "gravity", "pressure", "shear_modulus", "shear_viscosity", "bulk_modulus", "bulk_viscosity",
         "melt_fraction")


@pytest.fixture(scope="module")
def io():
    world = build_world("io")
    world.solve_eos()
    return world


def _single(source, key, radius):
    return getattr(source, f"get_{key}")(radius)


@pytest.mark.parametrize("which", ["world", "layer"])
def test_get_state_matches_the_single_getters(io, which):
    """get_state matches each single getter for a 2D radius array and a scalar radius."""
    source = io if which == "world" else io.mantle
    radii = np.linspace(source.radius_inner if which == "layer" else 0.0, io.mantle.radius_outer, 24).reshape(4, 6)
    state = source.get_state(radii)
    assert tuple(state) == _KEYS
    for key in _KEYS:
        assert state[key].shape == radii.shape
        np.testing.assert_array_equal(state[key], _single(source, key, radii), err_msg=key)
    scalar = source.get_state(float(radii[1, 2]))
    for key in _KEYS:
        assert scalar[key] == _single(source, key, float(radii[1, 2])) or (
            np.isnan(scalar[key]) and np.isnan(_single(source, key, float(radii[1, 2])))), key


def test_static_viscoelastics_order(io):
    """get_static_viscoelastics returns shear, shear viscosity, bulk, bulk viscosity in that order."""
    radii = np.linspace(0.5 * io.radius, io.radius, 7)
    shear, shear_visc, bulk, bulk_visc = io.get_static_viscoelastics(radii)
    np.testing.assert_array_equal(shear, io.get_shear_modulus(radii))
    np.testing.assert_array_equal(shear_visc, io.get_shear_viscosity(radii))
    np.testing.assert_array_equal(bulk, io.get_bulk_modulus(radii))
    np.testing.assert_array_equal(bulk_visc, io.get_bulk_viscosity(radii))
