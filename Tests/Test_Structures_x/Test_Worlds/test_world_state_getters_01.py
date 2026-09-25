"""The bundled profile getters read one evaluation of the solved state per radius and agree with the single ones.

``get_state`` and ``get_static_viscoelastics`` on a world and on a layer fill every quantity from one dense
evaluation, so they must return exactly what the one-quantity getters return, for a scalar radius and for an array
of any shape.
"""
import numpy as np
import pytest

from TidalPy.structures_x import build_world

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
    radii = np.linspace(0.5 * io.radius, io.radius, 7)
    shear, shear_visc, bulk, bulk_visc = io.get_static_viscoelastics(radii)
    np.testing.assert_array_equal(shear, io.get_shear_modulus(radii))
    np.testing.assert_array_equal(shear_visc, io.get_shear_viscosity(radii))
    np.testing.assert_array_equal(bulk, io.get_bulk_modulus(radii))
    np.testing.assert_array_equal(bulk_visc, io.get_bulk_viscosity(radii))
