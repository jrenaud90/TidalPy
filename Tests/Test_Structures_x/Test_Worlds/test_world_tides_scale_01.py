"""Tests for the per-layer ``tidal_scale_method`` (how the world's global tidal heating is
distributed to the layers in ``LayeredWorld.calc_tides``).

Methods:
  * ``user_provided``   : use the layer's ``tidal_scale`` field directly (default).
  * ``volume_fraction`` : layer volume / planet volume.
  * ``tidal_timescale`` : Maxwell-time bell curve (not yet wired; must raise loudly).

These use an analytic ``cpl`` tide model so no EOS solve is needed (the scale distribution
depends only on layer geometry and the model-supplied global heating).
"""
import math

import pytest

from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Tides_x.classes.tide import make_tide


_R         = 1.6e6
_HOST_MASS = 1.898e27
_SMA       = 4.2e8
_N         = 2.05e-5
_ECC       = 0.0041


def _cpl_two_layer(core_method="user_provided", mantle_method="user_provided",
                   core_scale=1.0, mantle_scale=1.0):
    r_core = 0.5 * _R   # core volume fraction = (0.5)^3 = 0.125
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0
    world = LayeredWorld("scaled", _R, mass)
    core = PhysicsLayer("core", 0, 0.0, r_core, 0.0,
                        tidal_scale=core_scale, tidal_scale_method=core_method)
    mantle = PhysicsLayer("mantle", 1, r_core, _R, 0.0,
                          tidal_scale=mantle_scale, tidal_scale_method=mantle_method)
    world.add_layer(core)
    world.add_layer(mantle)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=2, obliquity_truncation=0)
    return world


def _solve(world):
    world.calc_tides(orbital_frequency=_N, spin_frequency=_N, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST_MASS)


# =====================================================================================================================
# tidal_scale_method property
# =====================================================================================================================
def test_default_scale_method_is_user_provided():
    layer = PhysicsLayer("mantle", 0, 0.0, _R, 0.0)
    assert layer.tidal_scale_method == "user_provided_scale"


def test_scale_method_round_trips_via_aliases():
    layer = PhysicsLayer("mantle", 0, 0.0, _R, 0.0, tidal_scale_method="volume_fraction")
    assert layer.tidal_scale_method == "volume_fraction_scale"
    layer.tidal_scale_method = "timescale"
    assert layer.tidal_scale_method == "tidal_timescale_scale"


def test_unknown_scale_method_raises():
    with pytest.raises(ValueError):
        PhysicsLayer("mantle", 0, 0.0, _R, 0.0, tidal_scale_method="nonsense")


# =====================================================================================================================
# user_provided (default) distribution
# =====================================================================================================================
def test_user_provided_uses_tidal_scale_field():
    world = _cpl_two_layer(core_scale=0.2, mantle_scale=0.7)
    _solve(world)
    total = world.get_tidal_heating()
    assert math.isclose(world.get_layer_tidal_heating(0), total * 0.2, rel_tol=1.0e-12)
    assert math.isclose(world.get_layer_tidal_heating(1), total * 0.7, rel_tol=1.0e-12)


# =====================================================================================================================
# volume_fraction distribution
# =====================================================================================================================
def test_volume_fraction_distributes_by_volume():
    world = _cpl_two_layer(core_method="volume_fraction", mantle_method="volume_fraction")
    _solve(world)
    total = world.get_tidal_heating()
    # core fills the inner half-radius -> volume fraction (0.5)^3 = 0.125; mantle gets the rest.
    assert math.isclose(world.get_layer_tidal_heating(0), total * 0.125, rel_tol=1.0e-9)
    assert math.isclose(world.get_layer_tidal_heating(1), total * 0.875, rel_tol=1.0e-9)
    # Volume-fraction shares sum to the whole-body heating.
    assert math.isclose(world.get_layer_tidal_heating(0) + world.get_layer_tidal_heating(1),
                        total, rel_tol=1.0e-9)


def test_methods_can_differ_per_layer():
    """A volume_fraction core alongside a user_provided mantle each use their own rule."""
    world = _cpl_two_layer(core_method="volume_fraction",
                           mantle_method="user_provided", mantle_scale=0.5)
    _solve(world)
    total = world.get_tidal_heating()
    assert math.isclose(world.get_layer_tidal_heating(0), total * 0.125, rel_tol=1.0e-9)
    assert math.isclose(world.get_layer_tidal_heating(1), total * 0.5, rel_tol=1.0e-12)


# =====================================================================================================================
# tidal_timescale is not yet wired
# =====================================================================================================================
def test_tidal_timescale_method_raises_pending():
    world = _cpl_two_layer(core_method="tidal_timescale", mantle_method="user_provided")
    with pytest.raises(RuntimeError):
        _solve(world)
