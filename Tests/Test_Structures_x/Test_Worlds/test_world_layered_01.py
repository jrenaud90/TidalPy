"""
Tests for TidalPy.structures_x.worlds — LayeredWorld and GasGiantWorld (Phase 8b).

Covers layer ownership (add/aggregate/validate), geometry-continuity errors,
internal radiogenic heating, recursive binary round-trip (world + layers + each
layer's physics sub-models), config dict, and isinstance checks.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import math
import os
import tempfile

import pytest


# =====================================================================================================================
# Helpers
# =====================================================================================================================
def _import_layered():
    from TidalPy.structures_x.worlds import layered as _mod
    return _mod


def _import_gasgiant():
    from TidalPy.structures_x.worlds import gasgiant as _mod
    return _mod


def _import_layers():
    from TidalPy.structures_x.layers import base, physics, solidliquid, gas
    return base, physics, solidliquid, gas


def _import_radiogenics():
    from TidalPy.radiogenics_x import radiogenics as _mod
    return _mod


# Reference geometry (Earth-like two-layer, MKS)
_R_CMB   = 3.485e6
_R_SURF  = 6.371e6
_M_CORE  = 1.932e24
_M_MANT  = 4.040e24
_M_TOT   = _M_CORE + _M_MANT


def _make_two_layer_world():
    mod = _import_layered()
    base, physics, solidliquid, gas = _import_layers()
    world = mod.LayeredWorld("Earth", _R_SURF, _M_TOT, world_type="terrestrial")
    core = solidliquid.SolidLiquidLayer("core", 0, 0.0, _R_CMB, _M_CORE,
                                        material_name="iron")
    mantle = solidliquid.SolidLiquidLayer("mantle", 1, _R_CMB, _R_SURF, _M_MANT,
                                          material_name="perovskite")
    world.add_layer(core)
    world.add_layer(mantle)
    return world


# =====================================================================================================================
# Layer ownership
# =====================================================================================================================
def test_add_layers_and_count():
    world = _make_two_layer_world()
    assert world.num_layers == 2
    assert world.validate_layers() is True


def test_add_layer_consumes_wrapper():
    mod = _import_layered()
    _, _, solidliquid, _ = _import_layers()
    world = mod.LayeredWorld("W", _R_SURF, _M_TOT)
    layer = solidliquid.SolidLiquidLayer("core", 0, 0.0, _R_CMB, _M_CORE)
    world.add_layer(layer)
    with pytest.raises(ValueError):
        world.add_layer(layer)


def test_add_layer_discontinuity_raises():
    mod = _import_layered()
    _, _, solidliquid, _ = _import_layers()
    world = mod.LayeredWorld("W", _R_SURF, _M_TOT)
    # First layer must start at radius 0; this one starts at R_CMB -> gap.
    bad = solidliquid.SolidLiquidLayer("mantle", 0, _R_CMB, _R_SURF, _M_MANT)
    with pytest.raises(ValueError):
        world.add_layer(bad)
    # A rejected layer is NOT consumed and remains usable elsewhere.
    world2 = mod.LayeredWorld("W2", _R_SURF, _M_MANT)
    world2.add_layer(solidliquid.SolidLiquidLayer("core", 0, 0.0, _R_CMB, _M_CORE))
    world2.add_layer(bad)  # now continuous (inner R_CMB matches core's outer)
    assert world2.num_layers == 2


def test_calc_total_mass():
    world = _make_two_layer_world()
    assert world.calc_total_mass() == pytest.approx(_M_TOT, rel=1e-9)


def test_mixed_layer_types():
    mod = _import_layered()
    base, physics, solidliquid, gas = _import_layers()
    world = mod.LayeredWorld("Mixed", _R_SURF, _M_TOT)
    world.add_layer(base.BaseLayer("core", 0, 0.0, _R_CMB, _M_CORE))
    world.add_layer(physics.PhysicsLayer("mantle", 1, _R_CMB, _R_SURF, _M_MANT))
    assert world.num_layers == 2
    assert world.calc_total_mass() == pytest.approx(_M_TOT, rel=1e-9)


# =====================================================================================================================
# Internal radiogenic heating
# =====================================================================================================================
def test_internal_heating_zero_without_radiogenics():
    world = _make_two_layer_world()
    assert world.calc_internal_heating(0.0) == pytest.approx(0.0)


def test_internal_heating_with_radiogenics():
    mod   = _import_layered()
    _, _, solidliquid, _ = _import_layers()
    radio = _import_radiogenics()

    world = mod.LayeredWorld("Earth", _R_SURF, _M_TOT)
    core  = solidliquid.SolidLiquidLayer("core", 0, 0.0, _R_CMB, _M_CORE)
    mantle = solidliquid.SolidLiquidLayer("mantle", 1, _R_CMB, _R_SURF, _M_MANT)
    mantle.set_radiogenics(radio.FixedRadiogenics(fixed_heat_production_w_kg=1.0e-11))
    world.add_layer(core)
    world.add_layer(mantle)

    # Only the mantle has radiogenics: H = rate * mantle_mass.
    expected = 1.0e-11 * _M_MANT
    assert world.calc_internal_heating(0.0) == pytest.approx(expected, rel=1e-9)


# =====================================================================================================================
# Recursive binary round-trip (world + layers + sub-models)
# =====================================================================================================================
def test_layered_world_binary_roundtrip():
    mod   = _import_layered()
    _, _, solidliquid, gas = _import_layers()
    radio = _import_radiogenics()

    world = mod.LayeredWorld("Earth", _R_SURF, _M_TOT, world_type="terrestrial",
                             albedo=0.31, obliquity_rad=0.41)
    core  = solidliquid.SolidLiquidLayer("core", 0, 0.0, _R_CMB, _M_CORE,
                                         material_name="iron")
    mantle = solidliquid.SolidLiquidLayer("mantle", 1, _R_CMB, _R_SURF, _M_MANT,
                                          material_name="perovskite")
    mantle.set_radiogenics(radio.FixedRadiogenics(fixed_heat_production_w_kg=2.0e-11))
    world.add_layer(core)
    world.add_layer(mantle)
    heating_before = world.calc_internal_heating(0.0)

    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        world.save_binary(path)
        w2 = mod.LayeredWorld("placeholder", 1.0, 1.0)
        w2.load_binary(path)

        assert w2.name        == "Earth"
        assert w2.world_type  == "terrestrial"
        assert w2.albedo      == pytest.approx(0.31)
        assert w2.obliquity   == pytest.approx(0.41)
        assert w2.num_layers  == 2
        assert w2.calc_total_mass() == pytest.approx(_M_TOT, rel=1e-9)
        # The mantle's radiogenics model must survive the recursive round-trip.
        assert w2.calc_internal_heating(0.0) == pytest.approx(heating_before, rel=1e-12)
        # Per-layer config (geometry level) is restored in order.
        cfg = w2.get_config_dict()
        assert cfg["num_layers"] == 2
        assert cfg["layers"][0]["name"] == "core"
        assert cfg["layers"][1]["name"] == "mantle"
        assert cfg["layers"][1]["radius_outer_m"] == pytest.approx(_R_SURF)
    finally:
        os.unlink(path)


# =====================================================================================================================
# GasGiantWorld
# =====================================================================================================================
def test_gasgiant_construction_and_type():
    mod = _import_gasgiant()
    _, _, _, gas = _import_layers()
    gg = mod.GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
    assert gg.world_type == "gasgiant"
    gg.add_layer(gas.GasLayer("envelope", 0, 0.0, 7.0e7, 1.898e27))
    assert gg.num_layers == 1


def test_gasgiant_binary_roundtrip():
    mod = _import_gasgiant()
    _, _, _, gas = _import_layers()
    gg = mod.GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
    gg.add_layer(gas.GasLayer("envelope", 0, 0.0, 7.0e7, 1.898e27))
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name
    try:
        gg.save_binary(path)
        gg2 = mod.GasGiantWorld("placeholder", 1.0, 1.0)
        gg2.load_binary(path)
        assert gg2.name       == "Jupiter"
        assert gg2.world_type == "gasgiant"
        assert gg2.num_layers == 1
    finally:
        os.unlink(path)


# =====================================================================================================================
# isinstance
# =====================================================================================================================
def test_layered_is_base_world():
    from TidalPy.structures_x.worlds.base import BaseWorld
    assert isinstance(_make_two_layer_world(), BaseWorld)


def test_gasgiant_is_layered_world():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    mod = _import_gasgiant()
    gg = mod.GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
    assert isinstance(gg, LayeredWorld)
