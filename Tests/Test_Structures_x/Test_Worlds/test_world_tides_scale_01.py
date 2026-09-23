"""Tests for the per-layer ``tidal_scale_method`` (how the world's global tidal heating is
distributed to the layers in ``LayeredWorld.calc_tides``).

Methods:
  * ``user_provided``   : use the layer's ``tidal_scale`` field directly (default).
  * ``volume_fraction`` : layer volume / planet volume.
  * ``tidal_timescale`` : the layers using it share their combined volume fraction, split by volume times a
    Maxwell-time bell about the forcing period.

These use an analytic ``cpl`` tide model so no EOS solve is needed (the scale distribution
depends only on layer geometry and the model-supplied global heating).
"""
import math

import pytest

from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
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
# tidal_timescale (Maxwell-time bell curve about the orbital forcing period)
# =====================================================================================================================
def _cpl_timescale_world(layers, width=1.0):
    """A cpl world of concentric layers, each ``(outer radius fraction, shear modulus, shear viscosity, method)``."""
    mass = (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0
    world = LayeredWorld("ts", _R, mass)
    radius_inner = 0.0
    for index, (radius_fraction, shear_modulus, shear_viscosity, method) in enumerate(layers):
        layer = PhysicsLayer(f"layer_{index}", index, radius_inner, radius_fraction * _R, 0.0,
                             tidal_scale_method=method)
        layer.set_eos(ConstantDensityEOS(
            reference_density=4000.0, shear_modulus_static=shear_modulus, shear_viscosity_static=shear_viscosity))
        world.add_layer(layer)
        radius_inner = radius_fraction * _R
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=2, obliquity_truncation=0,
                          tidal_timescale_width_decades=width)
    return world


_FORCING_PERIOD = 2.0 * math.pi / _N
_MU = 6.0e10
_ETA_AT_PEAK = _FORCING_PERIOD * _MU   # tau = eta / mu = the forcing period


def test_tidal_timescale_single_layer_takes_its_volume_fraction():
    """A lone layer using the method gets the whole world's heating wherever its Maxwell time sits."""
    for viscosity in (_ETA_AT_PEAK, 10.0 * _ETA_AT_PEAK, 1.0e-3 * _ETA_AT_PEAK):
        world = _cpl_timescale_world([(1.0, _MU, viscosity, "tidal_timescale")])
        _solve(world)
        assert math.isclose(world.get_layer_tidal_heating(0), world.get_tidal_heating(), rel_tol=1.0e-12)


@pytest.mark.parametrize("width,decades", [(1.0, 1.0), (1.0, 2.0), (2.0, 1.0)])
def test_tidal_timescale_splits_by_volume_times_the_bell(width, decades):
    """Per unit volume, a layer ``decades`` off the forcing period gets exp(-(decades / width)^2 / 2) of one at it."""
    world = _cpl_timescale_world([
        (0.5, _MU, _ETA_AT_PEAK, "tidal_timescale"),
        (1.0, _MU, 10.0 ** decades * _ETA_AT_PEAK, "tidal_timescale")], width=width)
    _solve(world)
    total = world.get_tidal_heating()
    inner_volume, outer_volume = 0.125, 0.875   # fractions of the planet
    inner, outer = world.get_layer_tidal_heating(0) / total, world.get_layer_tidal_heating(1) / total
    weight = math.exp(-0.5 * (decades / width) ** 2)
    assert math.isclose((outer / outer_volume) / (inner / inner_volume), weight, rel_tol=1.0e-9)
    # Every layer uses the method, so the shares sum to the whole-body heating.
    assert math.isclose(inner + outer, 1.0, rel_tol=1.0e-12)


def test_tidal_timescale_equal_maxwell_times_give_volume_fractions():
    """With one Maxwell time throughout, the method reduces to volume_fraction."""
    world = _cpl_timescale_world([
        (0.5, _MU, 3.0 * _ETA_AT_PEAK, "tidal_timescale"),
        (1.0, _MU, 3.0 * _ETA_AT_PEAK, "tidal_timescale")])
    _solve(world)
    total = world.get_tidal_heating()
    assert math.isclose(world.get_layer_tidal_heating(0), 0.125 * total, rel_tol=1.0e-12)
    assert math.isclose(world.get_layer_tidal_heating(1), 0.875 * total, rel_tol=1.0e-12)


def test_tidal_timescale_group_shares_its_volume_fraction_beside_other_methods():
    """Beside a volume_fraction core, the method's layers together get their own volume fraction."""
    world = _cpl_timescale_world([
        (0.5, _MU, _ETA_AT_PEAK, "volume_fraction"),
        (0.8, _MU, _ETA_AT_PEAK, "tidal_timescale"),
        (1.0, _MU, 100.0 * _ETA_AT_PEAK, "tidal_timescale")])
    _solve(world)
    total = world.get_tidal_heating()
    assert math.isclose(world.get_layer_tidal_heating(0), 0.125 * total, rel_tol=1.0e-12)
    group_share = (world.get_layer_tidal_heating(1) + world.get_layer_tidal_heating(2)) / total
    assert math.isclose(group_share, 1.0 - 0.125, rel_tol=1.0e-12)
    assert world.get_layer_tidal_heating(1) > world.get_layer_tidal_heating(2)


def test_tidal_timescale_zero_without_moduli():
    """A layer with no shear modulus/viscosity has no Maxwell time, so its weight is 0."""
    world = _cpl_timescale_world([(1.0, 0.0, 0.0, "tidal_timescale")])
    _solve(world)
    assert world.get_layer_tidal_heating(0) == 0.0
    # Beside a usable layer, it gets nothing and the usable layer takes the group's whole share.
    world = _cpl_timescale_world([
        (0.5, 0.0, 0.0, "tidal_timescale"),
        (1.0, _MU, 10.0 * _ETA_AT_PEAK, "tidal_timescale")])
    _solve(world)
    assert world.get_layer_tidal_heating(0) == 0.0
    assert math.isclose(world.get_layer_tidal_heating(1), world.get_tidal_heating(), rel_tol=1.0e-12)


def test_tidal_timescale_uses_the_solved_viscosity_model():
    """A viscosity set by a model (no static constant) reaches the Maxwell time once the EOS is solved."""
    from TidalPy.structures_x.configs import build_world

    def solidliquid_layer(index, radius_fraction, maxwell_periods):
        return {
            "class": "solidliquid", "layer_index": index, "radius_fraction": radius_fraction, "is_tidal": True,
            "tidal_scale_method": "tidal_timescale",
            "material": {"model": "constant", "reference_density_kg_m3": 4000.0,
                         "shear_modulus_static_pa": _MU,
                         "shear_viscosity": {"model": "constant",
                                             "reference_viscosity_pas": maxwell_periods * _ETA_AT_PEAK},
                         "partial_melt": {"model": "off"}}}

    world = build_world({
        "schema_version": "0.2.0", "name": "ts-model", "type": "terrestrial", "radius_m": _R,
        "mass_kg": (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0,
        "tides": {"global_tidal_model": "cpl", "fixed_k": [0.3], "fixed_q": [50.0],
                  "eccentricity_trunc_lvl": 2, "obliquity_trunc_lvl": 0},
        # tau = the forcing period in the inner layer, ten forcing periods (one decade off the peak) in the outer
        "layers": {"inner": solidliquid_layer(0, 0.5, 1.0), "outer": solidliquid_layer(1, 1.0, 10.0)}})
    world.solve_eos()
    _solve(world)
    total = world.get_tidal_heating()
    inner, outer = world.get_layer_tidal_heating(0) / total, world.get_layer_tidal_heating(1) / total
    assert math.isclose((outer / 0.875) / (inner / 0.125), math.exp(-0.5), rel_tol=1.0e-6)
