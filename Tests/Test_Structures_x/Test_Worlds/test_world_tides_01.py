"""Tests for LayeredWorld global (1D) tidal dissipation (calc_tides).

Covers the analytic tide path end-to-end through the world builder: attaching a tide model
+ config from the ``[tides]`` table, the collapse matching the analytic CPL heating rate,
per-layer heating distribution by ``tidal_scale``, the orbital potential derivatives, the
per-world-family default models, and the error paths (no model / rheology-not-yet-wired).
"""
import math

import pytest

from TidalPy.constants import G
from TidalPy.structures_x import build_world
from TidalPy.structures_x.worlds.layered import LayeredWorld


# Io-Jupiter-like orbital state.
_HOST_MASS = 1.898e27
_SMA       = 4.2e8
_N         = 2.05e-5   # mean motion [rad s-1]
_ECC       = 0.0041


def _terrestrial_config(tides):
    return {
        "name": "test_terr", "type": "terrestrial", "radius_m": 1.6e6, "mass_kg": 8.9e22,
        "tides": tides,
        "layers": {
            "mantle": {
                "class": "physics", "type": "mantle_rock", "radius_fraction": 1.0,
                "tidal_scale": 0.8,
                "shear_modulus_static_pa": 6.0e10, "bulk_modulus_static_pa": 2.0e11,
            }
        },
    }


def _cpl_world(k2=0.3, q2=50.0, tidal_scale=0.8):
    cfg = _terrestrial_config({
        "global_tidal_model": "cpl", "max_degree_l": 2,
        "eccentricity_trunc_lvl": 2, "obliquity_trunc_lvl": "off",
        "fixed_k": [k2], "fixed_q": [q2],
    })
    cfg["layers"]["mantle"]["tidal_scale"] = tidal_scale
    return build_world(cfg)


def _solve(world):
    world.calc_tides(orbital_frequency=_N, spin_frequency=_N, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST_MASS)


def _analytic_cpl(k2, q2, radius, mass_host=_HOST_MASS):
    return (21.0 / 2.0) * (k2 / q2) * G * mass_host**2 * radius**5 * _N * _ECC**2 / _SMA**6


# =====================================================================================================================
# Analytic CPL world end-to-end
# =====================================================================================================================
def test_cpl_world_matches_analytic_heating():
    world = _cpl_world(k2=0.3, q2=50.0)
    assert world.tide_model_set
    _solve(world)
    assert world.tides_solved
    assert world.get_num_tidal_modes() > 0
    expected = _analytic_cpl(0.3, 50.0, 1.6e6)
    assert math.isclose(world.get_tidal_heating(), expected, rel_tol=5.0e-3)


def test_layer_heating_is_world_heating_times_tidal_scale():
    world = _cpl_world(tidal_scale=0.8)
    _solve(world)
    heat = world.get_tidal_heating()
    assert math.isclose(world.get_layer_tidal_heating(0), heat * 0.8, rel_tol=1.0e-9)


def test_potential_derivatives_present():
    world = _cpl_world()
    _solve(world)
    dUdM, dUdw, dUdO = world.get_tidal_potential_derivatives()
    assert abs(dUdM) > 0.0


def test_unsolved_world_returns_nan():
    world = _cpl_world()
    assert not world.tides_solved
    assert math.isnan(world.get_tidal_heating())
    assert math.isnan(world.get_layer_tidal_heating(0))
    assert world.get_num_tidal_modes() == 0


def test_heating_scales_with_k_over_q():
    base = _cpl_world(k2=0.3, q2=50.0); _solve(base)
    doubled = _cpl_world(k2=0.6, q2=50.0); _solve(doubled)
    assert math.isclose(doubled.get_tidal_heating(), 2.0 * base.get_tidal_heating(), rel_tol=1.0e-9)


# =====================================================================================================================
# CTL world
# =====================================================================================================================
def test_ctl_world_positive_heating():
    cfg = _terrestrial_config({
        "global_tidal_model": "ctl", "max_degree_l": 2,
        "eccentricity_trunc_lvl": 2, "obliquity_trunc_lvl": "off",
        "fixed_k": [0.3], "fixed_dt": [100.0],
    })
    world = build_world(cfg)
    _solve(world)
    assert world.get_tidal_heating() > 0.0


# =====================================================================================================================
# Per-world-family default models + error paths
# =====================================================================================================================
def test_terrestrial_default_model_is_rheology_requires_eos():
    """Terrestrial worlds default to the rheology model; calc_tides needs the EOS solved first.

    The rheology model runs the radial solver per tidal frequency, so without a prior
    solve_eos the world raises. The end-to-end rheology path is exercised in
    test_world_tides_rheology_01.py.
    """
    cfg = _terrestrial_config({})  # no global_tidal_model -> default (rheology)
    world = build_world(cfg)
    assert world.tide_model_set
    with pytest.raises(RuntimeError):
        _solve(world)  # EOS not solved yet


def test_calc_tides_without_model_raises():
    world = LayeredWorld(world_type="terrestrial", name="bare", radius_m=1.6e6, mass_kg=8.9e22)
    assert not world.tide_model_set
    with pytest.raises(RuntimeError):
        _solve(world)


# =====================================================================================================================
# Config-driven defaults (defaultc_x.py -> TidalPy_Configs_x.toml -> config_x)
# =====================================================================================================================
def test_config_x_carries_tide_defaults():
    """The _x config is the single source of truth for the tide defaults."""
    import TidalPy
    tides = (TidalPy.config_x or {}).get("tides")
    if not tides:
        pytest.skip("config_x has no [tides] block (stale config dir); regenerate to test.")
    assert tides["default_model"]["terrestrial"] == "rheology"
    assert tides["default_model"]["gasgiant"] == "fixed_dt"
    assert len(tides["fixed_k"]) >= 2 and tides["fixed_k"][0] > 0.0
    assert len(tides["fixed_q"]) >= 2 and tides["fixed_q"][0] > 0.0


def test_builder_uses_config_x_per_degree_defaults():
    """A cpl world that omits fixed_k/fixed_q picks them up from the config_x defaults."""
    import TidalPy
    tides = (TidalPy.config_x or {}).get("tides")
    if not tides:
        pytest.skip("config_x has no [tides] block (stale config dir); regenerate to test.")
    k2_default = tides["fixed_k"][0]
    q2_default = tides["fixed_q"][0]

    # No fixed_k / fixed_q in the world config -> defaults flow from config_x.
    cfg = _terrestrial_config({
        "global_tidal_model": "cpl", "max_degree_l": 2,
        "eccentricity_trunc_lvl": 2, "obliquity_trunc_lvl": "off",
    })
    world = build_world(cfg)
    _solve(world)
    expected = _analytic_cpl(k2_default, q2_default, 1.6e6)
    assert math.isclose(world.get_tidal_heating(), expected, rel_tol=5.0e-3)
