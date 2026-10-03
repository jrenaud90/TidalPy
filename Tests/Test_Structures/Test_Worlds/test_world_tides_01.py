"""BaseWorld global (1D) tides with analytic models: CPL heating, layer sharing, defaults, and error paths."""
import math

import pytest

from TidalPy.constants import G
from TidalPy.Structures import build_world
from TidalPy.Structures.worlds.base import BaseWorld


# Io-Jupiter-like orbital state.
_HOST_MASS = 1.898e27
_SMA       = 4.2e8
_N         = 2.05e-5   # mean motion [rad s-1]
_ECC       = 0.0041
_ORBIT = {"orbital_frequency": _N, "spin_frequency": _N, "eccentricity": _ECC, "obliquity": 0.0,
          "semi_major_axis": _SMA, "host_mass": _HOST_MASS}


def _terrestrial_config(tides):
    return {
        "name": "test_terr", "type": "terrestrial", "radius_m": 1.6e6, "mass_kg": 8.9e22,
        "tides": tides,
        "layers": {
            "mantle": {"radius_fraction": 1.0, "tidal_scale": 0.8, "material": {"solid": {
                "eos": {"model": "constant", "reference_density_kg_m3": 3500.0, "bulk_modulus_pa": 2.0e11},
                "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10}}}}
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
    world.calc_tides(**_ORBIT)


def _analytic_cpl(k2, q2, radius, mass_host=_HOST_MASS):
    return (21.0 / 2.0) * (k2 / q2) * G * mass_host**2 * radius**5 * _N * _ECC**2 / _SMA**6


def _config_tides():
    """The config's [tides] table (the packaged one is merged under any user file, so it is always there)."""
    import TidalPy
    return TidalPy.config["tides"]


def test_cpl_world_matches_analytic_heating():
    """A CPL world's heating matches the analytic low-eccentricity CPL rate."""
    world = _cpl_world(k2=0.3, q2=50.0)
    assert world.tide_model_set
    _solve(world)
    assert world.tides_solved
    assert world.get_num_tidal_modes() > 0
    expected = _analytic_cpl(0.3, 50.0, 1.6e6)
    assert math.isclose(world.get_tidal_heating(), expected, rel_tol=5.0e-3)


def test_analytic_layer_heating_sums_to_the_world_heating():
    """With an analytic model the tidal layers share all of the world heating."""
    world = _cpl_world(tidal_scale=0.8)
    _solve(world)
    heat = world.get_tidal_heating()
    layer_heating = [world.get_layer_tidal_heating(i) for i in range(world.num_layers)]
    assert math.isclose(sum(h for h in layer_heating if math.isfinite(h)), heat, rel_tol=1.0e-12)


@pytest.mark.parametrize("state, message", (
    ({"eccentricity": 1.2}, "eccentricity"),
    ({"eccentricity": -0.1}, "eccentricity"),
    ({"semi_major_axis": -4.2e8}, "semi-major axis"),
))
def test_unphysical_orbits_are_rejected(state, message):
    """Unphysical eccentricities and semi-major axes raise a ValueError."""
    world = _cpl_world()
    call = dict(_ORBIT)
    call.update(state)
    with pytest.raises(ValueError, match=message):
        world.calc_tides(**call)


@pytest.mark.parametrize("degrees", ((3, 2), (1, 2), (2, 11)))
def test_degree_range_is_checked_when_set(degrees):
    """An invalid degree range raises when set."""
    world = _cpl_world()
    with pytest.raises(ValueError, match="degree range"):
        world.set_tide_config(min_degree_l=degrees[0], max_degree_l=degrees[1])


def test_potential_derivatives_present():
    """calc_tides fills a nonzero dU/dM."""
    world = _cpl_world()
    _solve(world)
    dUdM, _, _ = world.get_tidal_potential_derivatives().values()
    assert abs(dUdM) > 0.0


def test_unsolved_world_returns_nan():
    """An unsolved world reports NaN heating and no modes."""
    world = _cpl_world()
    assert not world.tides_solved
    assert math.isnan(world.get_tidal_heating())
    assert math.isnan(world.get_layer_tidal_heating(0))
    assert world.get_num_tidal_modes() == 0


def test_heating_scales_with_k_over_q():
    """Doubling k2 doubles the CPL heating."""
    base = _cpl_world(k2=0.3, q2=50.0)
    _solve(base)
    doubled = _cpl_world(k2=0.6, q2=50.0)
    _solve(doubled)
    assert math.isclose(doubled.get_tidal_heating(), 2.0 * base.get_tidal_heating(), rel_tol=1.0e-9)


def test_ctl_world_positive_heating():
    """A CTL world produces positive heating."""
    cfg = _terrestrial_config({
        "global_tidal_model": "ctl", "max_degree_l": 2,
        "eccentricity_trunc_lvl": 2, "obliquity_trunc_lvl": "off",
        "fixed_k": [0.3], "fixed_dt_s": [100.0],
    })
    world = build_world(cfg)
    _solve(world)
    assert world.get_tidal_heating() > 0.0


def test_terrestrial_default_model_is_rheology_requires_eos():
    """Terrestrial worlds default to the rheology model, whose calc_tides needs a solved EOS."""
    world = build_world(_terrestrial_config({}))
    assert world.tide_model_set
    with pytest.raises(RuntimeError):
        _solve(world)


def test_calc_tides_without_model_raises():
    """calc_tides without a tide model raises."""
    world = BaseWorld(world_type="terrestrial", name="bare", radius=1.6e6, mass=8.9e22)
    assert not world.tide_model_set
    with pytest.raises(RuntimeError):
        _solve(world)


def test_config_carries_tide_defaults():
    """The configuration holds the tide defaults."""
    tides = _config_tides()
    assert tides["default_model"]["terrestrial"] == "rheology"
    assert tides["default_model"]["gasgiant"] == "fixed_dt"
    assert len(tides["fixed_k"]) >= 2 and tides["fixed_k"][0] > 0.0
    assert len(tides["fixed_q"]) >= 2 and tides["fixed_q"][0] > 0.0


def test_builder_uses_config_per_degree_defaults():
    """A CPL world without fixed_k and fixed_q takes them from the config defaults."""
    tides = _config_tides()
    k2_default = tides["fixed_k"][0]
    q2_default = tides["fixed_q"][0]
    cfg = _terrestrial_config({
        "global_tidal_model": "cpl", "max_degree_l": 2,
        "eccentricity_trunc_lvl": 2, "obliquity_trunc_lvl": "off",
    })
    world = build_world(cfg)
    _solve(world)
    expected = _analytic_cpl(k2_default, q2_default, 1.6e6)
    assert math.isclose(world.get_tidal_heating(), expected, rel_tol=5.0e-3)
