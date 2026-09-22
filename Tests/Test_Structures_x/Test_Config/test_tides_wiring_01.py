"""The ``[tides]`` table of a world configuration, key by key, from the file to the world it builds.

Every key has three sources (the world's own table, the ``[tides]`` block of ``TidalPy_Configs_x.toml``, and a
built-in fallback), and the builder hands some to the tide model and the rest to the world's tide configuration.
These tests follow each key to where it has to arrive.
"""
import copy

import pytest

import TidalPy
from TidalPy.structures_x import build_world


def _config(world_type="terrestrial", tides=None):
    config = {"name": "wired", "type": world_type, "radius_m": 2.0e6, "mass_kg": 1.0e23}
    if world_type == "star":
        config["luminosity_w"] = 3.8e26
    elif world_type == "gasgiant":
        config["layers"] = {"envelope": {"class": "gas", "type": "gas", "radius_fraction": 1.0}}
    else:
        config["layers"] = {"mantle": {"class": "physics", "type": "mantle_rock", "radius_fraction": 1.0}}
    if tides is not None:
        config["tides"] = copy.deepcopy(tides)
    return config


def _tides_of(world):
    return world.get_config_dict()["tides"]


def _assert_per_degree(found, given):
    """A per-degree list is stored for every degree the model can hold; the degrees left out are zero."""
    found = list(found)
    assert found[:len(given)] == given
    assert all(value == 0.0 for value in found[len(given):])


# =====================================================================================================================
# The model and its per-degree parameters
# =====================================================================================================================
@pytest.mark.parametrize("world_type", ["terrestrial", "layered", "gasgiant", "star"])
def test_default_model_follows_the_world_family(world_type):
    expected = TidalPy.config_x["tides"]["default_model"][world_type]
    assert _tides_of(build_world(_config(world_type)))["global_tidal_model"] == expected


@pytest.mark.parametrize("world_type", ["terrestrial", "gasgiant", "star"])
def test_named_model_overrides_the_family_default(world_type):
    tides = {"global_tidal_model": "fixed_q", "fixed_k": [0.25, 0.1], "fixed_q": [80.0, 90.0], "max_degree_l": 3}
    wired = _tides_of(build_world(_config(world_type, tides)))
    assert wired["global_tidal_model"] == "fixed_q"
    _assert_per_degree(wired["fixed_k"], [0.25, 0.1])
    _assert_per_degree(wired["fixed_q"], [80.0, 90.0])


def test_per_degree_time_lags_reach_the_model():
    tides = {"global_tidal_model": "fixed_dt", "fixed_k": [0.4], "fixed_dt": [250.0]}
    wired = _tides_of(build_world(_config("gasgiant", tides)))
    assert wired["global_tidal_model"] == "fixed_dt"
    _assert_per_degree(wired["fixed_k"], [0.4])
    _assert_per_degree(wired["fixed_dt"], [250.0])


# =====================================================================================================================
# The degree, truncation, and Love-method settings
# =====================================================================================================================
def test_every_setting_reaches_the_world():
    tides = {
        "min_degree_l": 2, "max_degree_l": 4, "eccentricity_trunc_lvl": 5, "obliquity_trunc_lvl": 2,
        "tidal_timescale_width_decades": 2.5, "love_method": "homogeneous", "love_fixed_q": 120.0,
        "love_fixed_dt": 45.0}
    world = build_world(_config("terrestrial", tides))
    found = world.get_tide_config()
    for key, value in tides.items():
        assert found[key] == value, key
    # The same settings ride along in the configuration dictionary, next to the model.
    wired = _tides_of(world)
    for key, value in tides.items():
        assert wired[key] == value, key


def test_a_key_left_out_takes_the_package_configuration():
    defaults = TidalPy.config_x["tides"]
    found = build_world(_config("terrestrial", {"max_degree_l": 3})).get_tide_config()
    assert found["max_degree_l"] == 3
    assert found["min_degree_l"] == defaults["min_degree_l"]
    assert found["eccentricity_trunc_lvl"] == defaults["eccentricity_trunc_lvl"]
    assert found["love_method"] == defaults.get("love_method", "radial_solver")
    # Unset lags stay unset: they are written only when the file gives them.
    assert "love_fixed_q" not in found and "love_fixed_dt" not in found


def test_no_tides_table_at_all_still_wires_the_world():
    world = build_world(_config("terrestrial"))
    assert world.tide_model_set
    assert world.get_tide_config()["max_degree_l"] == TidalPy.config_x["tides"]["max_degree_l"]


@pytest.mark.parametrize("spelling, level", [("off", 0), ("gen", 10), ("general", 10), (1, 1), (2, 2)])
def test_obliquity_truncation_spellings(spelling, level):
    found = build_world(_config("terrestrial", {"obliquity_trunc_lvl": spelling})).get_tide_config()
    assert found["obliquity_trunc_lvl"] == level


def test_an_unknown_tides_key_is_rejected():
    with pytest.raises(ValueError, match="fixed_kk"):
        build_world(_config("terrestrial", {"fixed_kk": [0.3]}))


def test_the_rheology_model_is_refused_by_a_star():
    with pytest.raises((ValueError, RuntimeError)):
        world = build_world(_config("star", {"global_tidal_model": "rheology"}))
        world.calc_tides(2.0e-5, 2.0e-5, 0.01, 0.0, 1.0e9, 1.0e30)
