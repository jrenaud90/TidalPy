"""The ``[tides]`` table of a world configuration, key by key, from the file to the world it builds.

Every key has three sources (the world's own table, the ``[tides]`` block of ``TidalPy_Configs_x.toml``, and a
built-in fallback), and the builder hands some to the tide model and the rest to the world's tide configuration.
These tests follow each key to where it has to arrive.
"""
import copy

import pytest

import TidalPy
from TidalPy.structures_x import build_world
from TidalPy.Tides_x.eccentricity import promote_eccentricity_truncation


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
    tides = {"global_tidal_model": "fixed_dt", "fixed_k": [0.4], "fixed_dt_s": [250.0]}
    wired = _tides_of(build_world(_config("gasgiant", tides)))
    assert wired["global_tidal_model"] == "fixed_dt"
    _assert_per_degree(wired["fixed_k"], [0.4])
    _assert_per_degree(wired["fixed_dt_s"], [250.0])


# =====================================================================================================================
# The degree, truncation, and Love-method settings
# =====================================================================================================================
# =====================================================================================================================
# A per-degree list that stops short of max_degree_l
# =====================================================================================================================
def _short_list_warnings(record):
    return [entry for entry in record if "stop short of max_degree_l" in str(entry.message)]


def test_a_short_list_the_model_reads_warns():
    import warnings
    tides = {"global_tidal_model": "fixed_q", "fixed_k": [0.3, 0.2], "fixed_q": [50.0], "max_degree_l": 3}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        world = build_world(_config("terrestrial", tides))
    found = _short_list_warnings(record)
    assert len(found) == 1
    message = str(found[0].message)
    assert "fixed_q (1)" in message and "fixed_k" not in message and "needs 2 entries" in message
    # The world still builds, with the missing degree zero-filled.
    _assert_per_degree(_tides_of(world)["fixed_q"], [50.0])


def test_a_short_list_the_model_does_not_read_is_silent():
    import warnings
    # A ctl model never reads fixed_q, and the rheology model reads no list at all.
    for tides in ({"global_tidal_model": "ctl", "fixed_k": [0.3, 0.2], "fixed_dt_s": [100.0, 90.0],
                   "fixed_q": [50.0], "max_degree_l": 3},
                  {"global_tidal_model": "rheology", "fixed_k": [0.3], "max_degree_l": 3}):
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            build_world(_config("terrestrial", tides))
        assert not _short_list_warnings(record)


def test_the_short_list_switch_silences_it(monkeypatch):
    import warnings
    config_x = dict(TidalPy.config_x)
    config_x["warnings"] = {"short_degree_list": False}
    monkeypatch.setattr(TidalPy, "config_x", config_x)
    tides = {"global_tidal_model": "fixed_q", "fixed_k": [0.3], "fixed_q": [50.0], "max_degree_l": 3}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        build_world(_config("terrestrial", tides))
    assert not _short_list_warnings(record)
    assert TidalPy.config_x["warnings"].get("short_degree_list") is False


def test_every_setting_reaches_the_world():
    tides = {
        "min_degree_l": 2, "max_degree_l": 4, "eccentricity_trunc_lvl": 10, "obliquity_trunc_lvl": 2,
        "layer_tidal_heating": False, "love_method": "homogeneous", "love_fixed_q": 120.0,
        "love_fixed_dt_s": 45.0}
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
    # A configuration file written before the levels changed can hold an untabulated level, which is promoted.
    assert found["eccentricity_trunc_lvl"] == promote_eccentricity_truncation(defaults["eccentricity_trunc_lvl"])
    assert found["love_method"] == defaults.get("love_method", "radial_solver")
    # Unset lags stay unset: they are written only when the file gives them.
    assert "love_fixed_q" not in found and "love_fixed_dt_s" not in found


def test_no_tides_table_at_all_still_wires_the_world():
    world = build_world(_config("terrestrial"))
    assert world.tide_model_set
    assert world.get_tide_config()["max_degree_l"] == TidalPy.config_x["tides"]["max_degree_l"]


@pytest.mark.parametrize("spelling, level", [("off", 0), ("gen", "gen"), ("general", "gen"), (2, 2), (4, 4),
                                             (1, 2), (10, "gen")])
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


def test_set_tide_config_changes_only_the_given_settings():
    """An omitted argument keeps its stored value instead of resetting to a default."""
    from TidalPy.structures_x.configs import build_world
    world = build_world("io")
    world.set_tide_config(min_degree_l=2, max_degree_l=3, eccentricity_truncation=10, obliquity_truncation=0,
                          layer_tidal_heating=False, love_method="homogeneous", love_fixed_q=80.0)
    before = world.get_tide_config()
    world.set_tide_config(eccentricity_truncation=10)
    after = world.get_tide_config()
    assert after["eccentricity_trunc_lvl"] == 10
    for key in ("min_degree_l", "max_degree_l", "obliquity_trunc_lvl", "layer_tidal_heating",
                "love_method", "love_fixed_q"):
        assert after[key] == before[key]
    world.set_tide_config(love_fixed_q=float("nan"))   # NaN clears it
    assert "love_fixed_q" not in world.get_tide_config()
    with pytest.raises(NotImplementedError):
        world.set_tide_config(eccentricity_truncation=7)


def test_exact_eccentricity_from_a_tides_table():
    """A [tides] table can ask for the exact eccentricity functions and their tolerance, and the world writes them back."""
    world = build_world(_config("terrestrial", {"eccentricity_trunc_lvl": "exact",
                                                "eccentricity_exact_tolerance": 1.0e-5}))
    tides = world.get_tide_config()
    assert tides["eccentricity_trunc_lvl"] == "exact"
    assert tides["eccentricity_exact_tolerance"] == 1.0e-5
    rebuilt = build_world(world.get_config_dict()).get_tide_config()
    assert rebuilt["eccentricity_trunc_lvl"] == "exact"
    assert rebuilt["eccentricity_exact_tolerance"] == 1.0e-5
    with pytest.raises(ValueError):
        build_world(_config("terrestrial", {"eccentricity_trunc_lvl": "exact", "eccentricity_exact_tolerance": 1.5}))
