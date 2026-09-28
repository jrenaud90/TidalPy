"""Stars take the ``[tides.star]`` parameters (Q' near 1e6), not the planet lists of the ``[tides]`` block."""
import copy

import pytest

import TidalPy
from TidalPy.configurations import find_unknown_config_keys, get_packaged_config
from TidalPy.Structures import build_world


def _config(world_type, tides=None):
    config = {"name": "tidal", "type": world_type, "radius_m": 7.0e8, "mass_kg": 2.0e30}
    if world_type == "star":
        config["luminosity_w"] = 3.8e26
    else:
        config["radius_m"] = 7.0e7
        config["mass_kg"] = 1.9e27
        config["layers"] = {"envelope": {"class": "gas", "type": "gas", "radius_fraction": 1.0}}
    if tides is not None:
        config["tides"] = copy.deepcopy(tides)
    return config


def _tides_of(world):
    return world.get_config_dict()["tides"]


def _modified_q(tides):
    """Q' = 3 Q / (2 k2) at degree 2."""
    return 3.0 * tides["fixed_q"][0] / (2.0 * tides["fixed_k"][0])


def test_a_star_takes_the_stellar_table():
    star_table = TidalPy.config["tides"]["star"]
    tides = _tides_of(build_world(_config("star")))
    assert tides["global_tidal_model"] == "fixed_q"
    assert list(tides["fixed_k"]) == pytest.approx(star_table["fixed_k"])
    assert list(tides["fixed_q"]) == pytest.approx(star_table["fixed_q"])
    # An n = 3 polytrope's k2.
    assert tides["fixed_k"][0] == pytest.approx(0.0289, rel=1e-3)
    assert _modified_q(tides) == pytest.approx(1.0e6, rel=1e-2)


def test_a_planet_keeps_the_planet_lists():
    defaults = TidalPy.config["tides"]
    tides = _tides_of(build_world(_config("gasgiant", {"global_tidal_model": "fixed_q"})))
    assert list(tides["fixed_k"]) == pytest.approx(defaults["fixed_k"])
    assert list(tides["fixed_q"]) == pytest.approx(defaults["fixed_q"])


def test_a_star_file_still_wins():
    tides = _tides_of(build_world(_config("star", {"fixed_k": [0.1], "fixed_q": [5.0e5], "max_degree_l": 2})))
    assert tides["fixed_k"][0] == 0.1
    assert tides["fixed_q"][0] == 5.0e5


@pytest.mark.parametrize("name, k2", [("sol", 0.0289), ("trappist1", 0.287)])
def test_bundled_stars_state_their_tides(name, k2):
    tides = _tides_of(build_world(name))
    assert tides["fixed_k"][0] == pytest.approx(k2, rel=1e-3)
    assert _modified_q(tides) == pytest.approx(1.0e6, rel=1e-2)


@pytest.mark.parametrize("config, unknown", [
    pytest.param({"tides": {"star": {"fixed_k": [0.03], "fixed_q": [2.0e4]}}}, [], id="star-table"),
    pytest.param({"tides": {"gasgiant": {"fixed_dt_s": [10.0]}}}, [], id="gasgiant-table"),
    pytest.param({"tides": {"star": {"fixed_kk": [0.03]}}}, ["tides.star.fixed_kk"], id="misspelled-key"),
])
def test_the_star_table_is_a_known_config_key(config, unknown):
    assert find_unknown_config_keys(config, get_packaged_config()) == unknown
