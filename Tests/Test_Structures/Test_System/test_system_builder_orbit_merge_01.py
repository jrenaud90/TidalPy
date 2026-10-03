"""A system file's tidal and stellar orbit elements, and a mutual pair's, merge element by element.

These need the builder to leave an orbit element the file does not give unset (``eccentricity`` absent is not 0).
"""
import pytest

from TidalPy.Structures.configs import build_system

AU = 1.495978707e11  # [m]


def _config(planet):
    return {"name": "t", "worlds": {"star": {"world": "sol", "is_star": True},
                                    "p": dict({"world": "earth_simple", "tidal_host": "star"}, **planet)}}


@pytest.mark.parametrize("planet, eccentricity", [
    ({"eccentricity": 0.1, "stellar_semi_major_axis_m": AU}, 0.1),
    ({"semi_major_axis_m": AU, "stellar_eccentricity": 0.3}, 0.3),
    ({"semi_major_axis_m": AU, "eccentricity": 0.1, "stellar_semi_major_axis_m": AU}, 0.1),
    ({"stellar_semi_major_axis_m": AU}, 0.0),
], ids=["tidal_e_stellar_a", "tidal_a_stellar_e", "both_a_tidal_e", "no_e"])
def test_star_hosted_world_keeps_every_element_its_file_gives(planet, eccentricity):
    system = build_system(_config(planet))
    assert system.get_eccentricity("p") == eccentricity
    assert system.get_stellar_eccentricity("p") == eccentricity
    assert system.get_semi_major_axis("p") == AU
    assert system.get_stellar_semi_major_axis("p") == AU


def test_star_hosted_world_with_two_eccentricities_is_refused():
    with pytest.raises(ValueError, match="different values"):
        build_system(_config({"semi_major_axis_m": AU, "eccentricity": 0.1, "stellar_eccentricity": 0.2}))


def test_mutual_pair_takes_each_element_from_the_member_that_gives_it():
    system = build_system({"name": "em", "worlds": {
        "earth": {"world": "earth_simple", "tidal_host": "moon", "eccentricity": 0.055},
        "moon": {"world": "luna", "tidal_host": "earth", "semi_major_axis_m": 3.844e8}}})
    for name in ("earth", "moon"):
        assert system.get_eccentricity(name) == 0.055
        assert system.get_semi_major_axis(name) == 3.844e8


def _mutual_pair_with_elements_on_both():
    return build_system({"name": "em", "worlds": {
        "earth": {"world": "earth_simple", "tidal_host": "moon", "semi_major_axis_m": 3.844e8, "eccentricity": 0.055},
        "moon": {"world": "luna", "tidal_host": "earth", "semi_major_axis_m": 3.844e8, "eccentricity": 0.055}}})


@pytest.mark.parametrize("member", ["earth", "moon"])
def test_setting_one_member_of_a_mutual_pair_sets_the_shared_orbit(member):
    system = _mutual_pair_with_elements_on_both()
    system.set_semi_major_axis(member, 3.85e8)
    system.set_eccentricity(member, 0.06)
    for name in ("earth", "moon"):
        assert system.get_semi_major_axis(name) == 3.85e8
        assert system.get_eccentricity(name) == 0.06


def test_a_rebuilt_mutual_pair_can_be_updated_from_one_side():
    from TidalPy.Structures.configs import build_system
    rebuilt = build_system(_mutual_pair_with_elements_on_both().get_config_dict())
    rebuilt.set_semi_major_axis("moon", 3.9e8)
    assert rebuilt.get_semi_major_axis("earth") == 3.9e8
