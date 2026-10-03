"""``build_system`` turns a bundled name, TOML path, or dict into a wired ``System``, and validates the schema."""
import math

import pytest

from TidalPy.Structures.system.system import System
from TidalPy.Structures.configs import (
    build_system,
    validate_system_config,
    available_systems,
)

AU = 1.495978707e11
SOLAR_CONSTANT = 1361.0   # W/m^2 at 1 AU


# =====================================================================================================================
# Building from the bundled example
# =====================================================================================================================
def test_build_bundled_sol_system():
    system = build_system("sol_system")
    assert system.name == "Sol System"
    assert system.num_worlds == 3
    # Worlds are identified by their [worlds.<name>] table keys.
    assert [w.name for w in system] == ["sun", "earth", "jupiter"]
    assert system.get_tidal_host("earth").name == "sun"
    assert system.get_tidal_host("jupiter").name == "sun"
    assert system.get_tidal_host("sun") is None
    assert system.star.name == "sun"
    assert math.isclose(system.get_semi_major_axis("earth"), AU, rel_tol=1e-9)
    assert math.isclose(system.get_eccentricity("earth"), 0.0167, rel_tol=1e-9)
    assert math.isclose(system.get_stellar_semi_major_axis("jupiter"), 7.785472e11, rel_tol=1e-9)


def test_bundled_sol_system_insolation():
    """The Sun's luminosity gives about the solar constant and a 255 K equilibrium temperature at Earth."""
    system = build_system("sol_system")
    assert system.get_star_luminosity() > 3.0e26
    assert math.isclose(system.calc_insolation_flux("earth"), SOLAR_CONSTANT, rel_tol=5e-3)
    assert 250.0 < system.calc_equilibrium_temperature("earth") < 260.0


# =====================================================================================================================
# Building from a dict or a path
# =====================================================================================================================
def _earth_moon_sun_config():
    """The Earth and the Moon host each other, and the Sun is the star."""
    return {
        "name": "earth_moon_sun",
        "worlds": {
            "sun":   {"world": "sol", "is_star": True},
            # Declared before the Moon, which it names all the same; the pair shares the orbit the Moon states.
            "earth": {"world": "earth_simple", "tidal_host": "moon",
                      "stellar_semi_major_axis_m": AU, "stellar_eccentricity": 0.0167},
            "moon":  {"world": "earth_simple",   # a stand-in body
                      "tidal_host": "earth",
                      "semi_major_axis_m": 3.844e8, "eccentricity": 0.0549,
                      "stellar_semi_major_axis_m": AU, "stellar_eccentricity": 0.0167},
        },
    }


def test_construct_from_dict_host_not_star():
    system = build_system(_earth_moon_sun_config())
    assert system.num_worlds == 3
    assert system.get_tidal_host("moon").name == "earth"
    assert system.get_tidal_host("earth").name == "moon"
    assert system.get_tidal_host_index("moon") == 1
    assert system.is_mutual_pair("earth") and system.is_mutual_pair("moon")
    assert not system.has_tidal_host("sun")
    assert system.star.name == "sun"
    assert system.star_index == 0
    # The Earth states no orbit of its own, so it takes the one it shares with the Moon.
    assert math.isclose(system.get_semi_major_axis("earth"), 3.844e8, rel_tol=1e-9)
    assert math.isclose(system.get_eccentricity("earth"), 0.0549, rel_tol=1e-9)
    assert math.isclose(system.get_semi_major_axis("moon"), 3.844e8, rel_tol=1e-9)
    assert math.isclose(system.get_stellar_semi_major_axis("moon"), AU, rel_tol=1e-9)
    assert math.isclose(system.calc_insolation_flux("moon"), SOLAR_CONSTANT, rel_tol=1e-2)


def test_build_from_toml_path(tmp_path):
    toml_text = (
        'schema_version = "0.2.0"\n'
        'name = "two_body"\n\n'
        '[worlds.star]\n'
        'world = "sol"\n'
        'is_star = true\n\n'
        '[worlds.planet]\n'
        'world = "earth_simple"\n'
        'tidal_host = "star"\n'
        'semi_major_axis_m = 1.2e11\n'
        'eccentricity = 0.02\n')
    path = tmp_path / "two_body.toml"
    path.write_text(toml_text)
    system = build_system(str(path))
    assert system.name == "two_body"
    assert [w.name for w in system] == ["star", "planet"]
    assert math.isclose(system.get_semi_major_axis("planet"), 1.2e11, rel_tol=1e-9)


def test_template_reuse_under_different_names():
    """One bundled world template can be added under different system keys."""
    config = {
        "name": "twins",
        "worlds": {
            "sun":     {"world": "sol", "is_star": True},
            "planet_a": {"world": "earth_simple", "tidal_host": "sun", "semi_major_axis_m": 1.0e11},
            "planet_b": {"world": "earth_simple", "tidal_host": "sun", "semi_major_axis_m": 2.0e11},
        },
    }
    system = build_system(config)
    assert [w.name for w in system] == ["sun", "planet_a", "planet_b"]
    assert system["planet_a"] is not system["planet_b"]
    assert math.isclose(system.get_semi_major_axis("planet_b"), 2.0e11, rel_tol=1e-9)


# =====================================================================================================================
# System.build and the config round trip
# =====================================================================================================================
def test_system_build_staticmethod():
    """System.build mirrors BaseWorld.build and keeps the source config."""
    system = System.build("sol_system")
    assert system.num_worlds == 3
    assert [w.name for w in system] == ["sun", "earth", "jupiter"]
    assert system.source_config is not None
    assert system.config is system.source_config


def test_save_to_toml_roundtrip(tmp_path):
    system = build_system("sol_system")
    path = tmp_path / "roundtrip.toml"
    system.save_to_toml(str(path))
    rebuilt = build_system(str(path))
    assert [w.name for w in rebuilt] == [w.name for w in system]
    assert rebuilt.get_tidal_host("earth").name == "sun" and rebuilt.star.name == "sun"
    assert math.isclose(rebuilt.get_semi_major_axis("earth"), system.get_semi_major_axis("earth"))
    assert math.isclose(
        rebuilt.get_stellar_eccentricity("jupiter"), system.get_stellar_eccentricity("jupiter"))


def test_get_config_dict_expanded():
    """get_config_dict inlines each world with its roles and orbital elements."""
    config = build_system("sol_system").get_config_dict()
    assert set(config["worlds"].keys()) == {"sun", "earth", "jupiter"}
    assert config["worlds"]["earth"]["world"]["name"] == "earth"
    assert config["worlds"]["earth"]["semi_major_axis_m"] > 0.0
    assert config["worlds"]["earth"]["tidal_host"] == "sun"
    assert "tidal_host" not in config["worlds"]["sun"]
    assert config["worlds"]["sun"]["is_star"] is True


def test_a_star_hosted_world_is_saved_without_stellar_elements(tmp_path):
    """A world whose tidal host is the star writes its one orbit as its tidal elements (add_world refuses stellar
    elements for it), while a moon keeps its separate orbit about the star; both round-trip."""
    system = build_system({"name": "s", "worlds": {
        "sun": {"world": "sol", "is_star": True},
        "planet": {"world": "earth_simple", "tidal_host": "sun", "semi_major_axis_m": AU, "eccentricity": 0.0167},
        "moon": {"world": "luna", "tidal_host": "planet", "semi_major_axis_m": 3.844e8, "eccentricity": 0.0549,
                 "stellar_semi_major_axis_m": AU, "stellar_eccentricity": 0.0167}}})
    for config in (system.get_config_dict(), system.get_save_config()):
        planet, moon = config["worlds"]["planet"], config["worlds"]["moon"]
        assert "stellar_semi_major_axis_m" not in planet and "stellar_eccentricity" not in planet
        assert planet["semi_major_axis_m"] == AU and planet["eccentricity"] == 0.0167
        assert moon["stellar_semi_major_axis_m"] == AU and moon["stellar_eccentricity"] == 0.0167

    path = tmp_path / "star_hosted.toml"
    system.save_to_toml(str(path))
    rebuilt = build_system(str(path))
    for name in ("planet", "moon"):
        assert rebuilt.get_stellar_semi_major_axis(name) == system.get_stellar_semi_major_axis(name)
        assert rebuilt.get_stellar_eccentricity(name) == system.get_stellar_eccentricity(name)
        assert rebuilt.calc_insolation_flux(name) == system.calc_insolation_flux(name)
    assert rebuilt.get_config_dict() == system.get_config_dict()


def test_save_expanded_roundtrip(tmp_path):
    """A system with no source config saves the self-contained expansion and rebuilds."""
    system = build_system("sol_system")
    # Forces the expanded, inlined save path.
    system.source_config = None
    path = tmp_path / "expanded.toml"
    system.save_to_toml(str(path))
    rebuilt = build_system(str(path))
    assert [w.name for w in rebuilt] == ["sun", "earth", "jupiter"]
    assert math.isclose(rebuilt.calc_insolation_flux("earth"), system.calc_insolation_flux("earth"), rel_tol=1e-9)


def test_available_systems_lists_bundled_systems_only():
    systems = available_systems()
    assert "sol_system" in systems
    assert "earth_simple" not in systems
    assert "sol" not in systems


# =====================================================================================================================
# Validation errors
# =====================================================================================================================
def _pair(world_a, world_b):
    return {"worlds": {"a": world_a, "b": world_b}}


@pytest.mark.parametrize("config, match", [
    pytest.param({"name": "empty"}, "at least one", id="no-worlds"),
    pytest.param({"worlds": {"a": {"world": "sol"}}, "bogus": 1}, "Unexpected system-level key", id="system-key"),
    pytest.param({"worlds": {"a": {"is_star": True}}}, "missing the required 'world'", id="no-world-source"),
    pytest.param({"worlds": {"a": {"world": "sol", "sma": 1.0}}}, "Unexpected key", id="world-key"),
    pytest.param(
        _pair({"world": "sol", "is_host": True}, {"world": "earth_simple", "semi_major_axis_m": 1.0e11}),
        "tidal_host",
        id="is-host-replaced"),
    *[
        pytest.param(
            _pair({"world": "sol"}, {"world": "earth_simple", "tidal_host": host, "semi_major_axis_m": 1.0e11}),
            "tidal host|tidal_host",
            id=f"tidal-host-{host}")
        for host in ("nobody", "b", 3)
    ],
    # Orbital elements are about the tidal host, so they mean nothing without one.
    pytest.param(
        _pair({"world": "sol"}, {"world": "earth_simple", "semi_major_axis_m": 1.0e11}),
        "no 'tidal_host'",
        id="orbit-without-host"),
    pytest.param(
        _pair({"world": "sol", "is_star": True}, {"world": "earth_simple", "is_star": True}),
        "star worlds",
        id="two-stars"),
])
def test_validate_rejects_a_bad_system(config, match):
    with pytest.raises(ValueError, match=match):
        validate_system_config(config)


def test_validate_a_stellar_orbit_needs_no_tidal_host():
    validate_system_config(
        _pair({"world": "sol", "is_star": True}, {"world": "earth_simple", "stellar_semi_major_axis_m": 1.0e11}))


def test_building_a_world_as_a_system_names_build_world():
    with pytest.raises(ValueError, match="world configuration.*build_world"):
        build_system("sol")
