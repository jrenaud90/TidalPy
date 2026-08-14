"""System builder (``TidalPy.structures_x.configs.build_system``).

Turns a system description (a bundled name, a TOML path, or a dict) into a wired ``System``: each
``[worlds.<name>]`` table names a world (bundled name, path, or inline config) plus its host/star roles
and its orbital elements about the tidal host and about the star. These tests cover building from the
bundled example and from dicts, the table-key-as-world-name behavior, template reuse, the star-is-not-
the-host case, the round-trip insolation physics, and the schema validation errors.
"""
import math

import pytest

from TidalPy.structures_x.system.system import System
from TidalPy.structures_x.configs import (
    build_system,
    construct_system,
    validate_system_config,
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
    assert system.host.name == "sun"
    assert system.star.name == "sun"
    assert math.isclose(system.get_semi_major_axis("earth"), AU, rel_tol=1e-9)
    assert math.isclose(system.get_eccentricity("earth"), 0.0167, rel_tol=1e-9)
    assert math.isclose(system.get_stellar_semi_major_axis("jupiter"), 7.785472e11, rel_tol=1e-9)


def test_bundled_sol_system_insolation():
    """The Sun's Stefan-Boltzmann luminosity gives ~the solar constant at Earth."""
    system = build_system("sol_system")
    assert system.get_star_luminosity() > 3.0e26   # ~3.83e26 W (solar)
    flux = system.calc_insolation_flux("earth")
    assert math.isclose(flux, SOLAR_CONSTANT, rel_tol=5e-3)
    # Earth's grey-body equilibrium temperature is ~255 K.
    assert 250.0 < system.calc_equilibrium_temperature("earth") < 260.0


# =====================================================================================================================
# Building from a dict
# =====================================================================================================================
def _earth_moon_sun_config():
    """Earth-Moon-Sun: the tidal host is the Earth but the star is the Sun (distinct orbits)."""
    return {
        "name": "earth_moon_sun",
        "worlds": {
            "sun":   {"world": "sol", "is_star": True},
            "earth": {"world": "earth_simple", "is_host": True,
                      "stellar_semi_major_axis_m": AU, "stellar_eccentricity": 0.0167},
            "moon":  {"world": "earth_simple",   # a stand-in body for the demo
                      "semi_major_axis_m": 3.844e8, "eccentricity": 0.0549,
                      "stellar_semi_major_axis_m": AU, "stellar_eccentricity": 0.0167},
        },
    }


def test_construct_from_dict_host_not_star():
    system = construct_system(_earth_moon_sun_config())
    assert system.num_worlds == 3
    assert system.host.name == "earth"
    assert system.star.name == "sun"
    assert system.host_index == 1
    assert system.star_index == 0
    # The moon's tidal orbit (about Earth) and stellar orbit (about the Sun) differ.
    assert math.isclose(system.get_semi_major_axis("moon"), 3.844e8, rel_tol=1e-9)
    assert math.isclose(system.get_stellar_semi_major_axis("moon"), AU, rel_tol=1e-9)
    # Insolation on the moon comes from the Sun via the moon's stellar orbit (~solar constant).
    assert math.isclose(system.calc_insolation_flux("moon"), SOLAR_CONSTANT, rel_tol=1e-2)


def test_build_from_toml_path(tmp_path):
    toml_text = (
        'schema_version = "0.2.0"\n'
        'name = "two_body"\n\n'
        '[worlds.star]\n'
        'world = "sol"\n'
        'is_host = true\n'
        'is_star = true\n\n'
        '[worlds.planet]\n'
        'world = "earth_simple"\n'
        'semi_major_axis_m = 1.2e11\n'
        'eccentricity = 0.02\n')
    path = tmp_path / "two_body.toml"
    path.write_text(toml_text)
    system = build_system(str(path))
    assert system.name == "two_body"
    assert [w.name for w in system] == ["star", "planet"]
    assert math.isclose(system.get_semi_major_axis("planet"), 1.2e11, rel_tol=1e-9)


def test_template_reuse_under_different_names():
    """The same bundled world template can be added under different system keys."""
    config = {
        "name": "twins",
        "worlds": {
            "sun":     {"world": "sol", "is_host": True, "is_star": True},
            "planet_a": {"world": "earth_simple", "semi_major_axis_m": 1.0e11},
            "planet_b": {"world": "earth_simple", "semi_major_axis_m": 2.0e11},
        },
    }
    system = construct_system(config)
    assert [w.name for w in system] == ["sun", "planet_a", "planet_b"]
    assert system["planet_a"] is not system["planet_b"]
    assert math.isclose(system.get_semi_major_axis("planet_b"), 2.0e11, rel_tol=1e-9)


# =====================================================================================================================
# System.build + config round-trip
# =====================================================================================================================
def test_system_build_staticmethod():
    """System.build mirrors BaseWorld.build and retains the source config."""
    system = System.build("sol_system")
    assert system.num_worlds == 3
    assert [w.name for w in system] == ["sun", "earth", "jupiter"]
    assert system.source_config is not None
    assert system.config is system.source_config


def test_save_to_toml_roundtrip(tmp_path):
    """build_system -> save_to_toml -> build_system reproduces the system (references preserved)."""
    system = build_system("sol_system")
    path = tmp_path / "roundtrip.toml"
    system.save_to_toml(str(path))
    rebuilt = build_system(str(path))
    assert [w.name for w in rebuilt] == [w.name for w in system]
    assert rebuilt.host.name == "sun" and rebuilt.star.name == "sun"
    assert math.isclose(rebuilt.get_semi_major_axis("earth"), system.get_semi_major_axis("earth"))
    assert math.isclose(rebuilt.get_stellar_eccentricity("jupiter"),
                        system.get_stellar_eccentricity("jupiter"))


def test_get_config_dict_expanded():
    """get_config_dict inlines each world with its roles + orbital elements (live-state expansion)."""
    system = build_system("sol_system")
    config = system.get_config_dict()
    assert set(config["worlds"].keys()) == {"sun", "earth", "jupiter"}
    # Each member is inlined as a full world config under the world's system name.
    assert config["worlds"]["earth"]["world"]["name"] == "earth"
    assert config["worlds"]["earth"]["semi_major_axis_m"] > 0.0
    assert config["worlds"]["sun"]["is_host"] is True
    assert config["worlds"]["sun"]["is_star"] is True


def test_save_expanded_roundtrip(tmp_path):
    """A directly-built system (no source_config) saves the self-contained expansion and rebuilds."""
    system = build_system("sol_system")
    system.source_config = None   # force the get_config_dict (expanded, inlined) save path
    path = tmp_path / "expanded.toml"
    system.save_to_toml(str(path))
    rebuilt = build_system(str(path))
    assert [w.name for w in rebuilt] == ["sun", "earth", "jupiter"]
    assert math.isclose(rebuilt.calc_insolation_flux("earth"),
                        system.calc_insolation_flux("earth"), rel_tol=1e-9)


# =====================================================================================================================
# Validation errors
# =====================================================================================================================
def test_validate_missing_worlds():
    with pytest.raises(ValueError, match="at least one"):
        validate_system_config({"name": "empty"})


def test_validate_unknown_system_key():
    with pytest.raises(ValueError, match="Unexpected system-level key"):
        validate_system_config({"worlds": {"a": {"world": "sol"}}, "bogus": 1})


def test_validate_missing_world_source():
    with pytest.raises(ValueError, match="missing the required 'world'"):
        validate_system_config({"worlds": {"a": {"is_host": True}}})


def test_validate_unknown_world_key():
    with pytest.raises(ValueError, match="Unexpected key"):
        validate_system_config({"worlds": {"a": {"world": "sol", "sma": 1.0}}})


def test_validate_multiple_hosts():
    config = {"worlds": {
        "a": {"world": "sol", "is_host": True},
        "b": {"world": "earth_simple", "is_host": True},
    }}
    with pytest.raises(ValueError, match="host worlds"):
        validate_system_config(config)


def test_validate_multiple_stars():
    config = {"worlds": {
        "a": {"world": "sol", "is_star": True},
        "b": {"world": "earth_simple", "is_star": True},
    }}
    with pytest.raises(ValueError, match="star worlds"):
        validate_system_config(config)
