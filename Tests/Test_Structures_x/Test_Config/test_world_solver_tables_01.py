"""A world file's ``[eos_solver]`` and ``[radial_solver]`` tables pin that world's solver settings.

A pinned key wins over the TidalPy configuration for every solve the world runs, a call's own argument wins over the
pinned key, an unpinned key keeps following the configuration, and the tables survive the config round trip.
"""
import copy

import pytest

from TidalPy.constants import G
from TidalPy.structures_x import build_world, build_world_from_dict
from TidalPy.structures_x.configs import EOS_SOLVER_KEYS, RADIAL_SOLVER_KEYS, validate_solver_table

_FREQUENCY = 2.0e-5


def _config(**tables):
    config = {
        "name": "pinned", "type": "terrestrial", "radius_m": 1.6e6, "mass_kg": 8.9e22,
        "layers": {
            "core":   {"class": "physics", "type": "iron", "radius_fraction": 0.5},
            "mantle": {"class": "physics", "type": "mantle_rock", "radius_fraction": 1.0},
        },
    }
    config.update(copy.deepcopy(tables))
    return config


def test_the_key_sets_are_the_configuration_sections():
    import TidalPy
    assert EOS_SOLVER_KEYS == frozenset(TidalPy.config_x["eos_solver"])
    assert RADIAL_SOLVER_KEYS == frozenset(TidalPy.config_x["radial_solver"])


def test_a_world_without_tables_pins_nothing():
    world = build_world(_config())
    assert world.get_solver_defaults() == {}
    assert "eos_solver" not in world.get_config_dict() and "radial_solver" not in world.get_config_dict()


def test_a_pinned_radial_key_reaches_the_love_solve_and_the_tide_paths():
    """A step cap far too small for the shooting method makes every radial solve fail, unless the call lifts it."""
    world = build_world(_config(radial_solver={"max_num_steps": 4}))
    world.solve_eos(G_to_use=G)
    result = world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2, warnings=False)
    assert result["success"] is False
    # The call's own argument wins over the pinned key.
    lifted = world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2, max_num_steps=500000)
    assert lifted["success"] is True, lifted["message"]
    # calc_tides builds its Love solves from the same pinned settings.
    with pytest.raises(RuntimeError):
        world.calc_tides(orbital_frequency=_FREQUENCY, spin_frequency=_FREQUENCY, eccentricity=0.05,
                         obliquity=0.0, semi_major_axis=4.2e8, host_mass=1.9e27)
    world.set_solver_defaults(radial_solver={})
    assert world.get_solver_defaults() == {}
    world.calc_tides(orbital_frequency=_FREQUENCY, spin_frequency=_FREQUENCY, eccentricity=0.05,
                     obliquity=0.0, semi_major_axis=4.2e8, host_mass=1.9e27)


def test_a_pinned_eos_key_reaches_the_eos_solve():
    """One central-pressure iteration cannot converge a two-layer world, unless the call allows more."""
    world = build_world(_config(eos_solver={"max_iters": 1, "pressure_tol": 1.0e-12}))
    assert world.solve_eos(G_to_use=G)["max_iters_hit"] is True
    assert world.solve_eos(G_to_use=G, max_iters=100, pressure_tol=1.0e-8)["max_iters_hit"] is False


def test_the_tables_round_trip_through_the_config_dict():
    tables = {
        "eos_solver": {"integration_method": "RK45", "rtol": 1.0e-8, "slices_per_layer": 40,
                       "solve_temperature": False},
        "radial_solver": {"use_kamata": True, "max_ram_mb": 250, "start_radius_tolerance": 1.0e-4}}
    world = build_world(_config(**tables))
    assert world.get_solver_defaults() == tables
    config = world.get_config_dict()
    assert config["eos_solver"] == tables["eos_solver"] and config["radial_solver"] == tables["radial_solver"]
    rebuilt = build_world_from_dict(config)
    assert rebuilt.get_solver_defaults() == tables
    # A method name comes back in the configuration's spelling whatever case it was given in.
    world.set_solver_defaults(eos_solver={"integration_method": "radau"})
    assert world.get_solver_defaults()["eos_solver"] == {"integration_method": "Radau"}
    # Only the table given changes.
    assert world.get_solver_defaults()["radial_solver"] == tables["radial_solver"]


def test_the_tables_survive_a_save(tmp_path):
    world = build_world(_config(radial_solver={"use_kamata": True}))
    path = str(tmp_path / "pinned.toml")
    world.save_to_toml(path)
    assert build_world(path).get_solver_defaults() == {"radial_solver": {"use_kamata": True}}


@pytest.mark.parametrize("table, section", [
    ({"no_such_key": 1.0}, "eos_solver"),
    ({"rtol": 0.0}, "eos_solver"),
    ({"rtol": "tight"}, "eos_solver"),
    ({"max_iters": 2.5}, "eos_solver"),
    ({"max_iters": True}, "eos_solver"),
    ({"slices_per_layer": 1}, "eos_solver"),
    ({"nondimensionalize": 1}, "eos_solver"),
    ({"use_kamata": "yes"}, "radial_solver"),
    ({"max_num_steps": 0}, "radial_solver"),
    ({"pressure_tol": 1.0e-8}, "radial_solver"),
])
def test_a_bad_table_is_refused(table, section):
    with pytest.raises(ValueError, match=section):
        build_world(_config(**{section: table}))
    with pytest.raises(ValueError, match=section):
        validate_solver_table(section, table, "test")


def test_an_unknown_method_and_a_non_table_are_refused():
    world = build_world(_config())
    with pytest.raises(ValueError, match="integration method"):
        world.set_solver_defaults(eos_solver={"integration_method": "Euler"})
    with pytest.raises(ValueError, match="must be a table"):
        build_world(_config(eos_solver=3))


def test_a_data_file_world_pins_rk45_unless_its_file_says_otherwise():
    from TidalPy.structures_x.configs.world_builder import DATA_FILE_EOS_INTEGRATION_METHOD
    assert DATA_FILE_EOS_INTEGRATION_METHOD == "RK45"
    world = build_world("earth_prem")
    assert world.get_solver_defaults() == {"eos_solver": {"integration_method": "RK45"}}
    config = dict(world.portable_config)
    config["eos_solver"] = {"integration_method": "DOP853", "rtol": 1.0e-9}
    config["radial_solver"] = {"use_kamata": True}
    pinned = build_world(config)
    assert pinned.get_solver_defaults() == {
        "eos_solver": {"integration_method": "DOP853", "rtol": 1.0e-9}, "radial_solver": {"use_kamata": True}}
    # A file that pins other EOS keys still takes RK45 for the method.
    config["eos_solver"] = {"rtol": 1.0e-9}
    assert build_world(config).get_solver_defaults()["eos_solver"] == {"integration_method": "RK45", "rtol": 1.0e-9}


def test_a_star_rejects_the_tables():
    config = {"name": "sun", "type": "star", "radius_m": 6.96e8, "mass_kg": 1.99e30, "luminosity_w": 3.8e26,
              "radial_solver": {"use_kamata": True}}
    with pytest.raises(ValueError, match="star"):
        build_world(config)
