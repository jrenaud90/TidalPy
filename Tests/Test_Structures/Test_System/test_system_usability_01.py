"""System conveniences: synchronous rotation, the stellar elements on add_world, named evolution rows, copies and
pickles, path-like binary files, and the one-line repr."""
import copy
import math
import pathlib
import pickle

import pytest

from TidalPy.Structures import System, TerrestrialWorld, build_system, build_world

IO_SEMI_MAJOR_AXIS = 4.217e8  # [m]
IO_ECCENTRICITY = 0.0041


def _jupiter_io(**io_options):
    """Jupiter hosting Io, with Io's add_world options."""
    system = System("jovian")
    jupiter = build_world("jupiter_simple")
    jupiter.name = "jupiter"
    io = build_world("io")
    io.name = "io"
    system.add_world(jupiter)
    system.add_world(io, **io_options)
    return system


# =====================================================================================================================
# Synchronous rotation
# =====================================================================================================================
def test_set_synchronous_rotation_sets_the_mean_motion():
    system = _jupiter_io(tidal_host="jupiter", semi_major_axis=IO_SEMI_MAJOR_AXIS, eccentricity=IO_ECCENTRICITY)
    mean_motion = system.calc_orbital_frequency("io")
    assert system["io"].spin_frequency != mean_motion
    assert system.set_synchronous_rotation("io") == mean_motion
    assert system["io"].spin_frequency == mean_motion
    # A world object names it as well as its name does.
    assert system.set_synchronous_rotation(system["io"]) == mean_motion


def test_set_synchronous_rotation_needs_an_orbit():
    system = _jupiter_io()
    with pytest.raises(ValueError, match="no mean motion"):
        system.set_synchronous_rotation("io")


def test_add_world_synchronous():
    system = _jupiter_io(
        tidal_host="jupiter", semi_major_axis=IO_SEMI_MAJOR_AXIS, eccentricity=IO_ECCENTRICITY, synchronous=True)
    assert system["io"].spin_frequency == system.calc_orbital_frequency("io")


@pytest.mark.parametrize("options", [
    pytest.param({"synchronous": True}, id="no-host"),
    pytest.param({"synchronous": True, "tidal_host": "jupiter"}, id="no-semi-major-axis"),
    pytest.param({"stellar_eccentricity": 1.5}, id="unbound-stellar-orbit"),
    pytest.param({"stellar_semi_major_axis": -1.0}, id="negative-stellar-axis"),
])
def test_a_refused_world_is_not_added(options):
    system = System("jovian")
    jupiter = build_world("jupiter_simple")
    jupiter.name = "jupiter"
    system.add_world(jupiter)
    with pytest.raises(ValueError):
        system.add_world(build_world("io"), **options)
    assert len(system) == 1


def test_add_world_stellar_elements():
    system = System("sol")
    sun = build_world("sol")
    sun.name = "sun"
    system.add_world(sun, is_star=True)
    jupiter = build_world("jupiter_simple")
    jupiter.name = "jupiter"
    system.add_world(jupiter, tidal_host="sun", semi_major_axis=7.785472e11, eccentricity=0.0489)
    io = build_world("io")
    io.name = "io"
    system.add_world(
        io,
        tidal_host="jupiter",
        semi_major_axis=IO_SEMI_MAJOR_AXIS,
        stellar_semi_major_axis=7.785472e11,
        stellar_eccentricity=0.0489)
    assert system.get_stellar_semi_major_axis("io") == 7.785472e11
    assert system.get_stellar_eccentricity("io") == 0.0489
    assert math.isfinite(system.calc_insolation_flux("io"))


def test_a_star_hosted_world_takes_no_stellar_elements():
    system = System("sol")
    sun = build_world("sol")
    sun.name = "sun"
    system.add_world(sun, is_star=True)
    with pytest.raises(ValueError, match="its orbit about the star is its tidal orbit"):
        system.add_world(
            build_world("earth_simple"), tidal_host="sun", semi_major_axis=1.5e11, stellar_semi_major_axis=1.5e11)
    assert len(system) == 1


# =====================================================================================================================
# Evolution rows
# =====================================================================================================================
def test_evolution_rows_carry_world_names():
    system = build_system("sol_system")
    system["earth"].solve_eos()
    rows = system.calc_system_evolution()
    assert [row["world_name"] for row in rows] == [world.name for world in system]
    assert system.calc_world_evolution("earth")["world_name"] == "earth"
    pair = system.calc_pair_evolution("earth")
    assert (pair["world_name"], pair["host_name"]) == ("earth", "sun")
    assert (pair["world"]["world_name"], pair["host"]["world_name"]) == ("earth", "sun")
    # The star has no tidal host to name.
    assert system.calc_pair_evolution("sun")["host_name"] is None


# =====================================================================================================================
# Copies, pickles, and binary files
# =====================================================================================================================
def test_copy_is_an_independent_system():
    system = build_system("sol_system")
    for duplicate in (system.copy(), copy.copy(system), copy.deepcopy(system), pickle.loads(pickle.dumps(system))):
        assert isinstance(duplicate, System)
        assert duplicate.get_config_dict() == system.get_config_dict()
        assert [type(world) for world in duplicate] == [type(world) for world in system]
        assert duplicate.source_config == system.source_config
        assert duplicate.source_config is not system.source_config
        assert duplicate["earth"].source_config == system["earth"].source_config
        assert duplicate.get_save_config() == system.get_save_config()
    duplicate = system.copy()
    duplicate.set_eccentricity("earth", 0.2)
    duplicate["earth"].set_spin_frequency(1.0e-5)
    assert system.get_eccentricity("earth") == 0.0167
    assert system["earth"].spin_frequency != 1.0e-5


def test_system_binary_paths_may_be_pathlike(tmp_path):
    system = build_system("sol_system")
    path = pathlib.Path(tmp_path) / "sol.tpyb"
    system.save_binary(path)
    loaded = System()
    loaded.load_binary(path)
    assert loaded.get_config_dict() == system.get_config_dict()
    with pytest.raises(TypeError):
        loaded.load_binary(42)
    with pytest.raises(FileNotFoundError):
        loaded.load_binary(pathlib.Path(tmp_path) / "missing.tpyb")


def test_system_repr():
    system = build_system("sol_system")
    assert repr(system) == "System('Sol System', worlds=['sun', 'earth', 'jupiter'], star='sun')"
    assert repr(System("empty")) == "System('empty', worlds=[], star=None)"
    assert "0x" not in repr(system)


def test_a_world_in_a_copied_system_is_its_own_class():
    system = _jupiter_io(tidal_host="jupiter", semi_major_axis=IO_SEMI_MAJOR_AXIS)
    duplicate = system.copy()
    assert isinstance(duplicate["io"], TerrestrialWorld)
    assert duplicate.get_tidal_host("io").name == "jupiter"
