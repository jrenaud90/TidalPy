"""Copying and pickling worlds through their binary record, binary paths given as os.PathLike, and the class names in
a wrong-class load error."""
import copy
import pathlib
import pickle

import pytest

from TidalPy.Structures import build_world
from TidalPy.Structures.worlds import BaseWorld, GasGiantWorld, StarWorld, TerrestrialWorld, world_from_bytes


@pytest.mark.parametrize("name", ["io", "earth_simple", "jupiter_simple", "sol"])
def test_copy_is_an_independent_world_of_the_same_class(name):
    world = build_world(name)
    copied = world.copy()
    assert type(copied) is type(world)
    assert copied is not world
    assert copied.get_config_dict() == world.get_config_dict()
    assert copied.source_config == world.source_config
    assert copied.source_config is not world.source_config
    assert copied.get_save_config() == world.get_save_config()
    if world.num_layers:
        copied[0].temperature = world[0].temperature + 100.0
        assert copied[0].temperature != world[0].temperature


def test_copy_module_and_pickle_use_copy(tmp_path):
    world = build_world("io")
    world.set_prescribed_heating("mantle", power=1.0e12)
    for duplicate in (copy.copy(world), copy.deepcopy(world), pickle.loads(pickle.dumps(world))):
        assert type(duplicate) is TerrestrialWorld
        assert duplicate.get_config_dict() == world.get_config_dict()
        assert duplicate.prescribed_heating == {"mantle": {"power": 1.0e12}}


def test_a_copy_solves_to_the_same_numbers():
    world = build_world("io")
    copied = world.copy()
    world.solve_eos()
    copied.solve_eos()
    assert not copied.tides_solved
    assert copied.planet_moi_eos == world.planet_moi_eos
    orbit = {"orbital_frequency": 4.11e-5, "spin_frequency": 4.11e-5, "eccentricity": 0.0041, "obliquity": 0.0,
             "semi_major_axis": 4.217e8, "host_mass": 1.898e27}
    assert copied.calc_tides(**orbit) == world.calc_tides(**orbit)


def test_a_star_copy_keeps_its_luminosity():
    star = build_world("sol")
    copied = star.copy()
    assert isinstance(copied, StarWorld)
    assert copied.luminosity == star.luminosity
    assert copied.effective_temperature == star.effective_temperature


def test_world_from_bytes_refuses_another_class():
    world = build_world("io")
    record = pickle.dumps(world)
    assert isinstance(pickle.loads(record), TerrestrialWorld)
    world_record = world.__reduce__()[1][1]
    with pytest.raises(TypeError, match="is a TerrestrialWorld, not a GasGiantWorld"):
        world_from_bytes(GasGiantWorld, world_record)
    with pytest.raises(TypeError, match="BaseWorld or a subclass"):
        world_from_bytes(dict, world_record)
    with pytest.raises(IOError):
        world_from_bytes(TerrestrialWorld, world_record[:-3])


def test_an_empty_wrapper_names_its_class():
    empty = BaseWorld.__new__(BaseWorld)
    with pytest.raises(RuntimeError, match="this BaseWorld holds no world"):
        empty.copy()


def test_binary_paths_may_be_pathlike(tmp_path):
    world = build_world("io")
    path = pathlib.Path(tmp_path) / "io.tpyb"
    world.save_binary(path)
    loaded = build_world("io")
    loaded.load_binary(path)
    assert loaded.get_config_dict() == world.get_config_dict()
    layer = loaded[0]
    standalone_path = pathlib.Path(tmp_path) / "material.tpyb"
    layer.material.save_binary(standalone_path)
    with pytest.raises(TypeError, match="binary file path"):
        world.save_binary(3)


def test_a_wrong_class_load_names_both_classes(tmp_path):
    path = pathlib.Path(tmp_path) / "jupiter.tpyb"
    build_world("jupiter_simple").save_binary(path)
    with pytest.raises(IOError, match="it is a GasGiantWorld file, not a TerrestrialWorld one"):
        build_world("io").load_binary(path)
