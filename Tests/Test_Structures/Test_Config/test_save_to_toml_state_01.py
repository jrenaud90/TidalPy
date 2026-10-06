"""save_to_toml writes a world or system as it is now, with its file references found from the folder saved into."""
import shutil
from pathlib import Path

import pytest
import toml

from TidalPy.Structures import build_system, build_world
from TidalPy.Structures.configs import worldpack

_NEW_OBLIQUITY = 0.123        # [rad]
_NEW_TEMPERATURE = 1555.0     # [K]
_NEW_SEMI_MAJOR_AXIS = 8.0e8  # [m]


def _bundled_file(name):
    """A bundled WorldPack file's path in the packaged folder."""
    return Path(worldpack.PACKAGED_WORLDPACK_DIR) / name


def test_a_change_after_the_build_is_saved(tmp_path):
    world = build_world("io")
    world.set_obliquity(_NEW_OBLIQUITY)
    world.mantle.temperature = _NEW_TEMPERATURE
    path = tmp_path / "io_changed.toml"
    world.save_to_toml(path)
    rebuilt = build_world(str(path))
    assert rebuilt.obliquity == pytest.approx(_NEW_OBLIQUITY)
    assert rebuilt.mantle.temperature == pytest.approx(_NEW_TEMPERATURE)


def test_an_unchanged_data_file_world_keeps_its_reference(tmp_path):
    """A solve sets the layer masses but changes nothing the user gave, so the compact form is still written."""
    world = build_world("earth_prem")
    world.solve_eos()
    path = tmp_path / "prem.toml"
    world.save_to_toml(path)
    saved = toml.load(path)
    assert "data_file" in saved
    assert "layers" not in saved or all("material" not in layer for layer in saved["layers"].values())
    assert len(build_world(str(path)).get_config_dict()["layers"]) == len(world.get_config_dict()["layers"])


def test_a_changed_data_file_world_is_saved_whole(tmp_path):
    world = build_world("earth_prem")
    world.set_obliquity(_NEW_OBLIQUITY)
    path = tmp_path / "prem_changed.toml"
    world.save_to_toml(path)
    saved = toml.load(path)
    assert "data_file" not in saved
    rebuilt = build_world(str(path))
    assert rebuilt.obliquity == pytest.approx(_NEW_OBLIQUITY)
    assert len(rebuilt.get_config_dict()["layers"]) == len(world.get_config_dict()["layers"])


def test_a_relative_data_file_is_found_from_the_new_folder(tmp_path):
    source_dir = tmp_path / "source"
    destination_dir = tmp_path / "destination"
    source_dir.mkdir()
    destination_dir.mkdir()
    shutil.copy(_bundled_file("PREM.csv"), source_dir / "my_profile.csv")
    config = toml.load(_bundled_file("earth_prem.toml"))
    config["data_file"] = "my_profile.csv"
    with open(source_dir / "my_earth.toml", "w", encoding="utf-8") as world_file:
        toml.dump(config, world_file)

    world = build_world(str(source_dir / "my_earth.toml"))
    world.save_to_toml(destination_dir / "copy.toml")
    saved = toml.load(destination_dir / "copy.toml")
    assert saved["data_file"] == "../source/my_profile.csv"
    assert build_world(str(destination_dir / "copy.toml")).name == world.name
    # Saved beside the file, the reference stays as given.
    world.save_to_toml(source_dir / "copy.toml")
    assert toml.load(source_dir / "copy.toml")["data_file"] == "my_profile.csv"


def test_a_system_saves_its_current_orbits_and_keeps_unchanged_references(tmp_path):
    system = build_system("sol_system")
    system.set_semi_major_axis("jupiter", _NEW_SEMI_MAJOR_AXIS)
    path = tmp_path / "system.toml"
    system.save_to_toml(path)
    saved = toml.load(path)
    assert saved["worlds"]["jupiter"]["semi_major_axis_m"] == pytest.approx(_NEW_SEMI_MAJOR_AXIS)
    assert all(isinstance(entry["world"], str) for entry in saved["worlds"].values())
    rebuilt = build_system(str(path))
    assert rebuilt.get_config_dict()["worlds"]["jupiter"]["semi_major_axis_m"] == pytest.approx(_NEW_SEMI_MAJOR_AXIS)


def test_a_changed_system_member_is_saved_inline(tmp_path):
    system = build_system("sol_system")
    system["earth"].set_obliquity(_NEW_OBLIQUITY)
    path = tmp_path / "system.toml"
    system.save_to_toml(path)
    saved = toml.load(path)
    assert isinstance(saved["worlds"]["earth"]["world"], dict)
    assert isinstance(saved["worlds"]["jupiter"]["world"], str)
    assert build_system(str(path))["earth"].obliquity == pytest.approx(_NEW_OBLIQUITY)


def test_a_relative_member_path_is_found_from_the_new_folder(tmp_path):
    source_dir = tmp_path / "source"
    destination_dir = tmp_path / "destination"
    source_dir.mkdir()
    destination_dir.mkdir()
    shutil.copy(_bundled_file("io.toml"), source_dir / "my_io.toml")
    system_config = {
        "name": "pair",
        "worlds": {
            "jupiter": {"world": "jupiter_simple"},
            "io": {"world": "my_io.toml", "tidal_host": "jupiter", "semi_major_axis_m": 4.217e8,
                   "eccentricity": 0.0041}}}
    with open(source_dir / "pair.toml", "w", encoding="utf-8") as system_file:
        toml.dump(system_config, system_file)

    system = build_system(str(source_dir / "pair.toml"))
    system.save_to_toml(destination_dir / "pair.toml")
    saved = toml.load(destination_dir / "pair.toml")
    assert saved["worlds"]["io"]["world"] == "../source/my_io.toml"
    assert saved["worlds"]["jupiter"]["world"] == "jupiter_simple"
    assert build_system(str(destination_dir / "pair.toml"))["io"].name == "io"


def test_an_inline_member_finds_its_data_file_beside_the_system_file(tmp_path):
    shutil.copy(_bundled_file("PREM.csv"), tmp_path / "beside.csv")
    world_config = toml.load(_bundled_file("earth_prem.toml"))
    world_config["data_file"] = "beside.csv"
    system_config = {"name": "solo", "worlds": {"earth": {"world": world_config}}}
    with open(tmp_path / "solo.toml", "w", encoding="utf-8") as system_file:
        toml.dump(system_config, system_file)
    system = build_system(str(tmp_path / "solo.toml"))
    assert len(system["earth"].get_config_dict()["layers"]) == 12


# =====================================================================================================================
# Prescribed heating
# =====================================================================================================================
_PRESCRIBED_POWER = 1.0e12        # [W]
_PRESCRIBED_RATE = 2.0e-12        # [W kg-1]


def _io_with_prescribed_heating():
    world = build_world("io")
    world.set_prescribed_heating("mantle", power=_PRESCRIBED_POWER)
    world.set_prescribed_heating("core", specific_rate=_PRESCRIBED_RATE)
    return world


def test_prescribed_heating_is_in_the_config_and_toml(tmp_path):
    world = _io_with_prescribed_heating()
    assert world.get_config_dict()["prescribed_heating"] == {
        "mantle": {"power_w": _PRESCRIBED_POWER}, "core": {"specific_rate_w_kg": _PRESCRIBED_RATE}}
    path = tmp_path / "io_heated.toml"
    world.save_to_toml(path)
    assert build_world(str(path)).prescribed_heating == world.prescribed_heating


def test_prescribed_heating_survives_a_binary_round_trip(tmp_path):
    world = _io_with_prescribed_heating()
    path = str(tmp_path / "io_heated.tpyb")
    world.save_binary(path)
    loaded = build_world("io")
    loaded.load_binary(path)
    assert loaded.prescribed_heating == world.prescribed_heating


@pytest.mark.parametrize("table, message", [
    ({"no_such_layer": {"power_w": 1.0}}, "names no layer"),
    ({"mantle": {"power_w": 1.0, "specific_rate_w_kg": 1.0}}, "exactly one of"),
    ({"mantle": {"power": 1.0}}, "exactly one of"),
    ({"mantle": {"power_w": "hot"}}, "finite number"),
], ids=["unknown_layer", "both", "unsuffixed_key", "not_a_number"])
def test_a_bad_prescribed_heating_table_is_refused(table, message):
    config = build_world("io").get_config_dict()
    config["prescribed_heating"] = table
    with pytest.raises(ValueError, match=message):
        build_world(config)
