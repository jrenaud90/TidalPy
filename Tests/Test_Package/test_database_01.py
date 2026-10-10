"""The pack database (TidalPy.database): the material, world, and system TOML files are read into memory at import, a
lookup by name comes from memory while the file is unchanged, and an edit, a new file, a deleted copy, or a new data
directory is noticed by the next lookup."""

import os

import pytest
import toml

from TidalPy import database
from TidalPy.database import MAT_PACK, WORLD_PACK
from TidalPy.Material import load_material, material_config, matpack
from TidalPy.schema import SCHEMA_VERSION
from TidalPy.Structures import build_world
from TidalPy.Structures.configs import toml_loader, worldpack
from TidalPy.Utilities.data_pack import tomllib


@pytest.fixture()
def materials_dir(tmp_path, monkeypatch):
    """A private materials data directory holding fresh copies of the packaged files."""
    monkeypatch.setattr(matpack, "get_materials_dir", lambda: str(tmp_path))
    monkeypatch.setattr(MAT_PACK, "p_warned_stale_copies", set())
    matpack.install_matpack()
    return tmp_path


@pytest.fixture()
def worlds_dir(tmp_path, monkeypatch):
    """A private worlds data directory holding fresh copies of the packaged files."""
    monkeypatch.setattr(worldpack, "get_worlds_dir", lambda: str(tmp_path))
    monkeypatch.setattr(WORLD_PACK, "p_warned_stale_copies", set())
    worldpack.install_worldpack()
    return tmp_path


def _rewrite(path, edit):
    """Rewrite a TOML file through ``edit(table)``, moving its modification time on by a second so the change is seen
    whatever the file system's time resolution."""
    status = os.stat(path)
    with open(path, "r", encoding="utf-8") as file:
        table = toml.load(file)
    edit(table)
    with open(path, "w", encoding="utf-8", newline="\n") as file:
        toml.dump(table, file)
    os.utime(path, ns=(status.st_atime_ns, status.st_mtime_ns + 1_000_000_000))


def _density(name):
    return material_config(name)["solid"]["eos"]["reference_density_kg_m3"]


def _set_density(value):
    def edit(table):
        table["solid"]["eos"]["reference_density_kg_m3"] = value
    return edit


# =====================================================================================================================
# Import
# =====================================================================================================================
@pytest.mark.parametrize("pack", [MAT_PACK, WORLD_PACK], ids=["materials", "worlds"])
def test_import_reads_every_packaged_toml(pack):
    packaged = {name.lower() for name in os.listdir(pack.packaged_dir) if name.lower().endswith(".toml")}
    # A lookup first, in case another test left the pack pointed at a directory since removed.
    assert pack.entry(sorted(packaged)[0]) is not None
    assert packaged <= set(pack.p_entries)
    assert not any(name.endswith(".csv") for name in pack.p_entries)
    if tomllib is not None:
        assert all(entry.p_table is not None for entry in pack.p_entries.values())


def test_a_lookup_by_name_comes_from_memory():
    first = MAT_PACK.entry("simple_rock.toml")
    assert MAT_PACK.entry("SIMPLE_ROCK.toml") is first
    assert matpack.p_resolve_preset("simple_rock", ()) is matpack.p_resolve_preset("simple_rock", ())
    assert WORLD_PACK.entry("io.toml") is WORLD_PACK.entry("io.toml")


def test_callers_get_copies_the_database_keeps_unchanged():
    table = material_config("simple_rock")
    table["solid"]["eos"]["reference_density_kg_m3"] = 1.0
    assert _density("simple_rock") == 3300.0
    config = toml_loader.load_toml(worldpack.resolve_world_path("io"))
    config["radius_m"] = 1.0
    assert WORLD_PACK.entry("io.toml").table()["radius_m"] == 1821490.0


def test_load_database_never_raises_for_a_broken_file(materials_dir):
    with open(os.path.join(str(materials_dir), "broken.toml"), "w", encoding="utf-8", newline="\n") as file:
        file.write("this is = = not toml\n")
    database.load_database()
    with pytest.raises(ValueError, match="could not parse"):
        load_material("broken")
    with pytest.raises(ValueError, match="could not parse"):
        load_material("broken")


# =====================================================================================================================
# Changes on disk
# =====================================================================================================================
def test_an_edited_material_is_read_again(materials_dir):
    assert _density("simple_rock") == 3300.0
    _rewrite(os.path.join(str(materials_dir), "simple_rock.toml"), _set_density(3000.0))
    assert _density("simple_rock") == 3000.0
    assert load_material("simple_rock").get_config_dict()["solid"]["eos"]["reference_density_kg_m3"] == 3000.0


def test_an_edited_preset_resolves_its_users_again(materials_dir):
    with open(os.path.join(str(materials_dir), "my_rock.toml"), "w", encoding="utf-8", newline="\n") as file:
        toml.dump({"schema_version": SCHEMA_VERSION, "description": "Mine.", "category": "rocky",
                   "preset": "simple_rock"}, file)
    assert _density("my_rock") == 3300.0
    _rewrite(os.path.join(str(materials_dir), "simple_rock.toml"), _set_density(2900.0))
    assert _density("my_rock") == 2900.0


def test_a_new_file_is_found_and_a_deleted_copy_falls_back_to_the_package(materials_dir):
    with pytest.raises(ValueError, match="no MatPack material named 'late_rock'"):
        load_material("late_rock")
    with open(os.path.join(str(materials_dir), "late_rock.toml"), "w", encoding="utf-8", newline="\n") as file:
        toml.dump({"schema_version": SCHEMA_VERSION, "description": "Late.", "category": "rocky",
                   "preset": "simple_rock", "solid": {"eos": {"reference_density_kg_m3": 3100.0}}}, file)
    assert _density("late_rock") == 3100.0

    _rewrite(os.path.join(str(materials_dir), "simple_rock.toml"), _set_density(3000.0))
    assert _density("simple_rock") == 3000.0
    os.remove(os.path.join(str(materials_dir), "simple_rock.toml"))
    assert _density("simple_rock") == 3300.0
    assert os.path.dirname(MAT_PACK.find("simple_rock.toml")) == matpack.PACKAGED_MATPACK_DIR


def test_a_new_data_directory_reads_the_database_again(materials_dir):
    assert os.path.dirname(MAT_PACK.find("water.toml")) == str(materials_dir)
    assert all(os.path.dirname(entry.path) == str(materials_dir) for entry in MAT_PACK.p_entries.values())


def test_a_forced_install_reads_the_database_again(materials_dir):
    _rewrite(os.path.join(str(materials_dir), "simple_rock.toml"), _set_density(3000.0))
    assert _density("simple_rock") == 3000.0
    matpack.install_matpack(force=True)
    assert _density("simple_rock") == 3300.0


def test_an_edited_world_is_built_again(worlds_dir):
    assert build_world("io").source_config["radius_m"] == 1821490.0

    def edit(table):
        table["radius_m"] = 1821000.0
    _rewrite(os.path.join(str(worlds_dir), "io.toml"), edit)
    assert build_world("io").source_config["radius_m"] == 1821000.0
