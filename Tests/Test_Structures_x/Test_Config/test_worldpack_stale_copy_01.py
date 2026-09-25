"""A data-directory copy of a bundled file that differs from the packaged one is reported, and never replaced."""
import os
import warnings

import pytest

import TidalPy
from TidalPy.structures_x import build_world
from TidalPy.structures_x.configs import worldpack


@pytest.fixture()
def data_dir(tmp_path, monkeypatch):
    """A private data directory holding fresh copies of the packaged files, and a clean warned-once record."""
    monkeypatch.setattr(worldpack, "get_worlds_x_dir", lambda: str(tmp_path))
    monkeypatch.setattr(worldpack, "_WARNED_STALE_COPIES", set())
    worldpack.install_worldpack_x()
    return tmp_path


def _edit(file_path):
    with open(file_path, "a", encoding="utf-8", newline="\n") as file:
        file.write("\n# A local edit.\n")


def _stale_warnings(record):
    return [entry for entry in record if "differs from the one packaged" in str(entry.message)]


def test_fresh_copies_do_not_warn(data_dir):
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        build_world("io")
        build_world("earth_prem")   # Also resolves its companion data file.
    assert not _stale_warnings(record)


def test_edited_world_warns_once_and_is_still_the_one_used(data_dir):
    copy_path = os.path.join(str(data_dir), "io.toml")
    _edit(copy_path)
    with open(copy_path, "rb") as file:
        edited = file.read()

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert worldpack.resolve_world_path("io") == copy_path
        build_world("io")
        build_world("io")
    found = _stale_warnings(record)
    assert len(found) == 1
    message = str(found[0].message)
    assert copy_path in message
    assert os.path.join(worldpack.PACKAGED_WORLDPACK_DIR, "io.toml") in message
    assert "install_worldpack_x(force=True)" in message
    assert "stale_worldpack_copy" in message

    # The copy is the user's: it is never rewritten.
    with open(copy_path, "rb") as file:
        assert file.read() == edited


def test_a_newline_only_difference_is_not_a_difference(data_dir):
    copy_path = os.path.join(str(data_dir), "io.toml")
    with open(copy_path, "rb") as file:
        contents = file.read()
    with open(copy_path, "wb") as file:
        file.write(contents.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
    assert worldpack.warn_if_stale_copy(copy_path) is False


def test_edited_data_file_warns(data_dir):
    _edit(os.path.join(str(data_dir), "PREM.csv"))
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        worldpack.resolve_data_file("PREM.csv")
    found = _stale_warnings(record)
    assert len(found) == 1 and "PREM.csv" in str(found[0].message)


def test_the_config_switch_silences_it(data_dir, monkeypatch):
    _edit(os.path.join(str(data_dir), "io.toml"))
    config_x = dict(TidalPy.config_x)
    config_x["warnings"] = {"stale_worldpack_copy": False}
    monkeypatch.setattr(TidalPy, "config_x", config_x)
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        build_world("io")
    assert not _stale_warnings(record)
    # The difference is still reported to a caller who asks.
    assert worldpack.warn_if_stale_copy(os.path.join(str(data_dir), "io.toml")) is True


def test_forced_install_replaces_the_copy_and_ends_the_warning(data_dir):
    copy_path = os.path.join(str(data_dir), "io.toml")
    _edit(copy_path)
    worldpack.install_worldpack_x(force=True)
    assert worldpack.warn_if_stale_copy(copy_path) is False


def test_the_switch_is_a_documented_default():
    assert TidalPy.config_x["warnings"]["stale_worldpack_copy"] is True


# =====================================================================================================================
# The other two switches of [warnings]
# =====================================================================================================================
def _switched_off(monkeypatch, name):
    config_x = dict(TidalPy.config_x)
    config_x["warnings"] = {name: False}
    monkeypatch.setattr(TidalPy, "config_x", config_x)


def test_schema_version_warning_has_a_switch(monkeypatch):
    from TidalPy.structures_x.configs import validate_schema_version
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert validate_schema_version({"name": "w"}) is True
    assert any("schema_version" in str(entry.message) for entry in record)

    _switched_off(monkeypatch, "schema_version")
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert validate_schema_version({"name": "w"}) is True
        assert validate_schema_version({"name": "w", "schema_version": "0.1.0"}) is True
    assert not record
    # A major difference is a refusal, not a warning, so the switch does not touch it.
    with pytest.raises(ValueError):
        validate_schema_version({"name": "w", "schema_version": "1.0.0"})


def test_truncation_promotion_warning_has_a_switch(monkeypatch):
    from TidalPy.structures_x.configs import world_builder
    monkeypatch.setattr(world_builder, "_WARNED_ECCENTRICITY_TRUNCATIONS", set())
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert world_builder._resolve_eccentricity_truncation(7) == 8
    assert any("not tabulated" in str(entry.message) for entry in record)

    monkeypatch.setattr(world_builder, "_WARNED_ECCENTRICITY_TRUNCATIONS", set())
    _switched_off(monkeypatch, "truncation_promotion")
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert world_builder._resolve_eccentricity_truncation(7) == 8
        assert world_builder._resolve_obliquity_truncation(3) == 4
    assert not record
