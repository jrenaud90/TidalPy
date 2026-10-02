"""Configuration values TidalPy cannot use: a file's fall back to the defaults with a warning (so the import never
fails), an override passed to reinit raises; and the default file is written atomically."""
import os

import pytest

import TidalPy
from TidalPy import paths
from TidalPy.configurations import find_invalid_config_values, get_default_config, get_packaged_config


@pytest.fixture
def packaged():
    return get_packaged_config()


@pytest.mark.parametrize("override", [
    {"logging": {"console_level": "verbose"}},
    {"logging": {"file_level": 9}},
    {"logging": {"file_level": True}},
    {"radial_solver": {"rtol": -1.0}},
    {"eos_solver": {"max_iters": 2.5}},
    {"eos_solver": {"solve_temperature": 1}},
    {"numerical": {"minimum_frequency": 0.0}},
    {"numerical": {"love_solve_threads": -1}},
    {"worlds": {"star": {"albedo": "bright"}}},
    {"tides": {"fixed_q": 100.0}},
], ids=["level_name", "level_number", "level_bool", "negative_rtol", "float_iterations", "int_for_bool",
        "zero_frequency_floor", "negative_threads", "per_type_world_string", "scalar_for_list"])
def test_unusable_values_are_found(packaged, override):
    assert len(find_invalid_config_values(override, packaged)) == 1


@pytest.mark.parametrize("override", [
    {"logging": {"console_level": "WARNING", "file_level": 1}},
    {"radial_solver": {"rtol": 1}},
    {"tides": {"eccentricity_trunc_lvl": "exact", "obliquity_trunc_lvl": 2}},
    {"tides": {"star": {"fixed_q": [1.0e4]}}},
    {"layers": {"material": {"solid": {"eos": {"model": "constant"}}}}},
    {"numerical": {"love_solve_threads": 0}},
    {"not_a_section": {"x": 1}},
], ids=["levels", "int_for_float", "truncation_names", "per_type_tides", "material_table", "zero_threads",
        "unknown_keys_left_to_the_key_check"])
def test_usable_values_pass(packaged, override):
    assert find_invalid_config_values(override, packaged) == []


def test_reinit_refuses_an_unusable_override_and_changes_nothing():
    before = TidalPy.config["radial_solver"]["rtol"]
    with pytest.raises(ValueError, match="radial_solver.rtol"):
        TidalPy.reinit({"radial_solver": {"rtol": -1.0}})
    assert TidalPy.config["radial_solver"]["rtol"] == before


@pytest.fixture
def own_data_dir(tmp_path, monkeypatch):
    """A fresh data directory for this test only; the configuration is reloaded from the suite's afterward."""
    monkeypatch.setenv(paths.DATA_DIR_ENVIRONMENT_VARIABLE, str(tmp_path))
    yield tmp_path
    monkeypatch.undo()
    TidalPy.reinit("default")


def _config_path(data_dir):
    return os.path.join(str(data_dir), paths.get_data_version(), "Config", "TidalPy_Configs.toml")


def test_a_bad_value_in_the_file_falls_back_with_a_warning(own_data_dir):
    get_default_config()
    path = _config_path(own_data_dir)
    text = open(path, encoding="utf-8").read()
    assert 'console_level = "info"' in text
    open(path, "w", encoding="utf-8").write(text.replace('console_level = "info"', 'console_level = "verbose"'))
    with pytest.warns(UserWarning, match="logging.console_level"):
        config = get_default_config()
    assert config["logging"]["console_level"] == get_packaged_config()["logging"]["console_level"]


def test_a_file_that_does_not_parse_falls_back_with_a_warning(own_data_dir):
    get_default_config()
    with open(_config_path(own_data_dir), "a", encoding="utf-8") as config_file:
        config_file.write("\nthis is = = not toml\n")
    with pytest.warns(UserWarning, match="could not be read"):
        config = get_default_config()
    assert config["radial_solver"] == get_packaged_config()["radial_solver"]


def test_the_default_file_is_written_whole(own_data_dir):
    get_default_config()
    path = _config_path(own_data_dir)
    assert os.path.isfile(path)
    # Nothing is left behind from the atomic write.
    assert [name for name in os.listdir(os.path.dirname(path)) if name.endswith(".tmp")] == []


def test_write_file_atomically_keeps_an_existing_file(tmp_path):
    path = tmp_path / "file.txt"
    paths.write_file_atomically(str(path), b"first")
    paths.write_file_atomically(str(path), b"second", keep_existing=True)
    assert path.read_bytes() == b"first"
    paths.write_file_atomically(str(path), b"third")
    assert path.read_bytes() == b"third"


def test_a_relative_data_directory_is_made_absolute(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(paths.DATA_DIR_ENVIRONMENT_VARIABLE, "relative_data")
    assert paths.get_data_dir() == os.path.join(str(tmp_path / "relative_data"), paths.get_data_version())
