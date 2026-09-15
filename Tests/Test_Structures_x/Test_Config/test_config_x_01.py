"""
Tests for the new ``_x`` configuration (``TidalPy_Configs_x.toml`` / ``TidalPy.config_x``)
built by ``TidalPy.defaultc_x`` and loaded during initialization.

Covers that the ``_x`` config is loaded with the expected numerical and per-material
layer sections, that ``update_constants_x`` populates the shared C++ config
singleton (mirrored on the ``TidalPy.constants`` module globals) from it, that the
user's file is merged over the packaged defaults, and that a configuration can be
provided through ``TidalPy.reinit`` and saved (with its version header) for reproduction.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import copy
import importlib.metadata
import math

import pytest
import toml

import TidalPy
import TidalPy.configurations as configurations
from TidalPy.configurations import (
    config_version_header, get_default_config_x, get_packaged_config_x, merge_configs, save_config_x, set_config_x)
from TidalPy.exceptions import InitializationError
from TidalPy.structures_x.configs.toml_loader import MATERIAL_TYPES


@pytest.fixture
def restore_config_x():
    """Restore ``TidalPy.config_x`` and the C++ numerical settings after a test changes them."""
    from TidalPy.constants import update_constants_x
    original = copy.deepcopy(TidalPy.config_x)
    original_path = TidalPy._config_x_path
    yield
    TidalPy.config_x = original
    TidalPy._config_x_path = original_path
    update_constants_x()


@pytest.fixture
def isolated_config_dir(tmp_path, monkeypatch):
    """Point the TidalPy Config directory at a temporary folder so the user's own files are never read or written."""
    config_dir = tmp_path / "Config"
    config_dir.mkdir()
    monkeypatch.setattr(configurations, "get_config_dir", lambda: str(config_dir))
    return config_dir


def test_config_x_loaded():
    assert TidalPy.config_x is not None
    assert TidalPy.config_x["schema_version"] == "0.2.0"


def test_config_x_has_numerical_section():
    numerical = TidalPy.config_x["numerical"]
    for key in ("minimum_frequency", "maximum_frequency", "min_spin_orbit_diff",
                "minimum_viscosity", "minimum_modulus", "minimum_layer_thickness",
                "test_constant"):
        assert key in numerical


@pytest.mark.parametrize("material_type", list(MATERIAL_TYPES))
def test_config_x_has_each_material_block(material_type):
    layers = TidalPy.config_x["layers"]
    assert material_type in layers
    # Every material block carries an EOS default.
    assert "eos" in layers[material_type]


def test_mantle_rock_defaults_present():
    mantle = TidalPy.config_x["layers"]["mantle_rock"]
    assert mantle["shear_rheology"]["model"] == "andrade"
    assert mantle["shear_rheology"]["zeta"] == 1.0
    assert mantle["cooling"]["model"] == "convection"
    assert mantle["radiogenics"]["model"] == "isotope"


def test_update_constants_x_populated_singleton():
    # update_constants_x runs during initialization and feeds the _x numerical
    # settings into the constants module globals.
    from TidalPy import constants
    numerical = TidalPy.config_x["numerical"]
    assert math.isclose(constants.min_viscosity, numerical["minimum_viscosity"])
    assert math.isclose(constants.min_modulus, numerical["minimum_modulus"])
    assert math.isclose(constants.min_thickness, numerical["minimum_layer_thickness"])
    assert math.isclose(constants.test_constant, numerical["test_constant"])


# =====================================================================================================================
# Merging, providing, and saving a configuration
# =====================================================================================================================
def test_merge_configs_overrides_only_the_given_values():
    base = {"numerical": {"a": 1.0, "b": 2.0}, "tides": {"fixed_q": [100.0, 100.0]}}
    merged = merge_configs(base, {"numerical": {"b": 3.0}, "tides": {"fixed_q": [50.0]}, "extra": 4})
    # Tables merge key by key; a list replaces the base list whole.
    assert merged == {"numerical": {"a": 1.0, "b": 3.0}, "tides": {"fixed_q": [50.0]}, "extra": 4}
    # Neither input is modified.
    assert base == {"numerical": {"a": 1.0, "b": 2.0}, "tides": {"fixed_q": [100.0, 100.0]}}


def test_merge_configs_replaces_a_model_table_only_when_the_model_changes():
    base = {"layers": {"ice": {"shear_viscosity": {"model": "arrhenius", "arrhenius_coeff": 1.0, "stress_pa": 1.0}}}}

    # The same model (names compare case-insensitively) merges its parameters.
    same = merge_configs(base, {"layers": {"ice": {"shear_viscosity": {"model": "Arrhenius", "stress_pa": 2.0}}}})
    assert same["layers"]["ice"]["shear_viscosity"] == {"model": "Arrhenius", "arrhenius_coeff": 1.0, "stress_pa": 2.0}

    # Parameters without a model name merge into the default model.
    params = merge_configs(base, {"layers": {"ice": {"shear_viscosity": {"stress_pa": 5.0}}}})
    expected = {"model": "arrhenius", "arrhenius_coeff": 1.0, "stress_pa": 5.0}
    assert params["layers"]["ice"]["shear_viscosity"] == expected

    # A different model keeps none of the default model's parameters.
    other = {"model": "constant", "reference_viscosity_pas": 1.0e14}
    changed = merge_configs(base, {"layers": {"ice": {"shear_viscosity": other}}})
    assert changed["layers"]["ice"]["shear_viscosity"] == other


def test_version_header_lists_the_package_versions():
    header = config_version_header("A title")
    lines = header.rstrip("\n").split("\n")
    assert all(line.startswith("#") for line in lines)
    assert "#  A title" in lines
    assert f"#  TidalPy version: {TidalPy.version}" in lines
    assert f"#  SciPy version: {importlib.metadata.version('scipy')}" in lines
    assert f"#  CyRK version: {importlib.metadata.version('cyrk')}" in lines


def test_missing_user_file_is_written_with_the_defaults_and_a_version_header(isolated_config_dir, restore_config_x):
    config = get_default_config_x()
    written = (isolated_config_dir / "TidalPy_Configs_x.toml").read_text(encoding="utf-8")
    assert f"#  TidalPy version: {TidalPy.version}" in written.split("\n")[:6]
    assert toml.loads(written) == get_packaged_config_x()
    assert config == get_packaged_config_x()
    assert TidalPy.config_x == config


def test_partial_user_file_is_merged_over_the_packaged_defaults(isolated_config_dir, restore_config_x):
    (isolated_config_dir / "TidalPy_Configs_x.toml").write_text(
        "[numerical]\nminimum_viscosity = 250.0\n\n[layers.ice.shear_rheology]\nmodel = \"andrade\"\n",
        encoding="utf-8")
    # A hand-written file has no version header; the existing check warns and loading continues.
    with pytest.warns(UserWarning, match="Could not determine version"):
        config = get_default_config_x()
    packaged = get_packaged_config_x()
    assert config["numerical"]["minimum_viscosity"] == 250.0
    assert config["numerical"]["minimum_modulus"] == packaged["numerical"]["minimum_modulus"]
    assert config["tides"] == packaged["tides"]
    # The ice rheology switched model, so its table is exactly the user's.
    assert config["layers"]["ice"]["shear_rheology"] == {"model": "andrade"}
    assert config["layers"]["mantle_rock"] == packaged["layers"]["mantle_rock"]


def test_reinit_merges_a_provided_config_dict_and_keeps_it(restore_config_x):
    from TidalPy import constants
    original = copy.deepcopy(TidalPy.config_x)
    TidalPy.reinit(provided_config_x={"numerical": {"minimum_viscosity": 321.0}})
    assert TidalPy.config_x["numerical"]["minimum_viscosity"] == 321.0
    assert math.isclose(constants.min_viscosity, 321.0)
    assert TidalPy.config_x["numerical"]["minimum_modulus"] == original["numerical"]["minimum_modulus"]
    assert TidalPy.config_x["layers"] == original["layers"]

    # A later reinit without a configuration keeps the override.
    TidalPy.reinit()
    assert TidalPy.config_x["numerical"]["minimum_viscosity"] == 321.0


def test_saved_config_reproduces_the_settings_through_reinit(tmp_path, isolated_config_dir, restore_config_x):
    TidalPy.reinit(provided_config_x={"numerical": {"minimum_modulus": 2.5e-3}, "tides": {"max_degree_l": 3}})
    expected = copy.deepcopy(TidalPy.config_x)
    saved_path = save_config_x(str(tmp_path / "run_config.toml"))
    saved = (tmp_path / "run_config.toml").read_text(encoding="utf-8")
    assert saved.startswith("# ")
    assert f"#  CyRK version: {importlib.metadata.version('cyrk')}" in saved.split("\n")[:6]
    assert toml.loads(saved) == expected

    # "default" discards the overrides (the isolated Config directory holds only the packaged defaults).
    TidalPy.reinit(provided_config_x="default")
    assert TidalPy.config_x == get_packaged_config_x()

    # Loading the saved file restores every setting.
    TidalPy.reinit(provided_config_x=saved_path)
    assert TidalPy.config_x == expected

    # Saving again without overwriting picks a new file name.
    second_path = save_config_x(saved_path, overwrite=False)
    assert second_path != saved_path
    assert second_path.endswith(".toml")
    assert toml.load(second_path) == expected


def test_set_config_x_rejects_a_missing_file_and_other_types(tmp_path):
    with pytest.raises(InitializationError):
        set_config_x(str(tmp_path / "does_not_exist.toml"))
    with pytest.raises(TypeError):
        set_config_x(42)


def test_save_config_x_requires_a_toml_path(tmp_path):
    with pytest.raises(ValueError):
        save_config_x(str(tmp_path / "run_config.txt"))
