"""TidalPy's configuration: packaged defaults, the C++ constants fed from it, merging, reinit, and saving."""
import copy
import importlib.metadata
import math

import numpy as np
import pytest
import toml

import TidalPy
import TidalPy.configurations as configurations
from TidalPy.configurations import (
    config_version_header,
    get_default_config,
    get_packaged_config,
    merge_configs,
    plain_config,
    save_config,
    set_config,
)
from TidalPy.exceptions import InitializationError
from TidalPy.Material import available_materials, load_material
from TidalPy.Structures.configs.toml_loader import validate_layer_config


@pytest.fixture
def restore_config():
    """Restore ``TidalPy.config`` and the C++ numerical settings after a test changes them."""
    from TidalPy.constants import update_constants
    original = copy.deepcopy(TidalPy.config)
    original_path = TidalPy._config_path
    yield
    TidalPy.config = original
    TidalPy._config_path = original_path
    update_constants()


@pytest.fixture
def isolated_config_dir(tmp_path, monkeypatch):
    """Point the TidalPy Config directory at a temporary folder so the user's own files are never touched."""
    config_dir = tmp_path / "Config"
    config_dir.mkdir()
    monkeypatch.setattr(configurations, "get_config_dir", lambda: str(config_dir))
    return config_dir


# =====================================================================================================================
# Contents of the loaded configuration
# =====================================================================================================================
@pytest.mark.parametrize("section, keys", [
    ("pathing", {"save_directory", "append_datetime"}),
    ("logging", {"use_cwd", "write_log_to_disk", "file_level", "console_level", "print_log_notebook",
                 "write_log_notebook"}),
    ("configs", {"save_configs_locally", "use_cwd_for_config"}),
])
def test_config_package_section_has_exactly_its_keys(section, keys):
    assert set(TidalPy.config[section]) == keys


@pytest.mark.parametrize("section, keys", [
    ("numerical", ("minimum_frequency", "maximum_frequency", "minimum_modulus", "minimum_solid_rigidity",
                   "minimum_zone_fraction", "minimum_layer_thickness", "numerical_floor", "layer_continuity_rtol",
                   "max_start_radius_fraction",
                   "frequency_match_rtol", "minimum_nusselt", "maximum_eos_mass_ratio", "eos_invert_rtol",
                   "eos_invert_max_iters", "test_constant")),
    ("eos_solver", ("integration_method", "rtol", "atol", "pressure_tol", "max_iters", "nondimensionalize",
                    "slices_per_layer")),
    ("radial_solver", ("integration_method", "rtol", "atol", "use_kamata", "start_radius_tolerance", "scale_rtols",
                       "max_num_steps", "expected_size", "max_ram_mb", "nondimensionalize")),
])
def test_config_section_has_its_keys(section, keys):
    missing = [key for key in keys if key not in TidalPy.config[section]]
    assert not missing


_DEFAULT_VALUES = [
    (("schema_version",), "0.2.0"),
    (("radiogenics", "known_isotope_data"), {}),
    (("eos_solver", "nondimensionalize"), True),
    (("radial_solver", "nondimensionalize"), True),
    (("layers", "material"), "simple_rock"),
    (("radiogenics", "isotopes"), "modern_day_chondritic"),
] + [
    (("warnings", key), True)
    for key in ("stale_worldpack_copy", "schema_version", "truncation_promotion", "short_degree_list",
                "unknown_config_key")
]


@pytest.mark.parametrize("path, expected", [
    pytest.param(path, expected, id=".".join(path)) for path, expected in _DEFAULT_VALUES])
def test_config_default_value(path, expected):
    value = TidalPy.config
    for key in path:
        value = value[key]
    assert value == expected
    assert type(value) is type(expected)


def test_layers_table_names_only_the_default_material():
    # Packaged defaults, not TidalPy.config: a user file may name another material.
    assert get_packaged_config()["layers"] == {"material": "simple_rock"}
    # The material a layer that names none takes is a MatPack material that builds.
    name = TidalPy.config["layers"]["material"]
    assert name in available_materials()
    assert load_material(name).get_config_dict()["solid"]["eos"]["model"] == "constant"


@pytest.mark.parametrize("constant_name, numerical_key", [
    ("min_frequency", "minimum_frequency"),
    ("min_modulus", "minimum_modulus"),
    ("minimum_solid_rigidity", "minimum_solid_rigidity"),
    ("minimum_zone_fraction", "minimum_zone_fraction"),
    ("min_thickness", "minimum_layer_thickness"),
    ("frequency_match_rtol", "frequency_match_rtol"),
    ("minimum_nusselt", "minimum_nusselt"),
    ("maximum_eos_mass_ratio", "maximum_eos_mass_ratio"),
    ("eos_invert_rtol", "eos_invert_rtol"),
    ("eos_invert_max_iters", "eos_invert_max_iters"),
    ("test_constant", "test_constant"),
])
def test_update_constants_populated_singleton(constant_name, numerical_key):
    from TidalPy import constants
    value = getattr(constants, constant_name)
    expected = TidalPy.config["numerical"][numerical_key]
    if isinstance(expected, int):
        assert value == expected
    else:
        assert math.isclose(value, expected)


# =====================================================================================================================
# Merging
# =====================================================================================================================
def test_merge_configs_overrides_only_the_given_values():
    base = {"numerical": {"a": 1.0, "b": 2.0}, "tides": {"fixed_q": [100.0, 100.0]}}
    merged = merge_configs(base, {"numerical": {"b": 3.0}, "tides": {"fixed_q": [50.0]}, "extra": 4})
    # Tables merge key by key; a list replaces the base list whole.
    assert merged == {"numerical": {"a": 1.0, "b": 3.0}, "tides": {"fixed_q": [50.0]}, "extra": 4}
    assert base == {"numerical": {"a": 1.0, "b": 2.0}, "tides": {"fixed_q": [100.0, 100.0]}}


def _model_table(table):
    return {"section": {"shear_viscosity": table}}


_ARRHENIUS = {"model": "arrhenius", "arrhenius_coeff": 1.0, "stress_pa": 1.0}
_CONSTANT_VISCOSITY = {"model": "constant", "reference_viscosity_pas": 1.0e14}


@pytest.mark.parametrize("base, override, expected", [
    pytest.param(
        _model_table(_ARRHENIUS),
        _model_table({"model": "Arrhenius", "stress_pa": 2.0}),
        _model_table({"model": "Arrhenius", "arrhenius_coeff": 1.0, "stress_pa": 2.0}),
        id="same-model-any-case"),
    pytest.param(
        _model_table(_ARRHENIUS),
        _model_table({"stress_pa": 5.0}),
        _model_table({"model": "arrhenius", "arrhenius_coeff": 1.0, "stress_pa": 5.0}),
        id="no-model-name"),
    pytest.param(
        _model_table(_ARRHENIUS),
        _model_table(_CONSTANT_VISCOSITY),
        # merge_configs has no model rule: a base key the new model does not read carries over.
        _model_table({**_ARRHENIUS, **_CONSTANT_VISCOSITY}),
        id="different-model"),
    pytest.param(
        {"radiogenics": {"model": "isotope", "isotopes": "modern_day_chondritic", "ref_time_s": 1.4e17}},
        {"radiogenics": {"model": "fixed", "fixed_heat_production_w_kg": 1.0e-11}},
        {"radiogenics": {"model": "fixed", "isotopes": "modern_day_chondritic", "ref_time_s": 1.4e17,
                         "fixed_heat_production_w_kg": 1.0e-11}},
        id="radiogenics-different-model"),
    pytest.param(
        {"radiogenics": {"model": "isotope", "isotopes": "modern_day_chondritic"}},
        {"radiogenics": {"model": "isotope", "ref_time_s": 0.0}},
        {"radiogenics": {"model": "isotope", "isotopes": "modern_day_chondritic", "ref_time_s": 0.0}},
        id="radiogenics-same-model"),
    # Nested tables merge key by key at every depth.
    pytest.param(
        {"material": {"solid": {"eos": {"model": "constant", "reference_density_kg_m3": 917.0},
                                "shear_viscosity": {"model": "arrhenius", "arrhenius_coeff": 1.0e14}}}},
        {"material": {"solid": {"eos": {"model": "birch_murnaghan", "reference_bulk_modulus_pa": 1.0e10}}}},
        {"material": {"solid": {"eos": {"model": "birch_murnaghan", "reference_density_kg_m3": 917.0,
                                        "reference_bulk_modulus_pa": 1.0e10},
                                "shear_viscosity": {"model": "arrhenius", "arrhenius_coeff": 1.0e14}}}},
        id="nested-model-change"),
])
def test_merge_configs_merges_model_tables_key_by_key(base, override, expected):
    assert merge_configs(base, override) == expected


@pytest.mark.parametrize("value", ("false", 0, 1.0))
def test_switches_must_be_booleans(value):
    with pytest.raises(ValueError, match="must be true or false"):
        validate_layer_config("mantle", {"radius_fraction": 1.0, "use_tides": value})


def test_numpy_values_are_written_as_numbers():
    config = {
        "a": np.float64(4.2e8),
        "b": [np.int64(3), np.int64(2)],
        "c": np.array([1.0, 2.0]),
        "d": {"e": np.bool_(True)}}
    loaded = toml.loads(toml.dumps(plain_config(config)))
    assert loaded == {"a": 4.2e8, "b": [3, 2], "c": [1.0, 2.0], "d": {"e": True}}


# =====================================================================================================================
# Providing and saving a configuration
# =====================================================================================================================
def test_version_header_lists_the_package_versions():
    lines = config_version_header("A title").rstrip("\n").split("\n")
    assert all(line.startswith("#") for line in lines)
    assert "#  A title" in lines
    assert f"#  TidalPy version: {TidalPy.version}" in lines
    assert f"#  SciPy version: {importlib.metadata.version('scipy')}" in lines
    assert f"#  CyRK version: {importlib.metadata.version('cyrk')}" in lines


def test_missing_user_file_is_written_with_the_defaults_and_a_version_header(isolated_config_dir, restore_config):
    config = get_default_config()
    written = (isolated_config_dir / "TidalPy_Configs.toml").read_text(encoding="utf-8")
    assert f"#  TidalPy version: {TidalPy.version}" in written.split("\n")[:6]
    assert toml.loads(written) == get_packaged_config()
    assert config == get_packaged_config()
    assert TidalPy.config == config


def test_partial_user_file_is_merged_over_the_packaged_defaults(isolated_config_dir, restore_config):
    (isolated_config_dir / "TidalPy_Configs.toml").write_text(
        "[numerical]\nmaximum_eos_mass_ratio = 25.0\n\n[layers]\nmaterial = \"simple_ice\"\n\n"
        "[layers.ice.shear_rheology]\nmodel = \"andrade\"\n",
        encoding="utf-8")
    # A hand-written file has no version header, so loading warns and continues; the retired [layers.ice] table of an
    # earlier build is dropped, with its own warning.
    with pytest.warns(UserWarning) as record:
        config = get_default_config()
    messages = [str(entry.message) for entry in record]
    assert any("Could not determine version" in message for message in messages)
    assert sum("no longer read" in message and "layers.ice" in message for message in messages) == 1
    packaged = get_packaged_config()
    assert config["numerical"]["maximum_eos_mass_ratio"] == 25.0
    assert config["numerical"]["minimum_modulus"] == packaged["numerical"]["minimum_modulus"]
    assert config["tides"] == packaged["tides"]
    assert config["layers"] == {"material": "simple_ice"}


def test_reinit_merges_a_provided_config_dict_and_keeps_it(restore_config):
    from TidalPy import constants
    original = copy.deepcopy(TidalPy.config)
    TidalPy.reinit(provided_config={"numerical": {"maximum_eos_mass_ratio": 32.0}})
    assert TidalPy.config["numerical"]["maximum_eos_mass_ratio"] == 32.0
    assert math.isclose(constants.maximum_eos_mass_ratio, 32.0)
    assert TidalPy.config["numerical"]["minimum_modulus"] == original["numerical"]["minimum_modulus"]
    assert TidalPy.config["layers"] == original["layers"]

    TidalPy.reinit()
    assert TidalPy.config["numerical"]["maximum_eos_mass_ratio"] == 32.0


def test_saved_config_reproduces_the_settings_through_reinit(tmp_path, isolated_config_dir, restore_config):
    TidalPy.reinit(provided_config={"numerical": {"minimum_modulus": 2.5e-3}, "tides": {"max_degree_l": 3}})
    expected = copy.deepcopy(TidalPy.config)
    saved_path = save_config(str(tmp_path / "run_config.toml"))
    saved = (tmp_path / "run_config.toml").read_text(encoding="utf-8")
    assert saved.startswith("# ")
    assert f"#  CyRK version: {importlib.metadata.version('cyrk')}" in saved.split("\n")[:6]
    assert toml.loads(saved) == expected

    # The isolated Config directory holds only the packaged defaults, so "default" discards the overrides.
    TidalPy.reinit(provided_config="default")
    assert TidalPy.config == get_packaged_config()

    TidalPy.reinit(provided_config=saved_path)
    assert TidalPy.config == expected

    # Saving again without overwriting picks a new file name.
    second_path = save_config(saved_path, overwrite=False)
    assert second_path != saved_path
    assert second_path.endswith(".toml")
    assert toml.load(second_path) == expected


@pytest.mark.parametrize("make_source, error", [
    pytest.param(lambda folder: str(folder / "does_not_exist.toml"), InitializationError, id="missing-file"),
    pytest.param(lambda folder: 42, TypeError, id="not-a-path"),
])
def test_set_config_rejects_a_bad_source(tmp_path, make_source, error):
    with pytest.raises(error):
        set_config(make_source(tmp_path))


def test_save_config_requires_a_toml_path(tmp_path):
    with pytest.raises(ValueError):
        save_config(str(tmp_path / "run_config.txt"))


# =====================================================================================================================
# Solver defaults flowing from the configuration
# =====================================================================================================================
def _homogeneous_inputs():
    """A homogeneous Maxwell sphere as standalone radial_solver inputs."""
    from TidalPy.Rheology import Maxwell
    frequency = 2.0 * math.pi / 86400.0
    num_slices = 20
    radius = np.linspace(0.0, 6.0e6, num_slices)
    density = 3500.0 * np.ones(num_slices)
    bulk = 1.0e11 * np.ones(num_slices, dtype=np.complex128)
    shear = Maxwell().calc_complex_modulus_vectorize_modulus(
        5.0e10 * np.ones(num_slices), 1.0e20 * np.ones(num_slices), frequency)
    return (radius, density, bulk, shear, frequency, 3500.0, ("solid",), (False,), (False,), np.array([6.0e6]))


def test_radial_solver_defaults_follow_the_config(restore_config):
    """The standalone solver takes its integration settings from [radial_solver]; an explicit argument wins."""
    from TidalPy.RadialSolver import radial_solver
    args = _homogeneous_inputs()
    config_rtol = TidalPy.config["radial_solver"]["rtol"]
    config_atol = TidalPy.config["radial_solver"]["atol"]
    default = radial_solver(*args)
    explicit = radial_solver(*args, integration_rtol=config_rtol, integration_atol=config_atol)
    assert default.success and explicit.success
    assert complex(default.k) == complex(explicit.k)
    assert default.steps_taken.max() == explicit.steps_taken.max()

    TidalPy.reinit(provided_config={"radial_solver": {"rtol": 1.0e-10, "atol": 1.0e-13}})
    tight = radial_solver(*args)
    assert tight.success
    assert tight.steps_taken.max() > default.steps_taken.max()
    assert math.isclose(complex(tight.k).real, complex(default.k).real, rel_tol=1.0e-4)


def test_world_love_defaults_follow_the_config(restore_config):
    """The world Love solve starts from [radial_solver] and the world's [tides] Love method."""
    from TidalPy.Structures import build_world
    world = build_world("earth_simple")
    world.solve_eos()
    frequency = 2.0 * math.pi / 86400.0
    config = TidalPy.config["radial_solver"]
    world.solve_love_numbers(frequency=frequency)
    k_default = world.love_number_k
    world.solve_love_numbers(
        frequency=frequency,
        rtol=config["rtol"],
        atol=config["atol"],
        use_kamata=config["use_kamata"],
        scale_rtols=config["scale_rtols"],
        start_radius_tol=config["start_radius_tolerance"],
        integration_method=config["integration_method"],
    )
    assert world.love_number_k == k_default

    TidalPy.reinit(provided_config={"radial_solver": {"rtol": 1.0e-10, "atol": 1.0e-13}})
    world.solve_love_numbers(frequency=frequency)
    k_tight = world.love_number_k
    assert k_tight != k_default
    assert math.isclose(k_tight.real, k_default.real, rel_tol=1.0e-4)

    world.set_tide_config(love_method="homogeneous")
    world.solve_love_numbers(frequency=frequency)
    assert world.love_method == "homogeneous"
    world.solve_love_numbers(frequency=frequency, love_method="radial_solver")
    assert world.love_method == "radial_solver"
