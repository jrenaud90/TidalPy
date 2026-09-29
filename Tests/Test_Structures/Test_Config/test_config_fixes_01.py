"""Layer defaults across a model change, the radial data-file reader's checks, and running without a data directory."""
import copy
import os
import subprocess
import sys
import textwrap
import warnings

import numpy as np
import pytest

import TidalPy
from TidalPy import paths
from TidalPy import configurations
from TidalPy.configurations import get_packaged_config, merge_configs
from TidalPy.PartialMelt.partial_melt import make_partial_melt
from TidalPy.Structures import build_world, build_world_from_dict, available_worlds
from TidalPy.Structures.configs import data_file, world_builder, worldpack


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
def fresh_data_dir_warning(monkeypatch):
    """Forget which data directories were already warned about, so each test sees its own warning."""
    monkeypatch.setattr(paths, "_WARNED_UNUSABLE_DATA_DIRS", set())


def _one_layer_world(layer_cfg):
    layer = {"class": "solidliquid", "layer_index": 0, "radius_outer_m": 1.0e6, **layer_cfg}
    return {"schema_version": "0.2.0", "name": "x", "type": "terrestrial", "radius_m": 1.0e6, "mass_kg": 1.0e22,
            "layers": {"shell": layer}}


def _two_layer_profile(**layer_tables):
    """A solid layer under a liquid one, with no row repeated at the boundary between them."""
    mapping = {
        "radius_km": [0.0, 500.0, 1000.0, 1100.0, 1500.0, 2000.0],
        "density":   [5000.0, 4900.0, 4800.0, 1100.0, 1050.0, 1000.0],
        "vp":        [8000.0, 7900.0, 7800.0, 1500.0, 1500.0, 1500.0],
        "vs":        [4000.0, 3900.0, 3800.0, 0.0, 0.0, 0.0],
    }
    config = {"schema_version": "0.2.0", "name": "two-layer", "type": "terrestrial", "radius_m": 2.0e6,
              "mass_kg": 5.0e22, "data": mapping}
    if layer_tables:
        config["layers"] = layer_tables
    return config


# =====================================================================================================================
# A model table fills in from the layer's material type, whichever model it names
# =====================================================================================================================
@pytest.mark.parametrize("material_type, solidus, liquidus, liquid_viscosity", [
    ("iron", 4000.0, 5000.0, 1.3e-2),
    ("ice", 250.0, 273.15, 8.9e-4),
    ("hp_ice", 270.0, 300.0, 8.9e-4),
    ("mantle_rock", 1600.0, 2000.0, 0.2),
])
def test_switching_the_melt_model_keeps_the_types_melt_properties(material_type, solidus, liquidus, liquid_viscosity):
    # Every type but mantle_rock defaults to the 'off' melt model; the henning table takes the type's values.
    world = build_world_from_dict(_one_layer_world(
        {"type": material_type, "material": {"partial_melt": {"model": "henning"}}}))
    melt = world.shell.get_config_dict()["material"]["partial_melt"]
    assert melt["model"] == "henning"
    assert melt["solidus_k"] == solidus
    assert melt["liquidus_k"] == liquidus
    assert melt["liquid_viscosity_pas"] == liquid_viscosity


def test_a_users_key_wins_over_the_types_default():
    world = build_world_from_dict(_one_layer_world(
        {"type": "ice", "material": {"partial_melt": {"model": "spohn", "solidus_k": 260.0}}}))
    melt = world.shell.get_config_dict()["material"]["partial_melt"]
    assert melt["solidus_k"] == 260.0
    assert melt["liquidus_k"] == 273.15


def test_switching_the_viscosity_model_keeps_the_types_activation_energy():
    table = {"model": "reference", "reference_viscosity_pas": 1.0e14, "reference_temperature_k": 260.0}
    world = build_world_from_dict(_one_layer_world({"type": "ice", "material": {"shear_viscosity": table}}))
    viscosity = world.shell.get_config_dict()["material"]["shear_viscosity"]
    assert viscosity["reference_viscosity_pas"] == 1.0e14
    # The ice value, not the silicate one (3e5).
    assert viscosity["molar_activation_energy_j_mol"] == 59.4e3


def test_a_layer_with_no_material_type_takes_only_its_own_keys():
    world = build_world_from_dict(_one_layer_world({
        "type": "none",
        "material": {"model": "constant", "reference_density_kg_m3": 1000.0, "partial_melt": {"model": "henning"}}}))
    melt = world.shell.get_config_dict()["material"]["partial_melt"]
    assert melt["solidus_k"] == make_partial_melt("henning", {}).get_config_dict()["solidus_k"]


def test_a_type_named_on_a_profile_layer_supplies_the_defaults():
    config = _two_layer_profile(layer_1={"layer_index": 1, "type": "ice",
                                         "material": {"partial_melt": {"model": "henning"}}})
    melt = build_world(config).layer_1.get_config_dict()["material"]["partial_melt"]
    assert melt["solidus_k"] == 250.0
    assert melt["liquid_density_kg_m3"] == 999.84


def test_merge_configs_keeps_the_packaged_melt_properties_across_a_model_change():
    packaged = get_packaged_config()
    merged = merge_configs(packaged, {"layers": {"ice": {"material": {"partial_melt": {"model": "spohn"}}}}})
    expected = dict(packaged["layers"]["ice"]["material"]["partial_melt"], model="spohn")
    assert merged["layers"]["ice"]["material"]["partial_melt"] == expected


@pytest.mark.parametrize("override_model, keeps_reference_time", [
    ("fixed", False),
    # An alias of the default model is the same model.
    ("isotopes", True),
])
def test_a_model_specific_key_is_dropped_only_on_a_real_model_change(override_model, keeps_reference_time):
    assert "ref_time_s" in configurations.MODEL_SPECIFIC_KEYS["radiogenics"]
    defaults = {"model": "isotope", "isotopes": "modern_day_chondritic", "ref_time_s": 1.0e16}
    merged = world_builder._merge_section(defaults, {"model": override_model}, "radiogenics")
    assert merged["model"] == override_model
    assert merged["isotopes"] == "modern_day_chondritic"
    assert ("ref_time_s" in merged) == keeps_reference_time


# =====================================================================================================================
# The radial data-file reader
# =====================================================================================================================
def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_a_header_row_that_does_not_fit_the_data_is_refused(tmp_path):
    path = _write(tmp_path, "wide_header.txt",
                  "radius density vp vs extra\n"
                  "0 9000 10000 0\n"
                  "1000 8000 10000 3000\n")
    with pytest.raises(ValueError, match=r"naming 5 column\(s\).*rows of 4"):
        data_file.load_radial_data(path)


def test_a_whitespace_header_keeps_a_bracketed_unit_with_its_name(tmp_path):
    path = _write(tmp_path, "depth_profile.txt",
                  "depth (km)  density  vp  vs\n"
                  "0     3000  8000  4000\n"
                  "1000  4000  9000  4500\n"
                  "2000  5000  10000  5000\n")
    arrays = data_file.load_radial_data(path, surface_radius=2.0e6)
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 1.0e6, 2.0e6])
    # The surface (depth 0) is the outermost row, not the innermost.
    np.testing.assert_allclose(arrays["density_kg_m3"], [5000.0, 4000.0, 3000.0])


def test_an_eta_column_is_not_a_viscosity(tmp_path):
    # In PREM and the IRIS tables eta is the anisotropy parameter.
    path = _write(tmp_path, "anisotropic.csv",
                  "radius_km,density,vp,vs,eta\n"
                  "0,9000,10000,3000,1.0\n"
                  "1000,8000,10000,3000,0.9\n")
    arrays = data_file.load_radial_data(path)
    assert arrays["shear_viscosity_pas"] is None
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 1.0e6])


@pytest.mark.parametrize("names", [("shear_modulus", "bulk_modulus"), ("mu", "k")])
def test_moduli_in_gpa_without_their_unit_are_refused(names):
    shear_name, bulk_name = names
    profile = {"radius_km": [0.0, 1000.0, 2000.0], "density": [5000.0, 4000.0, 3000.0],
               shear_name: [60.0, 60.0, 0.0], bulk_name: [100.0, 100.0, 100.0]}
    with pytest.raises(ValueError, match="bulk_modulus_gpa"):
        data_file.load_radial_data(profile)


def test_moduli_with_their_unit_are_converted():
    arrays = data_file.load_radial_data({
        "radius_km": [0.0, 1000.0, 2000.0], "density": [5000.0, 4000.0, 3000.0],
        "shear_modulus_gpa": [60.0, 60.0, 0.0], "bulk_modulus_gpa": [100.0, 100.0, 100.0]})
    # A liquid's zero shear modulus is not a unit mistake.
    np.testing.assert_allclose(arrays["shear_modulus_pa"], [6.0e10, 6.0e10, 0.0])
    np.testing.assert_allclose(arrays["bulk_modulus_pa"], 1.0e11)


def _interpolated_layer(radius):
    return {"type": "none", "material": {
        "model": "interpolate", "radius_m": radius, "density_kg_m3": [3000.0, 3000.0, 3000.0],
        "shear_modulus_pa": [5.0e10, 5.0e10, 5.0e10], "bulk_modulus_pa": [1.0e11, 1.0e11, 1.0e11]}}


def test_an_interpolated_table_that_misses_its_layer_is_refused():
    # A radius given in km under radius_m would otherwise build a uniform layer.
    with pytest.raises(ValueError, match="must cover its layer"):
        build_world_from_dict(_one_layer_world(_interpolated_layer([0.0, 500.0, 1000.0])))


def test_an_interpolated_table_that_covers_its_layer_builds():
    world = build_world_from_dict(_one_layer_world(_interpolated_layer([0.0, 5.0e5, 1.0e6])))
    assert world.shell.get_config_dict()["material"]["radius_m"][-1] == 1.0e6


def test_a_profile_without_a_repeated_boundary_row_spans_each_layer():
    world = build_world(_two_layer_profile())
    upper = world.get_config_dict()["layers"]["layer_1"]["material"]
    # The liquid layer's first row is repeated down at the solid layer's top.
    assert upper["radius_m"][:2] == [1.0e6, 1.1e6]
    assert upper["density_kg_m3"][:2] == [1100.0, 1100.0]
    # The saved configuration rebuilds.
    rebuilt = build_world_from_dict(world.get_config_dict())
    assert rebuilt.get_config_dict()["layers"]["layer_1"]["material"]["radius_m"] == upper["radius_m"]


# =====================================================================================================================
# Running without a data directory
# =====================================================================================================================
def test_the_data_directory_follows_the_environment_variable(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.DATA_DIR_ENVIRONMENT_VARIABLE, str(tmp_path))
    assert paths.get_data_dir() == os.path.join(str(tmp_path), paths.get_data_version())
    assert paths.get_config_dir() == os.path.join(paths.get_data_dir(), "Config")
    assert os.path.isdir(paths.get_config_dir())


def test_a_data_directory_that_cannot_be_created_is_none(tmp_path, monkeypatch, fresh_data_dir_warning):
    blocker = tmp_path / "blocker"
    blocker.write_text("")
    monkeypatch.setenv(paths.DATA_DIR_ENVIRONMENT_VARIABLE, str(blocker / "below_a_file"))
    with pytest.warns(UserWarning, match="cannot use its data directory"):
        assert paths.get_config_dir() is None
    # Once per session.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert paths.get_worlds_dir() is None


def test_an_unwritable_config_file_falls_back_to_the_packaged_defaults(
        tmp_path, monkeypatch, restore_config, fresh_data_dir_warning):
    blocker = tmp_path / "blocker"
    blocker.write_text("")
    # Writing a file below a file fails, as it does in a read-only directory.
    monkeypatch.setattr(configurations, "get_config_dir", lambda: str(blocker))
    with pytest.warns(UserWarning, match="cannot use its data directory"):
        config = configurations.get_default_config()
    assert config == get_packaged_config()
    assert TidalPy._config_path is None


def test_the_world_pack_is_read_from_the_package_without_a_data_directory(monkeypatch):
    monkeypatch.setattr(worldpack, "get_worlds_dir", lambda: None)
    assert worldpack.install_worldpack() is None
    assert os.path.dirname(worldpack.resolve_world_path("io")) == worldpack.PACKAGED_WORLDPACK_DIR
    assert os.path.dirname(worldpack.resolve_data_file("PREM.csv")) == worldpack.PACKAGED_WORLDPACK_DIR
    assert "io" in available_worlds()


def test_a_world_pack_directory_that_cannot_be_written_is_skipped(tmp_path, monkeypatch, fresh_data_dir_warning):
    blocker = tmp_path / "blocker"
    blocker.write_text("")
    monkeypatch.setattr(worldpack, "get_worlds_dir", lambda: str(blocker))
    with pytest.warns(UserWarning, match="cannot use its data directory"):
        path = worldpack.resolve_world_path("io")
    assert os.path.dirname(path) == worldpack.PACKAGED_WORLDPACK_DIR


def _import_in_subprocess(tmp_path, data_dir, script):
    environment = dict(os.environ)
    environment[paths.DATA_DIR_ENVIRONMENT_VARIABLE] = str(data_dir)
    # A neutral working directory, so the repository's in-tree files are not imported.
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(script)],
        cwd=str(tmp_path),
        env=environment,
        capture_output=True,
        text=True,
        timeout=300)


def test_import_uses_the_environment_variables_data_directory(tmp_path):
    result = _import_in_subprocess(tmp_path, tmp_path / "data", """
        import TidalPy
        print(TidalPy._config_path)
    """)
    assert result.returncode == 0, result.stderr
    config_path = result.stdout.strip().splitlines()[-1]
    assert config_path.startswith(str(tmp_path / "data"))
    assert os.path.isfile(config_path)


def test_import_succeeds_without_a_usable_data_directory(tmp_path):
    blocker = tmp_path / "blocker"
    blocker.write_text("")
    result = _import_in_subprocess(tmp_path, blocker / "below_a_file", """
        import warnings
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            import TidalPy
            from TidalPy.Structures import build_world
            world = build_world("io")
        print(sum("cannot use its data directory" in str(item.message) for item in caught))
        print(TidalPy._config_path, world.name)
    """)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[-2] == "1"
    assert lines[-1] == "None Io"


# =====================================================================================================================
# [configs] use_cwd_for_config
# =====================================================================================================================
def test_use_cwd_for_config_without_a_file_keeps_the_loaded_config(tmp_path, monkeypatch, restore_config):
    monkeypatch.chdir(tmp_path)
    TidalPy.config["configs"]["use_cwd_for_config"] = True
    rtol = TidalPy.config["radial_solver"]["rtol"]
    TidalPy.reinit()
    assert TidalPy.config["radial_solver"]["rtol"] == rtol


def test_use_cwd_for_config_merges_a_file_in_the_working_directory(tmp_path, monkeypatch, restore_config):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "TidalPy_Configs.toml").write_text("[radial_solver]\nrtol = 1.0e-7\n")
    TidalPy.config["configs"]["use_cwd_for_config"] = True
    TidalPy.reinit()
    assert TidalPy.config["radial_solver"]["rtol"] == 1.0e-7
