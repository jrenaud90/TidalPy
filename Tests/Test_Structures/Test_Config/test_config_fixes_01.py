"""A layer's material overrides across a model change, the radial data-file reader's checks, and running without a
data directory."""
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
from TidalPy.Material import merge_material_tables
from TidalPy.Structures import build_world, build_world_from_dict, available_worlds
from TidalPy.Structures.configs import data_file, worldpack


@pytest.fixture
def fresh_data_dir_warning(monkeypatch):
    """Forget which data directories were already warned about, so each test sees its own warning."""
    monkeypatch.setattr(paths, "_WARNED_UNUSABLE_DATA_DIRS", set())


def _one_layer_world(layer_cfg):
    layer = {"layer_index": 0, "radius_outer_m": 1.0e6, **layer_cfg}
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
# A preset's model table merges with an override naming its model, and is replaced by one naming another
# =====================================================================================================================
def _preset_layer(preset, **overrides):
    """A one-layer world whose material is a MatPack preset with overrides."""
    return _one_layer_world({"material": {"preset": preset, **overrides}})


def test_an_override_naming_the_presets_model_keeps_its_other_values():
    # "bm" is an alias of the preset's birch_murnaghan law, so the override merges into it.
    world = build_world_from_dict(_preset_layer(
        "peridotite", solid={"eos": {"model": "bm", "reference_density_kg_m3": 3400.0}}))
    eos = world.shell.get_config_dict()["material"]["solid"]["eos"]
    assert eos["model"] == "birch_murnaghan"
    assert eos["reference_density_kg_m3"] == 3400.0
    assert eos["reference_bulk_modulus_pa"] == 1.25e11
    assert eos["bulk_modulus_derivative"] == 4.9


def test_a_users_key_wins_over_the_presets_value():
    world = build_world_from_dict(_preset_layer("ice_ih", melting={"solidus": {"temperature_k": 260.0}}))
    melting = world.shell.get_config_dict()["material"]["melting"]
    assert melting["solidus"]["temperature_k"] == 260.0
    # The rest of the preset's curve and the other curve are kept.
    assert melting["solidus"]["simon_a_pa"] == -4.15e8
    assert melting["liquidus"]["temperature_k"] == 273.16


def test_switching_the_viscosity_model_replaces_the_presets_table():
    table = {"model": "reference", "reference_viscosity_pas": 1.0e14, "reference_temperature_k": 260.0}
    world = build_world_from_dict(_preset_layer("ice_ih", solid={"shear_viscosity": table}))
    viscosity = world.shell.get_config_dict()["material"]["solid"]["shear_viscosity"]
    assert viscosity["reference_viscosity_pas"] == 1.0e14
    # The reference law's own activation energy, not the ice preset's Arrhenius value (5.94e4): another model reads
    # its keys differently, so none carry over.
    assert viscosity["molar_activation_energy_j_mol"] == 3.0e5
    assert "arrhenius_coeff" not in viscosity


def test_a_full_material_table_takes_only_its_own_keys():
    world = build_world_from_dict(_one_layer_world({
        "material": {"solid": {"eos": {"model": "constant", "reference_density_kg_m3": 1000.0}}}}))
    material = world.shell.get_config_dict()["material"]
    # No preset: no liquid phase, melting curves, or laws the table does not name.
    assert set(material) == {"latent_heat_j_kg", "solid"}
    assert "shear_modulus" not in material["solid"] and "shear_viscosity" not in material["solid"]
    assert material["solid"]["eos"]["reference_density_kg_m3"] == 1000.0


def test_a_profile_layers_material_table_merges_over_its_slice():
    config = _two_layer_profile(layer_1={"layer_index": 1, "material": {
        "liquid": {"shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e-3}}}})
    liquid = build_world(config).layer_1.get_config_dict()["material"]["liquid"]
    assert liquid["shear_viscosity"] == {"model": "constant", "reference_viscosity_pas": 1.0e-3}
    # The profile's slice is kept.
    assert liquid["eos"]["model"] == "interpolate"
    assert liquid["eos"]["density_kg_m3"][-1] == 1000.0


def test_a_profile_layer_cannot_name_a_material():
    with pytest.raises(ValueError, match="must be a table of changes"):
        build_world(_two_layer_profile(layer_1={"layer_index": 1, "material": "water"}))


def test_merge_configs_merges_model_tables_key_by_key():
    # No model-change rule: a configuration table keeps the base's keys whichever model the override names.
    base = {"radiogenics": {"model": "isotope", "isotopes": "modern_day_chondritic", "ref_time_s": 1.0e16}}
    merged = merge_configs(base, {"radiogenics": {"model": "fixed"}})
    assert merged == {"radiogenics": {"model": "fixed", "isotopes": "modern_day_chondritic", "ref_time_s": 1.0e16}}
    assert base["radiogenics"]["model"] == "isotope"
    # The packaged [layers] table names the default material and nothing else.
    assert get_packaged_config()["layers"] == {"material": "simple_rock"}


@pytest.mark.parametrize("override_model, keeps_reference_viscosity", [
    ("reference", False),
    # An alias of the base's model is the same model.
    ("const", True),
])
def test_a_material_table_is_replaced_only_on_a_real_model_change(override_model, keeps_reference_viscosity):
    base = {"solid": {"eos": {"model": "constant"},
                      "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e21}}}
    merged = merge_material_tables(base, {"solid": {"shear_viscosity": {"model": override_model}}})
    assert merged["solid"]["shear_viscosity"]["model"] == override_model
    assert ("reference_viscosity_pas" in merged["solid"]["shear_viscosity"]) == keeps_reference_viscosity
    assert merged["solid"]["eos"] == {"model": "constant"}
    # None removes a slot, and the base is untouched.
    assert "shear_viscosity" not in merge_material_tables(base, {"solid": {"shear_viscosity": None}})["solid"]
    assert base["solid"]["shear_viscosity"]["reference_viscosity_pas"] == 1.0e21


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
    return {"material": {"solid": {
        "eos": {"model": "interpolate", "radius_m": radius, "density_kg_m3": [3000.0, 3000.0, 3000.0],
                "bulk_modulus_pa": [1.0e11, 1.0e11, 1.0e11]},
        "shear_modulus": {"model": "interpolate", "radius_m": radius, "shear_modulus_pa": [5.0e10, 5.0e10, 5.0e10]}}}}


def test_an_interpolated_table_that_misses_its_layer_is_refused():
    # A radius given in km under radius_m would otherwise build a uniform layer.
    with pytest.raises(ValueError, match="must cover its layer"):
        build_world_from_dict(_one_layer_world(_interpolated_layer([0.0, 500.0, 1000.0])))


def test_an_interpolated_table_that_covers_its_layer_builds():
    world = build_world_from_dict(_one_layer_world(_interpolated_layer([0.0, 5.0e5, 1.0e6])))
    assert world.shell.get_config_dict()["material"]["solid"]["eos"]["radius_m"][-1] == 1.0e6


def test_a_profile_without_a_repeated_boundary_row_spans_each_layer():
    world = build_world(_two_layer_profile())
    upper = world.get_config_dict()["layers"]["layer_1"]["material"]["liquid"]["eos"]
    # The liquid layer's first row is repeated down at the solid layer's top.
    assert upper["radius_m"][:2] == [1.0e6, 1.1e6]
    assert upper["density_kg_m3"][:2] == [1100.0, 1100.0]
    # The saved configuration rebuilds.
    rebuilt = build_world_from_dict(world.get_config_dict())
    assert rebuilt.get_config_dict()["layers"]["layer_1"]["material"]["liquid"]["eos"]["radius_m"] == \
        upper["radius_m"]


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
