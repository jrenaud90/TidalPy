"""
Tests for the structures_x TOML loader and configuration validation
(``TidalPy.structures_x.configs.toml_loader``).

Covers schema-version compatibility, world/layer validation (required keys,
unknown-key rejection, the layer ``class`` / material ``type`` split, and
model-on-wrong-class rejection), default merging, and loading from a file path vs a
dict.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import warnings

import pytest

from TidalPy.structures_x.configs import toml_loader as tl


# =====================================================================================================================
# Helpers: minimal valid configurations
# =====================================================================================================================
def _valid_terrestrial():
    return {
        "schema_version": "0.2.0",
        "name": "T",
        "type": "terrestrial",
        "radius_m": 6.0e6,
        "mass_kg": 5.0e24,
        "layers": {
            "core": {
                "class": "physics",
                "type": "iron",
                "layer_index": 0,
                "radius_outer_m": 3.0e6,
                "eos": {"model": "constant", "reference_density_kg_m3": 9000.0},
            },
            "mantle": {
                "class": "solidliquid",
                "type": "mantle_rock",
                "layer_index": 1,
                "radius_outer_m": 6.0e6,
                "eos": {"model": "constant", "reference_density_kg_m3": 4000.0},
                "cooling": {"model": "convection"},
            },
        },
    }


def _valid_star():
    return {
        "schema_version": "0.2.0",
        "name": "S",
        "type": "star",
        "radius_m": 7.0e8,
        "mass_kg": 2.0e30,
        "effective_temperature_k": 5772.0,
    }


# =====================================================================================================================
# Schema version
# =====================================================================================================================
def test_schema_version_constant():
    assert tl.SCHEMA_VERSION == "0.2.0"


def test_schema_version_matching_allowed_silently():
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning becomes an error
        assert tl.validate_schema_version({"schema_version": "0.2.0"}) is True


@pytest.mark.parametrize("patch_version", ["0.2.1", "0.2.9", "0.2.123"])
def test_schema_version_patch_difference_allowed_silently(patch_version):
    # Patch (0.0.X) differences are allowed with no warning.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert tl.validate_schema_version({"schema_version": patch_version}) is True


@pytest.mark.parametrize("minor_version", ["0.1.0", "0.3.0", "0.0.5"])
def test_schema_version_minor_difference_warns_but_allowed(minor_version):
    # Minor (0.X.0) differences are allowed but warn.
    with pytest.warns(UserWarning, match="functionality may break"):
        assert tl.validate_schema_version({"schema_version": minor_version}) is True


@pytest.mark.parametrize("major_version", ["1.2.0", "2.0.0", "1.0.0"])
def test_schema_version_major_difference_raises(major_version):
    # Major (X.0.0) differences are not allowed.
    with pytest.raises(ValueError, match="major versions differ"):
        tl.validate_schema_version({"schema_version": major_version})


def test_schema_version_force_bypasses_everything():
    # force=True accepts any version silently, even an incompatible major version.
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning becomes an error
        assert tl.validate_schema_version({"schema_version": "9.9.9"}, force=True) is True


def test_schema_version_missing_warns_but_allowed():
    with pytest.warns(UserWarning):
        assert tl.validate_schema_version({}) is True


# =====================================================================================================================
# Default merging
# =====================================================================================================================
def test_merge_adds_schema_version():
    merged = tl.merge_with_defaults({"name": "X", "type": "star",
                                     "radius_m": 1.0, "mass_kg": 1.0})
    assert merged["schema_version"] == tl.SCHEMA_VERSION


def test_merge_preserves_existing_schema_version():
    merged = tl.merge_with_defaults({"schema_version": "0.2.5"})
    assert merged["schema_version"] == "0.2.5"


def test_merge_does_not_mutate_input():
    original = {"name": "X", "type": "star", "radius_m": 1.0, "mass_kg": 1.0}
    tl.merge_with_defaults(original)
    assert "schema_version" not in original


# =====================================================================================================================
# World validation: success
# =====================================================================================================================
def test_valid_terrestrial_passes():
    tl.validate_world_config(_valid_terrestrial())


def test_valid_star_passes():
    tl.validate_world_config(_valid_star())


# =====================================================================================================================
# World validation: failures
# =====================================================================================================================
def test_missing_type_raises():
    config = _valid_star()
    del config["type"]
    with pytest.raises(ValueError, match="type"):
        tl.validate_world_config(config)


def test_unknown_world_type_raises():
    config = _valid_star()
    config["type"] = "blackhole"
    with pytest.raises(ValueError, match="Unknown world type"):
        tl.validate_world_config(config)


@pytest.mark.parametrize("missing", ["name", "radius_m", "mass_kg"])
def test_missing_required_world_key_raises(missing):
    config = _valid_star()
    del config[missing]
    with pytest.raises(ValueError, match=missing):
        tl.validate_world_config(config)


def test_unknown_world_key_raises():
    config = _valid_star()
    config["surface_gravity"] = 9.8
    with pytest.raises(ValueError, match="Unexpected world-level key"):
        tl.validate_world_config(config)


def test_star_with_layers_raises():
    config = _valid_star()
    config["layers"] = {"x": {"class": "base", "radius_outer_m": 1.0}}
    with pytest.raises(ValueError, match="must not declare any layers"):
        tl.validate_world_config(config)


def test_layered_without_layers_raises():
    config = _valid_terrestrial()
    config["layers"] = {}
    with pytest.raises(ValueError, match="at least one"):
        tl.validate_world_config(config)


def test_tides_table_validated():
    # A well-formed [tides] table passes; unknown keys are rejected (typo protection).
    config = _valid_terrestrial()
    config["tides"] = {"global_tidal_model": "fixed_q", "fixed_q": [100.0], "max_degree_l": 2}
    tl.validate_world_config(config)

    config["tides"] = {"fixed_qq": [100.0]}
    with pytest.raises(ValueError):
        tl.validate_world_config(config)

    config["tides"] = "fixed_q"
    with pytest.raises(ValueError):
        tl.validate_world_config(config)


# =====================================================================================================================
# Layer validation: class / material type
# =====================================================================================================================
def test_layer_missing_class_raises():
    with pytest.raises(ValueError, match="class"):
        tl.validate_layer_config("L", {"radius_outer_m": 1.0})


def test_layer_unknown_class_raises():
    with pytest.raises(ValueError, match="unknown class"):
        tl.validate_layer_config("L", {"class": "plasma", "radius_outer_m": 1.0})


def test_layer_material_type_optional():
    # A layer needs a class but the material type is optional.
    tl.validate_layer_config("L", {"class": "base", "radius_outer_m": 1.0})


@pytest.mark.parametrize("material_type", list(tl.MATERIAL_TYPES))
def test_layer_known_material_types_pass(material_type):
    tl.validate_layer_config("L", {"class": "solidliquid", "type": material_type,
                                   "radius_outer_m": 1.0})


def test_layer_unknown_material_type_raises():
    with pytest.raises(ValueError, match="unknown material type"):
        tl.validate_layer_config("L", {"class": "solidliquid", "type": "cheese",
                                       "radius_outer_m": 1.0})


# =====================================================================================================================
# Layer validation: outer-radius specifier (inner radius is derived, never supplied)
# =====================================================================================================================
def test_layer_rejects_radius_inner_m():
    # The inner radius is derived from the previous layer; the user must not supply it.
    with pytest.raises(ValueError, match="radius_inner_m"):
        tl.validate_layer_config("L", {"class": "base", "radius_inner_m": 0.0,
                                       "radius_outer_m": 1.0})


def test_layer_missing_outer_radius_spec_raises():
    with pytest.raises(ValueError, match="exactly one outer-radius"):
        tl.validate_layer_config("L", {"class": "base"})


@pytest.mark.parametrize("spec_key, spec_value", [
    ("radius_outer_m", 1.0e6),
    ("radius_fraction", 0.5),
    ("volume_fraction", 0.3),
])
def test_layer_each_outer_radius_spec_passes(spec_key, spec_value):
    tl.validate_layer_config("L", {"class": "base", spec_key: spec_value})


def test_layer_multiple_outer_radius_specs_raises():
    with pytest.raises(ValueError, match="multiple outer-radius"):
        tl.validate_layer_config("L", {"class": "base", "radius_outer_m": 1.0,
                                       "radius_fraction": 0.5})


def test_layer_unknown_scalar_key_raises():
    with pytest.raises(ValueError, match="Unexpected key"):
        tl.validate_layer_config("L", {"class": "base", "radius_outer_m": 1.0, "bogus": 1.0})


def test_base_layer_rejects_rheology():
    # A geometry-only base layer cannot hold a rheology model.
    cfg = {"class": "base", "radius_outer_m": 1.0,
           "shear_rheology": {"model": "maxwell"}}
    with pytest.raises(ValueError, match="cannot hold"):
        tl.validate_layer_config("L", cfg)


def test_physics_layer_rejects_cooling():
    # cooling/radiogenics are solidliquid-only.
    cfg = {"class": "physics", "radius_outer_m": 1.0,
           "cooling": {"model": "convection"}}
    with pytest.raises(ValueError, match="cannot hold"):
        tl.validate_layer_config("L", cfg)


def test_solidliquid_layer_accepts_all_models():
    cfg = {
        "class": "solidliquid", "radius_outer_m": 1.0,
        "eos": {"model": "constant"},
        "shear_rheology": {"model": "maxwell"},
        "bulk_rheology": {"model": "elastic"},
        "shear_viscosity": {"model": "constant"},
        "bulk_viscosity": {"model": "constant"},
        "partial_melt": {"model": "off"},
        "cooling": {"model": "convection"},
        "radiogenics": {"model": "off"},
    }
    tl.validate_layer_config("L", cfg)


def test_model_table_missing_model_key_raises():
    cfg = {"class": "physics", "radius_outer_m": 1.0,
           "shear_rheology": {"alpha": 0.3}}
    with pytest.raises(ValueError, match="missing the required 'model'"):
        tl.validate_layer_config("L", cfg)


def test_unknown_model_table_raises():
    cfg = {"class": "physics", "radius_outer_m": 1.0,
           "magnetics": {"model": "dynamo"}}
    with pytest.raises(ValueError, match="unknown model table"):
        tl.validate_layer_config("L", cfg)


# =====================================================================================================================
# load_toml
# =====================================================================================================================
def test_load_toml_from_dict_returns_copy():
    source = _valid_star()
    loaded = tl.load_toml(source)
    assert loaded == source
    assert loaded is not source


def test_load_toml_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        tl.load_toml("definitely_not_a_real_file_12345.toml")


def test_load_toml_bad_type_raises():
    with pytest.raises(TypeError):
        tl.load_toml(12345)


def test_load_toml_from_file(tmp_path):
    import toml
    path = tmp_path / "world.toml"
    with open(path, "w") as handle:
        toml.dump(_valid_star(), handle)
    loaded = tl.load_toml(str(path))
    assert loaded["name"] == "S"
    assert loaded["type"] == "star"
