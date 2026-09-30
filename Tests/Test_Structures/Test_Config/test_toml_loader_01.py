"""The Structures TOML loader: schema versions, default merging, world and layer validation, and ``load_toml``."""
import warnings

import pytest
import toml

from TidalPy.Structures.configs import toml_loader as tl


_DELETE = object()


def _valid_terrestrial():
    return {
        "schema_version": "0.2.0",
        "name": "T",
        "type": "terrestrial",
        "radius_m": 6.0e6,
        "mass_kg": 5.0e24,
        "layers": {
            "core": {"class": "base", "type": "iron", "layer_index": 0, "radius_outer_m": 3.0e6,
                     "material": {"model": "constant", "reference_density_kg_m3": 9000.0}},
            "mantle": {"class": "solidliquid", "type": "mantle_rock", "layer_index": 1, "radius_outer_m": 6.0e6,
                       "material": {"model": "constant", "reference_density_kg_m3": 4000.0},
                       "cooling": {"model": "convection"}},
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


def _edited(make_config, **changes):
    """A valid config with top-level keys replaced; ``_DELETE`` removes one."""
    config = make_config()
    for key, value in changes.items():
        if value is _DELETE:
            del config[key]
        else:
            config[key] = value
    return config


def _layer(layer_class="base", **keys):
    """A layer table with a unit outer radius; ``layer_class=None`` leaves out the class."""
    config = {} if layer_class is None else {"class": layer_class}
    config["radius_outer_m"] = 1.0
    config.update(keys)
    return config


# =====================================================================================================================
# Schema version
# =====================================================================================================================
def test_schema_version_constant():
    assert tl.SCHEMA_VERSION == "0.2.0"


@pytest.mark.parametrize("version, force", [
    ("0.2.0", False),
    ("0.2.1", False),
    ("0.2.9", False),
    ("0.2.123", False),
    # force accepts even an incompatible major version.
    ("9.9.9", True),
])
def test_schema_version_allowed_silently(version, force):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert tl.validate_schema_version({"schema_version": version}, force=force) is True


@pytest.mark.parametrize("config, match", [
    pytest.param({"schema_version": "0.1.0"}, "functionality may break", id="0.1.0"),
    pytest.param({"schema_version": "0.3.0"}, "functionality may break", id="0.3.0"),
    pytest.param({"schema_version": "0.0.5"}, "functionality may break", id="0.0.5"),
    pytest.param({}, None, id="missing"),
])
def test_schema_version_warns_but_allowed(config, match):
    with pytest.warns(UserWarning, match=match):
        assert tl.validate_schema_version(config) is True


@pytest.mark.parametrize("major_version", ["1.2.0", "2.0.0", "1.0.0"])
def test_schema_version_major_difference_raises(major_version):
    with pytest.raises(ValueError, match="major versions differ"):
        tl.validate_schema_version({"schema_version": major_version})


# =====================================================================================================================
# Default merging
# =====================================================================================================================
@pytest.mark.parametrize("config, expected", [
    pytest.param({"name": "X", "type": "star", "radius_m": 1.0, "mass_kg": 1.0}, tl.SCHEMA_VERSION, id="added"),
    pytest.param({"schema_version": "0.2.5"}, "0.2.5", id="preserved"),
])
def test_merge_sets_the_schema_version(config, expected):
    assert tl.merge_with_defaults(config)["schema_version"] == expected


def test_merge_does_not_mutate_input():
    original = {"name": "X", "type": "star", "radius_m": 1.0, "mass_kg": 1.0}
    tl.merge_with_defaults(original)
    assert "schema_version" not in original


# =====================================================================================================================
# World validation
# =====================================================================================================================
@pytest.mark.parametrize("config", [
    pytest.param(_valid_terrestrial(), id="terrestrial"),
    pytest.param(_valid_star(), id="star"),
    pytest.param(
        _edited(_valid_terrestrial, tides={"global_tidal_model": "fixed_q", "fixed_q": [100.0], "max_degree_l": 2}),
        id="tides-table"),
])
def test_valid_world_passes(config):
    tl.validate_world_config(config)


@pytest.mark.parametrize("config, match", [
    pytest.param(_edited(_valid_star, type=_DELETE), "type", id="no-type"),
    pytest.param(_edited(_valid_star, type="blackhole"), "Unknown world type", id="unknown-type"),
    *[pytest.param(_edited(_valid_star, **{key: _DELETE}), key, id=f"no-{key}")
      for key in ("name", "radius_m", "mass_kg")],
    pytest.param(_edited(_valid_star, surface_gravity=9.8), "Unexpected world-level key", id="unknown-key"),
    pytest.param(
        _edited(_valid_star, layers={"x": {"class": "base", "radius_outer_m": 1.0}}),
        "must not declare any layers",
        id="star-with-layers"),
    pytest.param(_edited(_valid_terrestrial, layers={}), "at least one", id="no-layers"),
    # Unknown [tides] keys are rejected as typo protection.
    pytest.param(_edited(_valid_terrestrial, tides={"fixed_qq": [100.0]}), None, id="tides-unknown-key"),
    pytest.param(_edited(_valid_terrestrial, tides="fixed_q"), None, id="tides-not-a-table"),
])
def test_invalid_world_raises(config, match):
    with pytest.raises(ValueError, match=match):
        tl.validate_world_config(config)


# =====================================================================================================================
# Layer validation
# =====================================================================================================================
@pytest.mark.parametrize("config", [
    # A layer needs a class, but the material type is optional.
    pytest.param(_layer(), id="no-material-type"),
    *[pytest.param(_layer("solidliquid", type=material_type), id=f"type-{material_type}")
      for material_type in tl.MATERIAL_TYPES],
    *[pytest.param({"class": "base", spec_key: spec_value}, id=f"spec-{spec_key}")
      for spec_key, spec_value in (("radius_outer_m", 1.0e6), ("radius_fraction", 0.5), ("volume_fraction", 0.3))],
    *[pytest.param(_layer(layer_class, is_solid=False, is_static=False, is_incompressible=True),
                   id=f"flags-{layer_class}")
      for layer_class in ("base", "solidliquid", "gas")],
    pytest.param(
        _layer(
            "solidliquid",
            material={"model": "constant", "shear_viscosity": {"model": "constant"},
                      "bulk_viscosity": {"model": "constant"}, "partial_melt": {"model": "off"}},
            shear_rheology={"model": "maxwell"},
            bulk_rheology={"model": "elastic"},
            cooling={"model": "convection"},
            radiogenics={"model": "off"},
        ),
        id="solidliquid-all-models"),
])
def test_valid_layer_passes(config):
    tl.validate_layer_config("L", config)


@pytest.mark.parametrize("config, match", [
    pytest.param(_layer(None), "class", id="no-class"),
    pytest.param(_layer("plasma"), "unknown class", id="unknown-class"),
    pytest.param(_layer("solidliquid", type="cheese"), "unknown material type", id="unknown-material-type"),
    # The inner radius is derived from the previous layer, never supplied.
    pytest.param(_layer(radius_inner_m=0.0), "radius_inner_m", id="radius-inner"),
    pytest.param({"class": "base"}, "exactly one outer-radius", id="no-outer-radius"),
    pytest.param(_layer(radius_fraction=0.5), "multiple outer-radius", id="two-outer-radii"),
    pytest.param(_layer(bogus=1.0), "Unexpected key", id="unknown-key"),
    # Cooling and radiogenics are solidliquid-only.
    pytest.param(_layer("base", cooling={"model": "convection"}), "cannot hold", id="base-cooling"),
    pytest.param(
        _layer("base", shear_rheology={"alpha": 0.3}), "missing the required 'model'", id="model-table-no-model"),
    pytest.param(_layer("base", magnetics={"model": "dynamo"}), "unknown model table", id="unknown-model-table"),
    *[pytest.param(_layer("base", **{table: {"model": "constant"}}), r"layers\.L\.material", id=f"moved-{table}")
      for table in ("eos", "shear_viscosity", "bulk_viscosity", "partial_melt")],
])
def test_invalid_layer_raises(config, match):
    with pytest.raises(ValueError, match=match):
        tl.validate_layer_config("L", config)


@pytest.mark.parametrize("layer_class", ["base", "solidliquid", "gas"])
@pytest.mark.parametrize("key, value", [("temperature_k", 1600.0), ("use_thermal_eos", True), ("use_heating", True)])
def test_layer_state_keys_are_schema_keys(key, value, layer_class):
    tl.validate_layer_config("L", _layer(layer_class, **{key: value}))


@pytest.mark.parametrize("layer_class", ["base", "solidliquid", "gas"])
@pytest.mark.parametrize("key, value", [
    ("shear_modulus_static_pa", 6.0e10),
    ("shear_viscosity_static_pas", 1.0e21),
    ("shear_modulus_pressure_derivative", 1.4),
    ("shear_modulus_temperature_derivative_pa_k", -8.0e6),
    ("shear_modulus_reference_temperature_k", 1600.0),
])
def test_material_keys_belong_to_the_material_table(key, value, layer_class):
    """On the layer they are rejected with a pointer to the material; in the material table they pass."""
    with pytest.raises(ValueError, match=r"layers\.L\.material"):
        tl.validate_layer_config("L", _layer(layer_class, **{key: value}))
    tl.validate_layer_config("L", _layer(layer_class, material={key: value}))


@pytest.mark.parametrize("old_key, new_key", [
    ("thermal_conductivity_ref_w_mk", "thermal_conductivity_w_mk"),
    ("thermal_expansion_ref_1_k", "thermal_expansion_1_k"),
    ("heat_capacity_ref_j_kgk", "heat_capacity_j_kgk"),
])
def test_thermal_keys_left_on_the_layer_name_their_material_key(old_key, new_key):
    with pytest.raises(ValueError, match=new_key + r".*layers\.L\.material"):
        tl.validate_layer_config("L", _layer("solidliquid", **{old_key: 1.0}))
    tl.validate_layer_config("L", _layer("solidliquid", material={new_key: 1.0}))


# =====================================================================================================================
# load_toml
# =====================================================================================================================
def test_load_toml_from_dict_returns_copy():
    source = _valid_star()
    loaded = tl.load_toml(source)
    assert loaded == source
    assert loaded is not source


def test_load_toml_from_file(tmp_path):
    path = tmp_path / "world.toml"
    with open(path, "w") as handle:
        toml.dump(_valid_star(), handle)
    loaded = tl.load_toml(str(path))
    assert loaded["name"] == "S"
    assert loaded["type"] == "star"


@pytest.mark.parametrize("source, error", [
    pytest.param("definitely_not_a_real_file_12345.toml", FileNotFoundError, id="missing-file"),
    pytest.param(12345, TypeError, id="bad-type"),
])
def test_load_toml_bad_source_raises(source, error):
    with pytest.raises(error):
        tl.load_toml(source)
