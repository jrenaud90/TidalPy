"""The Structures TOML loader: schema versions, default merging, world and layer validation, and ``load_toml``."""
import re
import warnings

import pytest
import toml

from TidalPy.Structures.configs import toml_loader as tl
from TidalPy.Structures.configs import world_builder


_DELETE = object()


def _valid_terrestrial():
    return {
        "schema_version": "0.2.0",
        "name": "T",
        "type": "terrestrial",
        "radius_m": 6.0e6,
        "mass_kg": 5.0e24,
        "layers": {
            "core": {"layer_index": 0, "radius_outer_m": 3.0e6,
                     "material": {"solid": {"eos": {"model": "constant", "reference_density_kg_m3": 9000.0}}}},
            "mantle": {"layer_index": 1, "radius_outer_m": 6.0e6, "material": "peridotite",
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


def _layer(**keys):
    """A layer table with a unit outer radius and the given keys."""
    config = {"radius_outer_m": 1.0}
    config.update(keys)
    return config


# A full material table: both phases, every law slot, and melting curves.
_FULL_MATERIAL = {
    "latent_heat_j_kg": 4.0e5,
    "solid": {
        "thermal_conductivity_w_mk": 3.3,
        "heat_capacity_j_kgk": 1200.0,
        "eos": {"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 1.3e11,
                "thermal_expansion_1_k": 3.0e-5},
        "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10},
        "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e21},
        "bulk_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e22},
        "shear_rheology": {"model": "andrade"},
    },
    "liquid": {"eos": {"model": "constant", "reference_density_kg_m3": 2800.0, "bulk_modulus_pa": 2.0e10},
               "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.2}},
    "melting": {"solidus": {"model": "constant", "temperature_k": 1600.0},
                "liquidus": {"model": "constant", "temperature_k": 2000.0},
                "weakening": {"model": "henning"}},
}

# A valid value for every scalar layer key.
_SCALAR_VALUES = {
    "mass_kg": 1.0e20,
    "use_tides": False,
    "tidal_scale": 0.5,
    "is_volume_fixed": False,
    "state": "liquid",
    "is_static": False,
    "is_incompressible": True,
    "temperature_k": 1600.0,
    "use_thermal_expansion": True,
    "use_melting": True,
    "use_pressure_melting": True,
    "use_melt_density": True,
    "use_heating": True,
}


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
        _edited(_valid_star, layers={"x": {"radius_outer_m": 1.0}}),
        "have to fill the world",
        id="star-layers-short-of-the-surface"),
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
def test_every_scalar_key_has_a_test_value():
    assert set(_SCALAR_VALUES) == set(tl.LAYER_SCALAR_KEYS)


@pytest.mark.parametrize("config", [
    # A layer needs neither a material (the configuration's [layers] material is taken) nor any model table.
    pytest.param(_layer(), id="no-material"),
    pytest.param(_layer(material="peridotite"), id="material-name"),
    pytest.param(_layer(material={"preset": "peridotite", "latent_heat_j_kg": 3.0e5}), id="material-preset"),
    pytest.param(_layer(material=_FULL_MATERIAL), id="material-table"),
    *[pytest.param({spec_key: spec_value}, id=f"spec-{spec_key}")
      for spec_key, spec_value in (("radius_outer_m", 1.0e6), ("radius_fraction", 0.5), ("volume_fraction", 0.3))],
    *[pytest.param(_layer(state=state, is_static=False, is_incompressible=True), id=f"state-{state}")
      for state in tl.LAYER_STATES],
    *[pytest.param(_layer(**{key: value}), id=f"scalar-{key}") for key, value in sorted(_SCALAR_VALUES.items())],
    pytest.param(_layer(layer_index=3), id="layer-index"),
    pytest.param(
        _layer(
            material=_FULL_MATERIAL,
            shear_rheology={"model": "maxwell"},
            bulk_rheology={"model": "elastic"},
            cooling={"model": "convection"},
            radiogenics={"model": "off"},
        ),
        id="all-models"),
])
def test_valid_layer_passes(config):
    tl.validate_layer_config("L", config)


@pytest.mark.parametrize("config, match", [
    # The inner radius is derived from the previous layer, never supplied.
    pytest.param(_layer(radius_inner_m=0.0), "radius_inner_m", id="radius-inner"),
    pytest.param({"material": "peridotite"}, "exactly one outer-radius", id="no-outer-radius"),
    pytest.param(_layer(radius_fraction=0.5), "multiple outer-radius", id="two-outer-radii"),
    pytest.param(_layer(bogus=1.0), "Unexpected key", id="unknown-key"),
    pytest.param("peridotite", "table of key-value pairs", id="not-a-table"),
    pytest.param(_layer(material=3300.0), "MatPack name or a material table", id="material-not-name-or-table"),
    pytest.param(_layer(state="plasma"), "'state' must be one of", id="unknown-state"),
    pytest.param(_layer(use_tides="yes"), "must be true or false", id="switch-not-a-boolean"),
    pytest.param(_layer(shear_rheology={"alpha": 0.3}), "with a 'model' key", id="model-table-no-model"),
    pytest.param(_layer(cooling="convection"), "with a 'model' key", id="model-table-not-a-table"),
    pytest.param(_layer(magnetics={"model": "dynamo"}), r"unknown table '\[magnetics\]'", id="unknown-model-table"),
    # The material laws live in the material's tables, not on the layer, and the error says where.
    *[pytest.param(_layer(**{table: {"model": "constant"}}), rf"sets '{table}'.*material", id=f"material-{table}")
      for table in ("eos", "shear_viscosity", "bulk_viscosity", "partial_melt")],
])
def test_invalid_layer_raises(config, match):
    with pytest.raises(ValueError, match=match):
        tl.validate_layer_config("L", config)


@pytest.mark.parametrize("key, value", [
    ("class", "solidliquid"),
    ("type", "mantle_rock"),
    ("material_name", "rock"),
    ("is_tidal", True),
    ("is_solid", False),
    ("use_thermal_eos", True),
    ("mean_molecular_weight_kg_mol", 2.3e-3),
    ("adiabatic_index", 1.4),
    ("reference_temperature_k", 1000.0),
])
def test_a_retired_layer_key_names_its_replacement(key, value):
    """A key of the earlier schema is rejected with what replaced it."""
    assert key in tl.RETIRED_LAYER_KEYS
    with pytest.raises(ValueError, match=f"sets '{key}'.*" + re.escape(tl.RETIRED_LAYER_KEYS[key])):
        tl.validate_layer_config("L", _layer(**{key: value}))


@pytest.mark.parametrize("key, value, law, law_key", [
    ("shear_modulus_static_pa", 6.0e10, "shear_modulus", "shear_modulus_pa"),
    ("shear_viscosity_static_pas", 1.0e21, "shear_viscosity", "reference_viscosity_pas"),
    ("bulk_modulus_static_pa", 1.0e11, "eos", "bulk_modulus_pa"),
    ("reference_density_kg_m3", 3000.0, "eos", "reference_density_kg_m3"),
])
def test_material_keys_belong_to_the_material_table(key, value, law, law_key):
    """On the layer they are rejected; a material reads them from its phase's law tables."""
    with pytest.raises(ValueError, match="Unexpected key"):
        tl.validate_layer_config("L", _layer(**{key: value}))
    # The flat placement of the earlier schema is rejected when the material is built, naming the layer's table.
    flat = _layer(material={"solid": {"eos": {"model": "constant"}, key: value}})
    tl.validate_layer_config("L", flat)
    with pytest.raises(ValueError, match=r"\[layers\.L\.material\].*" + key):
        world_builder.construct_layer("L", flat, 0, 0.0, 1.0)
    # In its law table it builds.
    if law == "eos":
        phase = {"eos": {"model": "constant", law_key: value}}
    else:
        phase = {"eos": {"model": "constant"}, law: {"model": "constant", law_key: value}}
    layer = world_builder.construct_layer("L", _layer(material={"solid": phase}), 0, 0.0, 1.0)
    assert layer.get_config_dict()["material"]["solid"][law][law_key] == value


@pytest.mark.parametrize("old_key, new_key, table", [
    ("thermal_conductivity_ref_w_mk", "thermal_conductivity_w_mk", "phase"),
    ("thermal_expansion_ref_1_k", "thermal_expansion_1_k", "eos"),
    ("heat_capacity_ref_j_kgk", "heat_capacity_j_kgk", "phase"),
])
def test_thermal_keys_belong_to_the_phase(old_key, new_key, table):
    """Thermal keys on the layer are rejected; the phase (or its EOS, for the expansivity) holds them."""
    with pytest.raises(ValueError, match="Unexpected key"):
        tl.validate_layer_config("L", _layer(**{old_key: 1.0e-5}))
    phase = {"eos": {"model": "constant"}}
    if table == "eos":
        phase["eos"][new_key] = 1.0e-5
    else:
        phase[new_key] = 1.0e-5
    layer = world_builder.construct_layer("L", _layer(material={"solid": phase}), 0, 0.0, 1.0)
    held = layer.get_config_dict()["material"]["solid"]
    assert (held["eos"] if table == "eos" else held)[new_key] == 1.0e-5


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
