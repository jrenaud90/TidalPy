"""A key of ``TidalPy_Configs_x.toml`` that nothing reads is reported, once, by its dotted path."""
import warnings

import pytest

import TidalPy
from TidalPy.configurations import (
    find_unknown_config_x_keys, get_packaged_config_x, set_config_x, warn_unknown_config_x_keys)


@pytest.fixture()
def packaged():
    return get_packaged_config_x()


def test_a_file_that_only_uses_known_keys_is_clean(packaged):
    assert find_unknown_config_x_keys(packaged, packaged) == []
    assert find_unknown_config_x_keys({}, packaged) == []


def test_unknown_keys_are_named_by_their_path(packaged):
    overrides = {
        "numerical": {"minimum_viscosty": 1.0, "minimum_viscosity": 2.0},
        "eos_solver": {"rtol": 1.0e-9, "tolerance": 1.0e-9},
        "no_such_section": {"a": 1},
        "warnings": {"stale_worldpack_copy": False, "stale_copy": False},
        "graphics": {"interior": {"gravity_color": "g", "gravity_colour": "g"}, "maps": {}},
    }
    assert find_unknown_config_x_keys(overrides, packaged) == [
        "numerical.minimum_viscosty", "eos_solver.tolerance", "no_such_section", "warnings.stale_copy",
        "graphics.interior.gravity_colour", "graphics.maps"]


def test_a_material_type_of_the_users_own_is_allowed(packaged):
    overrides = {"layers": {
        "my_rock": {"is_tidal": True, "material": {"model": "constant", "reference_density_kg_m3": 3300.0},
                    "shear_rheology": {"model": "andrade", "alpha": 0.2}},
        "iron": {"is_tidl": False, "cooling": {"model": "off"}},
    }}
    # A model table's contents are the factory's business: only the layer-level keys are checked.
    assert find_unknown_config_x_keys(overrides, packaged) == ["layers.iron.is_tidl"]


def test_world_and_tide_tables_follow_the_world_schema(packaged):
    overrides = {
        "worlds": {"albedo": 0.2, "albdo": 0.2, "star": {"luminosity_w": 1.0, "luminsity_w": 1.0},
                   "moon": {"albedo": 0.1}},
        "tides": {"max_degree_l": 3, "max_degree": 3, "default_model": {"star": "cpl", "planet": "cpl"}},
    }
    assert find_unknown_config_x_keys(overrides, packaged) == [
        "worlds.albdo", "worlds.star.luminsity_w", "worlds.moon", "tides.max_degree", "tides.default_model.planet"]


def test_the_warning_names_the_source_and_the_switch_silences_it(packaged):
    overrides = {"numerical": {"minimum_viscosty": 1.0}}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        found = warn_unknown_config_x_keys(overrides, packaged, "The test file")
    assert found == ["numerical.minimum_viscosty"]
    assert len(record) == 1
    message = str(record[0].message)
    assert "The test file" in message and "numerical.minimum_viscosty" in message and "unknown_config_key" in message

    silenced = {"numerical": {"minimum_viscosty": 1.0}, "warnings": {"unknown_config_key": False}}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert warn_unknown_config_x_keys(silenced, packaged, "The test file") == ["numerical.minimum_viscosty"]
    assert not record


def test_set_config_x_checks_an_override(packaged):
    original = TidalPy.config_x
    try:
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            set_config_x({"numerical": {"test_constant": 42.0, "test_constnt": 1.0}})
        assert any("numerical.test_constnt" in str(entry.message) for entry in record)
    finally:
        TidalPy.config_x = original
        from TidalPy.constants import update_constants_x
        update_constants_x()


def test_the_switch_is_a_documented_default():
    assert TidalPy.config_x["warnings"]["unknown_config_key"] is True
