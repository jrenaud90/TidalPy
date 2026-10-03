"""Tests that a ``TidalPy_Configs.toml`` key nothing reads is reported, once, by its dotted path."""
import warnings

import pytest

import TidalPy
from TidalPy.configurations import (
    RETIRED_LAYER_BLOCKS_REASON, drop_retired_config_keys, find_unknown_config_keys, get_packaged_config, set_config,
    warn_unknown_config_keys)


@pytest.fixture()
def packaged():
    return get_packaged_config()


@pytest.mark.parametrize(
    "overrides, expected",
    [
        pytest.param(None, [], id="packaged_config_is_clean"),
        pytest.param({}, [], id="empty_is_clean"),
        pytest.param(
            {
                "numerical": {"minimum_modulos": 1.0, "minimum_modulus": 2.0},
                "eos_solver": {"rtol": 1.0e-9, "tolerance": 1.0e-9},
                "no_such_section": {"a": 1},
                "warnings": {"stale_worldpack_copy": False, "stale_copy": False},
                "graphics": {"interior": {"gravity_color": "g", "gravity_colour": "g"}, "maps": {}},
            },
            ["numerical.minimum_modulos", "eos_solver.tolerance", "no_such_section", "warnings.stale_copy",
             "graphics.interior.gravity_colour", "graphics.maps"],
            id="named_by_path"),
        # [layers] holds the default material and nothing else.
        pytest.param(
            {"layers": {"material": "peridotite", "matrial": "ice_ih"}},
            ["layers.matrial"],
            id="layers_default_material"),
        pytest.param(
            {
                "worlds": {"albedo": 0.2, "albdo": 0.2, "star": {"luminosity_w": 1.0, "luminsity_w": 1.0},
                           "moon": {"albedo": 0.1}},
                "tides": {"max_degree_l": 3, "max_degree": 3, "default_model": {"star": "cpl", "planet": "cpl"}},
            },
            ["worlds.albdo", "worlds.star.luminsity_w", "worlds.moon", "tides.max_degree",
             "tides.default_model.planet"],
            id="world_and_tide_tables_follow_world_schema"),
        # The datasets are named by the user; only the section's own keys are checked.
        pytest.param(
            {"radiogenics": {
                "known_isotope_data": {"my_dataset": {"ref_time": 0.0, "U238": {"hpr": 9.5e-5, "half_life": 4470.0}}},
                "known_isotopes": {}}},
            ["radiogenics.known_isotopes"],
            id="user_isotope_datasets_allowed"),
    ])
def test_find_unknown_config_keys(packaged, overrides, expected):
    """Unknown keys are listed by dotted path; known and user-named keys are not."""
    if overrides is None:
        overrides = packaged
    assert find_unknown_config_keys(overrides, packaged) == expected


def test_the_warning_names_the_source_and_the_switch_silences_it(packaged):
    """The warning names the source file and key, and the warnings switch silences it."""
    overrides = {"numerical": {"minimum_modulos": 1.0}}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        found = warn_unknown_config_keys(overrides, packaged, "The test file")
    assert found == ["numerical.minimum_modulos"]
    assert len(record) == 1
    message = str(record[0].message)
    assert "The test file" in message and "numerical.minimum_modulos" in message and "unknown_config_key" in message

    silenced = {"numerical": {"minimum_modulos": 1.0}, "warnings": {"unknown_config_key": False}}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert warn_unknown_config_keys(silenced, packaged, "The test file") == ["numerical.minimum_modulos"]
    assert not record


def test_set_config_checks_an_override(packaged):
    """set_config warns about an unknown key in an override."""
    original = TidalPy.config
    try:
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            set_config({"numerical": {"test_constant": 42.0, "test_constnt": 1.0}})
        assert any("numerical.test_constnt" in str(entry.message) for entry in record)
    finally:
        TidalPy.config = original
        from TidalPy.constants import update_constants
        update_constants()


def test_retired_layer_blocks_are_dropped_with_one_warning():
    """The per-material [layers.<type>] tables of earlier builds are removed, with one warning naming each."""
    overrides = {"layers": {"material": "peridotite", "ice": {"cooling": {"model": "off"}},
                            "mantle_rock": {"shear_rheology": {"model": "andrade"}}},
                 "numerical": {"minimum_modulus": 2.0}}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        removed = drop_retired_config_keys(overrides, "The test file")
    assert removed == ["layers.ice", "layers.mantle_rock"]
    assert overrides == {"layers": {"material": "peridotite"}, "numerical": {"minimum_modulus": 2.0}}
    assert len(record) == 1
    message = str(record[0].message)
    assert "The test file" in message and "layers.ice" in message and "layers.mantle_rock" in message
    assert RETIRED_LAYER_BLOCKS_REASON in message

    # Nothing to drop is silent, and the switch silences a drop.
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert drop_retired_config_keys({"layers": {"material": "peridotite"}}, "The test file") == []
        silenced = {"layers": {"iron": {}}, "warnings": {"unknown_config_key": False}}
        assert drop_retired_config_keys(silenced, "The test file") == ["layers.iron"]
    assert not record


def test_a_material_table_under_layers_is_kept():
    """[layers] material may be a material table; it is the default, not a retired block."""
    material = {"preset": "simple_rock", "solid": {"eos": {"reference_density_kg_m3": 3.0e3}}}
    overrides = {"layers": {"material": material}}
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert drop_retired_config_keys(overrides, "The test file") == []
    assert not record
    assert overrides["layers"]["material"]["preset"] == "simple_rock"


def test_set_config_drops_a_retired_layer_block(packaged):
    """set_config drops a [layers.<type>] table, which then is not also reported as unknown."""
    original = TidalPy.config
    try:
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            set_config({"layers": {"material": "simple_ice", "ice": {"cooling": {"model": "off"}}}})
        messages = [str(entry.message) for entry in record]
        assert sum("layers.ice" in message for message in messages) == 1
        assert not any("does not read" in message for message in messages)
        assert TidalPy.config["layers"] == {"material": "simple_ice"}
    finally:
        TidalPy.config = original
        from TidalPy.constants import update_constants
        update_constants()


def test_the_switch_is_a_documented_default():
    """The unknown key warning is on by default."""
    assert TidalPy.config["warnings"]["unknown_config_key"] is True
