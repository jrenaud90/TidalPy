"""Tests that a ``make_*`` factory given no config takes its model's own defaults, except the two families with a
configuration table: tide models read ``[tides]`` and radiogenics reads ``[radiogenics]`` (its isotope dataset)."""
import copy

import pytest

import TidalPy
from TidalPy.Cooling.cooling import make_cooling
from TidalPy.Radiogenics.radiogenics import RADIOGENICS_CONFIG_KEYS, make_radiogenics
from TidalPy.Rheology.rheology import _FAMILY as RHEOLOGY_FAMILY
from TidalPy.Rheology.rheology import make_rheology
from TidalPy.Tides.classes.tide import TIDE_CONFIG_KEYS, make_tide
from TidalPy.Utilities.classes import factory_defaults
from TidalPy.Viscosity import make_viscosity
from TidalPy.Viscosity.viscosity import _FAMILY as VISCOSITY_FAMILY


@pytest.fixture()
def config(monkeypatch):
    """A private copy of the configuration that a test may edit."""
    private = copy.deepcopy(TidalPy.config)
    monkeypatch.setattr(TidalPy, "config", private)
    return private


def spec_defaults(model) -> dict:
    """A model's parameters as its parameter spec states their defaults, by config key."""
    return {entry["key"]: entry["default"] for entry in model.get_parameter_info()}


@pytest.mark.parametrize(
    "args, expected",
    [
        pytest.param(("radiogenics", RADIOGENICS_CONFIG_KEYS), {"isotopes": "modern_day_chondritic"},
                     id="radiogenics"),
        pytest.param(("radiogenics", ("ref_time_s",)), {}, id="only_accepted_keys"),
        pytest.param(("tides", ("fixed_q",)), {"fixed_q": [100.0] * 9}, id="tides"),
        pytest.param(("no_such_table", ("alpha",)), {}, id="unknown_table"),
        # A top-level scalar is not a table.
        pytest.param(("schema_version", ("alpha",)), {}, id="scalar_section"),
    ])
def test_factory_defaults_reads_a_top_level_table(args, expected):
    """factory_defaults returns the accepted keys of a top-level table, and nothing for a missing or scalar one."""
    assert factory_defaults(*args) == expected


def test_factory_defaults_special_cases(config):
    """The model key is never a parameter, and the tide table's per-family and default-model tables are left out."""
    config["radiogenics"]["model"] = "fixed"
    assert "model" not in factory_defaults("radiogenics", ("model", "isotopes"))
    tides = factory_defaults("tides", TIDE_CONFIG_KEYS)
    assert set(tides) == {"fixed_k", "fixed_q", "fixed_dt_s"}
    # The per-layer families have no top-level table.
    for section in ("shear_rheology", "viscosity", "cooling"):
        assert factory_defaults(section, ("alpha", "reference_viscosity_pas", "critical_rayleigh")) == {}


@pytest.mark.parametrize("model_name", RHEOLOGY_FAMILY.model_names())
def test_a_rheology_without_config_takes_its_spec_defaults(model_name):
    """make_rheology(name) builds the model its parameter spec describes, the same as an empty config."""
    model = make_rheology(model_name)
    assert {key: value for key, value in model.get_config_dict().items() if key != "model"} == spec_defaults(model)
    assert model.get_config_dict() == make_rheology(model_name, {}).get_config_dict()


@pytest.mark.parametrize("model_name", ["arrhenius", "reference", "constant", "interpolate"])
def test_a_viscosity_without_config_takes_its_spec_defaults(model_name):
    """make_viscosity(name) builds the model its parameter spec describes, the same as an empty config."""
    assert model_name in VISCOSITY_FAMILY.model_names()
    model = make_viscosity(model_name)
    assert model.parameters == make_viscosity(model_name, {}).parameters
    for key, default in spec_defaults(model).items():
        assert model.get_parameter(key) == default


@pytest.mark.parametrize("model_name", ["off", "convection", "conduction"])
def test_a_cooling_model_without_config_takes_its_own_defaults(model_name):
    """make_cooling(name) is the same as an empty config: the C++ model's defaults."""
    assert make_cooling(model_name).get_config_dict() == make_cooling(model_name, {}).get_config_dict()


def test_the_convection_defaults():
    """The convection model's own defaults."""
    convection = make_cooling("convection").get_config_dict()
    assert convection["convection_alpha"] == 1.0
    assert convection["convection_beta"] == pytest.approx(1.0 / 3.0, rel=1.0e-15)
    assert convection["critical_rayleigh"] == 1100.0


def test_a_spec_family_ignores_the_configuration(config):
    """Tables a configuration might carry for the per-layer families do not reach their factories."""
    config["shear_rheology"] = {"model": "andrade", "alpha": 0.21, "zeta": 7.0}
    config["viscosity"] = {"model": "reference", "reference_viscosity_pas": 4.0e20}
    config["cooling"] = {"model": "convection", "critical_rayleigh": 999.0}
    config["layers"]["shear_rheology"] = {"model": "andrade", "alpha": 0.21}
    andrade = make_rheology("andrade")
    assert andrade.alpha == 0.3 and andrade.zeta == 1.0
    assert make_rheology("andrade", {"alpha": 0.5}).alpha == 0.5
    assert make_viscosity("reference").reference_viscosity == 1.0e22
    assert make_cooling("convection").get_config_dict()["critical_rayleigh"] == 1100.0


def test_the_tide_and_radiogenics_factories_follow_an_edited_configuration(config):
    """make_tide follows an edited [tides] table, make_radiogenics an edited [radiogenics] isotopes dataset."""
    config["tides"]["fixed_q"] = [33.0]
    assert make_tide("fixed_q").get_config_dict()["fixed_q"][0] == 33.0
    # A given config wins over the table for the keys it holds.
    assert make_tide("fixed_q", {"fixed_q": [44.0]}).get_config_dict()["fixed_q"][0] == 44.0

    chondritic = make_radiogenics("isotope").get_config_dict()
    config["radiogenics"]["isotopes"] = "llri_and_slri"
    short_lived = make_radiogenics("isotope").get_config_dict()
    assert short_lived["isotope_names"] == ["U238", "U235", "Th232", "K40", "Al26", "Fe60", "Mn53"]
    assert short_lived["isotope_names"] != chondritic["isotope_names"]
    # A config naming a dataset or giving arrays is used as given; one naming neither, whatever else it holds, takes
    # the configured dataset, so an isotope model is never left without isotopes by accident.
    assert make_radiogenics("isotope", {"isotopes": "modern_day_chondritic"}).get_config_dict() == chondritic
    assert make_radiogenics("isotope", {}).get_config_dict() == short_lived
    given_time = make_radiogenics("isotope", {"ref_time_s": 1.0e16}).get_config_dict()
    assert given_time["isotope_names"] == short_lived["isotope_names"] and given_time["ref_time_s"] == 1.0e16
    no_isotopes = {key: [] for key in ("heat_production_w_kg", "half_lives_s", "mass_fracs", "concentrations",
                                       "isotope_names")}
    assert make_radiogenics("isotope", no_isotopes).get_config_dict()["isotope_names"] == []


def test_the_isotope_dataset_does_not_reach_a_fixed_model():
    """The [radiogenics] isotope dataset does not set a fixed model's reference time."""
    fixed = make_radiogenics("fixed").get_config_dict()
    assert fixed["ref_time_s"] == 0.0
    assert fixed["fixed_heat_production_w_kg"] == 0.0


def test_a_model_ignores_the_other_model_keys_of_a_merged_table():
    """Each model takes only its own keys from a table that carries keys for several models."""
    merged = {"isotopes": "modern_day_chondritic", "fixed_heat_production_w_kg": 2.0e-12}
    fixed = make_radiogenics("constant", dict(merged)).get_config_dict()
    assert fixed["ref_time_s"] == 0.0 and fixed["fixed_heat_production_w_kg"] == 2.0e-12
    assert "isotope_names" not in fixed
    isotope = make_radiogenics("isotopes", dict(merged)).get_config_dict()
    assert isotope["isotope_names"] == ["U238", "U235", "Th232", "K40"] and isotope["ref_time_s"] > 1.0e17


def test_the_world_builder_gives_a_fixed_layer_the_fixed_keys_only():
    """A layer naming a fixed radiogenics model gets no dataset reference time."""
    from TidalPy.Structures import build_world
    config = {
        "name": "fixed_rock", "type": "terrestrial", "radius_m": 2.0e6, "mass_kg": 1.0e23,
        "layers": {"mantle": {
            "radius_fraction": 1.0,
            "radiogenics": {"model": "fixed", "fixed_heat_production_w_kg": 3.0e-12}}}}
    layer_config = build_world(config).get_config_dict()["layers"]["mantle"]["radiogenics"]
    assert layer_config["model"] == "fixed"
    assert layer_config["ref_time_s"] == 0.0 and layer_config["fixed_heat_production_w_kg"] == 3.0e-12


@pytest.mark.parametrize("name", ["isotope", "isotopes", "Isotope"])
def test_the_default_radiogenics_are_the_chondritic_isotopes(name):
    """The [radiogenics] dataset is taken by canonical name, alias, or another case."""
    model = make_radiogenics(name)
    assert model.get_config_dict()["isotope_names"] == ["U238", "U235", "Th232", "K40"]
