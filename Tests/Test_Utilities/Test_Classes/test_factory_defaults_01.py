"""Tests that a ``make_*`` factory given no config takes the world builder's defaults for the same model."""
import copy

import pytest

import TidalPy
from TidalPy.Cooling.cooling import make_cooling
from TidalPy.Material.eos.material_eos import make_material_eos
from TidalPy.PartialMelt.partial_melt import make_partial_melt
from TidalPy.Radiogenics import radiogenics as radiogenics_module
from TidalPy.Radiogenics.radiogenics import make_radiogenics
from TidalPy.Rheology.rheology import make_rheology
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Utilities.classes import factory_defaults
from TidalPy.Viscosity import make_viscosity


@pytest.fixture()
def config(monkeypatch):
    """A private copy of the configuration that a test may edit."""
    private = copy.deepcopy(TidalPy.config)
    monkeypatch.setattr(TidalPy, "config", private)
    return private


@pytest.mark.parametrize(
    "args, expected",
    [
        pytest.param(("shear_rheology", ("alpha", "zeta")), {"alpha": 0.3, "zeta": 1.0}, id="shear_rheology"),
        pytest.param(("no_such_table", ("alpha",)), {}, id="unknown_table"),
        pytest.param(("material.no_such_table", ("alpha",)), {}, id="unknown_nested_table"),
        pytest.param(("shear_rheology", ("alpha",), "Andrade"), {"alpha": 0.3}, id="named_model_any_case"),
        pytest.param(("shear_rheology", ("alpha",), "maxwell"), {}, id="other_model"),
        # Without the family's resolver an alias does not match; with it, it does.
        pytest.param(("radiogenics", ("isotopes",), "isotopes"), {}, id="alias_without_resolver"),
        pytest.param(
            ("radiogenics", ("isotopes",), "isotopes", radiogenics_module._same_model),
            {"isotopes": "modern_day_chondritic"},
            id="alias_with_resolver"),
        pytest.param(
            ("radiogenics", ("isotopes",), "constant", radiogenics_module._same_model),
            {},
            id="other_model_with_resolver"),
    ])
def test_factory_defaults_reads_the_default_layer_tables(args, expected):
    """factory_defaults returns a table's keys only when the table applies to the named model."""
    assert factory_defaults(*args) == expected


def test_factory_defaults_special_cases():
    """Nested tables are read, the model key is never a parameter, and a table with no model key always applies."""
    viscosity = factory_defaults("material.shear_viscosity", ("reference_viscosity_pas", "reference_temperature_k"))
    assert viscosity["reference_viscosity_pas"] == 1.0e22
    assert "model" not in factory_defaults("shear_rheology", ("model", "alpha"))
    assert factory_defaults("tides", ("fixed_q",), "fixed_q") == factory_defaults("tides", ("fixed_q",))
    # An unknown table name is a mismatch for the family's resolver.
    with pytest.raises(ValueError):
        radiogenics_module._same_model("no_such_model", "isotope")


def test_a_factory_without_config_follows_an_edited_configuration(config):
    """A factory given no config follows an edited default table; an explicit config, even empty, does not."""
    config["layers"]["default"]["shear_rheology"]["alpha"] = 0.21
    config["layers"]["default"]["shear_rheology"]["zeta"] = 7.0
    andrade = make_rheology("andrade")
    assert andrade.alpha == 0.21 and andrade.zeta == 7.0
    assert make_rheology("andrade", {}).alpha == 0.3
    assert make_rheology("andrade", {"alpha": 0.5}).alpha == 0.5


def test_every_family_takes_its_own_table(config):
    """Each factory family reads its own default table."""
    default = config["layers"]["default"]
    default["material"]["shear_viscosity"]["reference_viscosity_pas"] = 4.0e20
    default["material"]["partial_melt"]["solidus_k"] = 1234.0
    default["cooling"]["critical_rayleigh"] = 999.0
    default["material"]["shear_modulus_static_pa"] = 7.7e10
    config["tides"]["fixed_q"] = [33.0]

    assert make_viscosity("reference").reference_viscosity == 4.0e20
    assert make_partial_melt("henning").solidus == 1234.0
    assert make_cooling("convection").critical_rayleigh == 999.0
    material = make_material_eos("constant")
    assert material.shear_modulus_static == 7.7e10
    # The material's nested model tables are attached as the world builder attaches them.
    assert material.shear_viscosity_set and material.partial_melt_set
    assert make_tide("fixed_q").get_config_dict()["fixed_q"][0] == 33.0


def test_a_model_the_table_does_not_name_keeps_its_own_defaults(config):
    """A parameter of the named model never carries over to another model of the family."""
    config["layers"]["default"]["material"]["shear_viscosity"]["reference_viscosity_pas"] = 5.0e21
    # The table names the "reference" law.
    assert make_viscosity("reference").reference_viscosity == 5.0e21
    assert make_viscosity("constant").reference_viscosity == 1.0e22
    # The isotope table's dataset must not set a fixed model's reference time.
    fixed = make_radiogenics("fixed").get_config_dict()
    assert fixed["ref_time_s"] == 0.0


def test_a_model_ignores_the_other_model_keys_of_a_merged_table():
    """Each model takes only its own keys from a table that carries keys for several models."""
    merged = {"isotopes": "modern_day_chondritic", "fixed_heat_production_w_kg": 2.0e-12}
    fixed = make_radiogenics("constant", dict(merged)).get_config_dict()
    assert fixed["ref_time_s"] == 0.0 and fixed["fixed_heat_production_w_kg"] == 2.0e-12
    assert "isotope_names" not in fixed
    isotope = make_radiogenics("isotopes", dict(merged)).get_config_dict()
    assert isotope["isotope_names"] == ["U238", "U235", "Th232", "K40"] and isotope["ref_time_s"] > 1.0e17


def test_the_world_builder_gives_a_fixed_layer_the_fixed_keys_only():
    """A layer naming a fixed radiogenics model over isotope defaults gets no dataset reference time."""
    from TidalPy.Structures import build_world
    config = {
        "name": "fixed_rock", "type": "terrestrial", "radius_m": 2.0e6, "mass_kg": 1.0e23,
        "layers": {"mantle": {
            "class": "solidliquid", "type": "mantle_rock", "radius_fraction": 1.0,
            "radiogenics": {"model": "fixed", "fixed_heat_production_w_kg": 3.0e-12}}}}
    layer_config = build_world(config).get_config_dict()["layers"]["mantle"]["radiogenics"]
    assert layer_config["model"] == "fixed"
    assert layer_config["ref_time_s"] == 0.0 and layer_config["fixed_heat_production_w_kg"] == 3.0e-12


@pytest.mark.parametrize("name", ["isotope", "isotopes", "Isotope"])
def test_the_default_radiogenics_are_the_chondritic_isotopes(name):
    """The default radiogenics table is found by canonical name, alias, or another case."""
    model = make_radiogenics(name)
    assert model.get_config_dict()["isotope_names"] == ["U238", "U235", "Th232", "K40"]
