"""A ``make_*`` factory called with no config takes the defaults the world builder would give the same model.

The second tier of the world builder's defaults (``[layers.default.<section>]`` of ``TidalPy_Configs_x.toml``, and
``[tides]`` for a tide model) also serves a model built directly through its factory, so one configuration file
describes both paths. An explicit config, even an empty one, is used as given.
"""
import copy

import pytest

import TidalPy
from TidalPy.cooling_x.cooling import make_cooling
from TidalPy.Material_x.eos.material_eos import make_material_eos
from TidalPy.partial_melt_x.partial_melt import make_partial_melt
from TidalPy.radiogenics_x.radiogenics import make_radiogenics
from TidalPy.rheology_x.rheology import make_rheology
from TidalPy.Tides_x.classes.tide import make_tide
from TidalPy.Utilities_x.classes_x import factory_defaults
from TidalPy.viscosity_x import make_viscosity


@pytest.fixture()
def config_x(monkeypatch):
    """A private copy of the configuration that a test may edit."""
    private = copy.deepcopy(TidalPy.config_x)
    monkeypatch.setattr(TidalPy, "config_x", private)
    return private


def test_factory_defaults_reads_the_default_layer_tables():
    rheology = factory_defaults("shear_rheology", ("alpha", "zeta"))
    assert rheology == {"alpha": 0.3, "zeta": 1.0}
    viscosity = factory_defaults("material.shear_viscosity", ("reference_viscosity_pas", "reference_temperature_k"))
    assert viscosity["reference_viscosity_pas"] == 1.0e22
    # The table's model name is never a parameter, and an unknown table is empty.
    assert "model" not in factory_defaults("shear_rheology", ("model", "alpha"))
    assert factory_defaults("no_such_table", ("alpha",)) == {}
    assert factory_defaults("material.no_such_table", ("alpha",)) == {}
    # Named, the model has to be the one the table names (in any case); a table with no model key always applies.
    assert factory_defaults("shear_rheology", ("alpha",), "Andrade") == {"alpha": 0.3}
    assert factory_defaults("shear_rheology", ("alpha",), "maxwell") == {}
    assert factory_defaults("tides", ("fixed_q",), "fixed_q") == factory_defaults("tides", ("fixed_q",))
    # A family's own resolver makes an alias match, and an unknown table name a mismatch.
    from TidalPy.radiogenics_x import radiogenics as radiogenics_module
    same = radiogenics_module._same_model
    assert factory_defaults("radiogenics", ("isotopes",), "isotopes") == {}
    assert factory_defaults("radiogenics", ("isotopes",), "isotopes", same) == {"isotopes": "modern_day_chondritic"}
    assert factory_defaults("radiogenics", ("isotopes",), "constant", same) == {}
    with pytest.raises(ValueError):
        same("no_such_model", "isotope")


def test_a_factory_without_config_follows_an_edited_configuration(config_x):
    config_x["layers"]["default"]["shear_rheology"]["alpha"] = 0.21
    config_x["layers"]["default"]["shear_rheology"]["zeta"] = 7.0
    andrade = make_rheology("andrade")
    assert andrade.alpha == 0.21 and andrade.zeta == 7.0
    # An explicit config, even an empty one, asks for the model's own defaults instead.
    assert make_rheology("andrade", {}).alpha == 0.3
    assert make_rheology("andrade", {"alpha": 0.5}).alpha == 0.5


def test_every_family_takes_its_own_table(config_x):
    default = config_x["layers"]["default"]
    default["material"]["shear_viscosity"]["reference_viscosity_pas"] = 4.0e20
    default["material"]["partial_melt"]["solidus_k"] = 1234.0
    default["cooling"]["critical_rayleigh"] = 999.0
    default["material"]["shear_modulus_static_pa"] = 7.7e10
    config_x["tides"]["fixed_q"] = [33.0]

    assert make_viscosity("reference").reference_viscosity == 4.0e20
    assert make_partial_melt("henning").solidus == 1234.0
    assert make_cooling("convection").critical_rayleigh == 999.0
    material = make_material_eos("constant")
    assert material.shear_modulus_static == 7.7e10
    # The nested model tables of the material come along, attached exactly as the world builder attaches them.
    assert material.shear_viscosity_set and material.partial_melt_set
    assert make_tide("fixed_q").get_config_dict()["fixed_q"][0] == 33.0


def test_a_model_the_table_does_not_name_keeps_its_own_defaults(config_x):
    """A parameter of the named model never carries over to another model of the family."""
    config_x["layers"]["default"]["material"]["shear_viscosity"]["reference_viscosity_pas"] = 5.0e21
    # The table names the "reference" law; a constant law is built from its own defaults.
    assert make_viscosity("reference").reference_viscosity == 5.0e21
    assert make_viscosity("constant").reference_viscosity == 1.0e22
    # The case that showed why: the isotope table's dataset must not set a fixed model's reference time.
    fixed = make_radiogenics("fixed").get_config_dict()
    assert fixed["ref_time_s"] == 0.0


def test_a_model_ignores_the_other_model_keys_of_a_merged_table():
    """A table merged family-wide can carry both models' keys; each model is built from its own."""
    merged = {"isotopes": "modern_day_chondritic", "fixed_heat_production_w_kg": 2.0e-12}
    fixed = make_radiogenics("constant", dict(merged)).get_config_dict()
    assert fixed["ref_time_s"] == 0.0 and fixed["fixed_heat_production_w_kg"] == 2.0e-12
    assert "isotope_names" not in fixed
    isotope = make_radiogenics("isotopes", dict(merged)).get_config_dict()
    assert isotope["isotope_names"] == ["U238", "U235", "Th232", "K40"] and isotope["ref_time_s"] > 1.0e17


def test_the_world_builder_gives_a_fixed_layer_the_fixed_keys_only():
    """A rock layer that names a fixed model over the material's isotope defaults gets no dataset reference time."""
    from TidalPy.structures_x import build_world
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
    """The default table is found under the model's canonical name, an alias, or another case."""
    model = make_radiogenics(name)
    assert model.get_config_dict()["isotope_names"] == ["U238", "U235", "Th232", "K40"]
