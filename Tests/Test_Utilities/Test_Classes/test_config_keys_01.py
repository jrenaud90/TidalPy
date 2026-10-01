"""Tests for physics model config key checking (``TidalPy.Utilities.classes.check_config_keys``)."""
import importlib

import pytest

from TidalPy.Utilities.classes import check_config_keys


# Each family: factory module, factory name, family label used in the error, and (model name, required config).
_FAMILIES = [
    ("TidalPy.Rheology.rheology", "make_rheology", "rheology",
     [("elastic", None), ("viscous", None), ("voigt", None), ("maxwell", None), ("burgers", None),
      ("andrade", None), ("sundberg", None)]),
    ("TidalPy.Viscosity.viscosity", "make_viscosity", "viscosity",
     [("arrhenius", None), ("reference", None), ("constant", None)]),
    ("TidalPy.PartialMelt.partial_melt", "make_partial_melt", "partial-melt",
     [("off", None), ("spohn", None), ("henning", None)]),
    ("TidalPy.Cooling.cooling", "make_cooling", "cooling",
     [("off", None), ("conduction", None), ("convection", None)]),
    ("TidalPy.Radiogenics.radiogenics", "make_radiogenics", "radiogenics",
     [("off", None), ("fixed", None), ("isotope", {"isotopes": "modern_day_chondritic"})]),
    ("TidalPy.Material.eos.material_eos", "make_material_eos", "material EOS",
     [("constant", None), ("bm", None), ("vinet", None),
      ("interpolate", {"radius_m": [0.0, 1.0e6, 2.0e6], "density_kg_m3": [5000.0, 4000.0, 3000.0],
                       "shear_modulus_pa": [1.0e11, 8.0e10, 6.0e10]})]),
    ("TidalPy.Stellar.luminosity", "make_luminosity", "luminosity",
     [("fixed", None), ("mass_to_luminosity", None), ("power_law", None)]),
    ("TidalPy.Tides.classes.tide", "make_tide", "tide",
     [("rheology", None), ("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]}),
      ("ctl", {"fixed_k": [0.3], "fixed_dt_s": [600.0]}),
      ("ctl_q", {"fixed_k": [0.3], "fixed_q": [50.0], "fixed_dt_s": [600.0]})]),
]


def _factory(module_path, factory_name):
    return getattr(importlib.import_module(module_path), factory_name)


def _model_cases():
    return [pytest.param(module_path, factory_name, model_name, base_config, id=f"{family}-{model_name}")
            for module_path, factory_name, family, models in _FAMILIES
            for model_name, base_config in models]


def _family_cases():
    return [pytest.param(module_path, factory_name, family, models[0][0], models[0][1], id=family)
            for module_path, factory_name, family, models in _FAMILIES]


@pytest.mark.parametrize(
    "config",
    [None, {}, {"model": "andrade", "alpha": 0.3}],
    ids=["none", "empty", "model_and_known_key"])
def test_empty_config_and_model_key_accepted(config):
    """An empty config and the model key are accepted without raising."""
    check_config_keys(config, {"alpha"}, "rheology")


def test_unrecognized_key_names_closest_accepted_key():
    """An unrecognized key raises ValueError naming the closest accepted key."""
    expected = r"unrecognized partial-melt config key.*'solidus' \(did you mean 'solidus_k'\?\)"
    with pytest.raises(ValueError, match=expected):
        check_config_keys({"solidus": 1500.0}, {"solidus_k", "liquidus_k"}, "partial-melt")


@pytest.mark.parametrize("module_path,factory_name,model_name,base_config", _model_cases())
def test_factory_accepts_its_own_config_dict(module_path, factory_name, model_name, base_config):
    """Every key a model emits is accepted back, so get_config_dict() round trips through the factory."""
    factory = _factory(module_path, factory_name)
    config = factory(model_name, base_config).get_config_dict()
    rebuilt = factory(config["model"], config)
    assert rebuilt.get_config_dict() == config


@pytest.mark.parametrize("module_path,factory_name,family,model_name,base_config", _family_cases())
def test_factory_rejects_unrecognized_key(module_path, factory_name, family, model_name, base_config):
    """A key no model in the family reads raises ValueError naming the family and the key."""
    factory = _factory(module_path, factory_name)
    config = dict(base_config or {})
    config["not_a_real_parameter"] = 1.0
    # A spec-driven family names the model; the rest name the family.
    with pytest.raises(
            ValueError,
            match=f"(unrecognized {family} config key|{family} model .* has no parameter).*not_a_real_parameter"):
        factory(model_name, config)
