"""Every physics family's factory takes back its models' config dicts and refuses a key its model does not read."""
import importlib

import pytest


# Each family: factory module, factory name, family label used in the error, and (model name, required config).
_FAMILIES = [
    ("TidalPy.Rheology.rheology", "make_rheology", "rheology",
     [("elastic", None), ("viscous", None), ("voigt", None), ("maxwell", None), ("burgers", None),
      ("andrade", None), ("sundberg", None)]),
    ("TidalPy.Viscosity.viscosity", "make_viscosity", "viscosity",
     [("arrhenius", None), ("reference", None), ("constant", None)]),
    ("TidalPy.Cooling.cooling", "make_cooling", "cooling",
     [("off", None), ("conduction", None), ("convection", None)]),
    ("TidalPy.Radiogenics.radiogenics", "make_radiogenics", "radiogenics",
     [("off", None), ("fixed", None), ("isotope", {"isotopes": "modern_day_chondritic"})]),
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


@pytest.mark.parametrize("module_path,factory_name,model_name,base_config", _model_cases())
def test_factory_accepts_its_own_config_dict(module_path, factory_name, model_name, base_config):
    """Every key a model emits is accepted back, so get_config_dict() round trips through the factory."""
    factory = _factory(module_path, factory_name)
    config = factory(model_name, base_config).get_config_dict()
    rebuilt = factory(config["model"], config)
    assert rebuilt.get_config_dict() == config


@pytest.mark.parametrize("module_path,factory_name,family,model_name,base_config", _family_cases())
def test_factory_rejects_unrecognized_key(module_path, factory_name, family, model_name, base_config):
    """A key the model does not read raises ValueError naming the family and the key."""
    factory = _factory(module_path, factory_name)
    config = dict(base_config or {})
    config["not_a_real_parameter"] = 1.0
    with pytest.raises(ValueError, match=f"{family} model .* has no parameter 'not_a_real_parameter'"):
        factory(model_name, config)


@pytest.mark.parametrize("module_path,factory_name,model_name,base_config", _model_cases())
def test_factory_names_the_closest_key(module_path, factory_name, model_name, base_config):
    """A misspelled key names the closest one the model reads."""
    factory = _factory(module_path, factory_name)
    info = factory(model_name, base_config).get_parameter_info()
    if not info:
        pytest.skip(f"the model '{model_name}' takes no parameters")
    key = info[0]["key"]
    config = {name: value for name, value in (base_config or {}).items() if name != key}
    config[key + "x"] = info[0]["default"]
    with pytest.raises(ValueError, match=f"did you mean '{key}'"):
        factory(model_name, config)
