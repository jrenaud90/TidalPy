"""Every spec-driven physics model, tested generically from its own parameter table.

A family built on parameter specs joins by adding one row to SPEC_FAMILIES: its factory, its model-name listing, and
its base class. Each of its models is then built from defaults and checked for its config and binary round trips,
unknown-key and bounds errors, with_parameters, parameter access, and copies. Physics is tested per family.
"""

import copy
import math

import pytest

from TidalPy.Material import laws
from TidalPy.Material import material as material_module
from TidalPy.PartialMelt import melting
from TidalPy.Rheology import rheology
from TidalPy.Utilities.classes.classes import PhysicsBase, TidalPyBaseClass
from TidalPy.Viscosity import viscosity as viscosity_module

# (family, factory, model-name listing, family base class)
SPEC_FAMILIES = [
    ("viscosity", viscosity_module.make_viscosity, viscosity_module.viscosity_model_names,
     viscosity_module.ViscosityBase),
    ("rheology", rheology.make_rheology, rheology.rheology_model_names, rheology.RheologyBase),
    ("eos", laws.make_eos, laws.eos_model_names, laws.EOSBase),
    ("shear_modulus", laws.make_shear_modulus, laws.shear_modulus_model_names, laws.ShearModulusBase),
    ("melting_curve", melting.make_melting_curve, melting.melting_curve_model_names, melting.MeltingCurveBase),
    ("melt_weakening", melting.make_melt_weakening, melting.melt_weakening_model_names, melting.MeltWeakeningBase),
    ("bulk_modulus_mixing", melting.make_bulk_modulus_mixing, melting.bulk_modulus_mixing_model_names,
     melting.BulkModulusMixingBase),
    ("bulk_viscosity_mixing", melting.make_bulk_viscosity_mixing, melting.bulk_viscosity_mixing_model_names,
     melting.BulkViscosityMixingBase),
    ("phase", material_module._PHASE_FAMILY.make, material_module._PHASE_FAMILY.model_names, material_module.Phase),
    ("material", material_module._MATERIAL_FAMILY.make, material_module._MATERIAL_FAMILY.model_names,
     material_module.Material),
]

CASES = [
    pytest.param(family, factory, base, name, id=f"{family}-{name}")
    for family, factory, names, base in SPEC_FAMILIES
    for name in names()
]


def _same(first, second):
    """Equality that counts two NaNs (an unset parameter) as equal."""
    if isinstance(first, float) and isinstance(second, float) and math.isnan(first) and math.isnan(second):
        return True
    return first == second


def _valid_change(entry):
    """A value inside the parameter's bounds that differs from its default."""
    if entry["kind"] == "bool":
        return not entry["default"]
    if entry["kind"] == "int":
        return entry["default"] + 1
    if entry["kind"] == "list[float]":
        return None
    default = entry["default"]
    if entry["bounds"] == "unit interval":
        return 0.5 if default != 0.5 else 0.25
    if default is None or not math.isfinite(default) or default == 0.0:
        return 1.0
    return default * 2.0


def _invalid_value(entry):
    """A value outside the parameter's bounds, or None when every value is allowed."""
    return {
        "finite": math.inf,
        "positive": -1.0,
        "non-negative": -1.0,
        "unit interval": 2.0,
        "positive or infinite": -1.0,
    }.get(entry["bounds"])


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_defaults_build(family, factory, base, name):
    model = factory(name, {})
    assert model.model_name == name
    assert isinstance(model, base)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_parameter_info_is_consistent(family, factory, base, name):
    model = factory(name, {})
    info = model.get_parameter_info()
    names = [entry["name"] for entry in info]
    keys = [entry["key"] for entry in info]
    assert len(set(names)) == len(names)
    assert len(set(keys)) == len(keys)
    config = model.get_config_dict()
    for entry in info:
        assert entry["doc"], entry["key"]
        by_name = model.get_parameter(entry["name"])
        assert _same(by_name, model.get_parameter(entry["key"]))
        assert _same(getattr(model, entry["name"]), by_name)
        assert _same(model.parameters[entry["name"]], by_name)
        default_unset = entry["kind"] == "float" and isinstance(entry["default"], float) and \
            math.isnan(entry["default"])
        if entry["kind"] != "list[float]" and not default_unset:
            assert config[entry["key"]] == by_name


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_config_round_trip(family, factory, base, name):
    model = factory(name, {})
    config = model.get_config_dict()
    rebuilt = factory(name, config)
    assert rebuilt.get_config_dict() == config


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_binary_round_trip(family, factory, base, name, tmp_path):
    model = factory(name, {})
    changed = model
    for entry in model.get_parameter_info():
        value = _valid_change(entry)
        if value is not None:
            changed = changed.with_parameters(**{entry["key"]: value})
    path = str(tmp_path / f"{family}_{name}.tpyb")
    changed.save_binary(path)
    reloaded = factory(name, {})
    reloaded.load_binary(path)
    assert reloaded.get_config_dict() == changed.get_config_dict()


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_unknown_key_names_closest(family, factory, base, name):
    info = factory(name, {}).get_parameter_info()
    if not info:
        pytest.skip("model has no parameters")
    key = info[0]["key"]
    misspelled = key[:-1] + ("x" if key[-1] != "x" else "y")
    with pytest.raises(ValueError, match=f"did you mean '{key}'"):
        factory(name, {misspelled: 1.0})


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_out_of_bounds_rejected(family, factory, base, name):
    model = factory(name, {})
    for entry in model.get_parameter_info():
        invalid = _invalid_value(entry)
        if entry["kind"] == "bool":
            invalid = 0.5
        if invalid is None:
            continue
        with pytest.raises(ValueError, match=entry["key"]):
            factory(name, {entry["key"]: invalid})


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_with_parameters_returns_a_new_model(family, factory, base, name):
    model = factory(name, {})
    before = model.get_config_dict()
    for entry in model.get_parameter_info():
        value = _valid_change(entry)
        if value is None:
            continue
        changed = model.with_parameters(**{entry["name"]: value})
        assert type(changed) is type(model)
        assert changed.get_parameter(entry["name"]) == value
        assert model.get_config_dict() == before


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_copies(family, factory, base, name):
    model = factory(name, {})
    assert copy.copy(model) is model
    duplicate = copy.deepcopy(model)
    assert duplicate is not model
    assert type(duplicate) is type(model)
    assert duplicate.get_config_dict() == model.get_config_dict()


@pytest.mark.parametrize("family,factory,base,name", CASES)
def test_unknown_model_name_names_closest(family, factory, base, name):
    with pytest.raises(ValueError, match=f"did you mean '{name}'"):
        factory(name + "x", {})
