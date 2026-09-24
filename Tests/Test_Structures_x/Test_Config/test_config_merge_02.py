"""Merging a model table over its defaults when the model changes, and the checks on switches and saved values.

A default model's own parameters never reach a different model, while nested tables and a material's
law-independent properties still merge; switches must be real booleans; numpy values are written as numbers.
"""
import numpy as np
import pytest
import toml

from TidalPy.configurations import merge_configs, plain_config
from TidalPy.structures_x.configs.toml_loader import validate_layer_config


def test_material_model_change_keeps_the_material():
    base = {"material": {
        "model": "constant",
        "reference_density_kg_m3": 917.0,
        "shear_modulus_static_pa": 3.3e9,
        "shear_viscosity": {"model": "arrhenius", "reference_viscosity_pas": 1.0e14}}}
    merged = merge_configs(base, {"material": {"model": "birch_murnaghan", "reference_bulk_modulus_pa": 1.0e10}})
    material = merged["material"]
    assert material["model"] == "birch_murnaghan"
    assert material["reference_density_kg_m3"] == 917.0
    assert material["shear_modulus_static_pa"] == 3.3e9
    assert material["shear_viscosity"]["reference_viscosity_pas"] == 1.0e14


def test_model_change_drops_the_old_models_parameters():
    base = {"radiogenics": {"model": "isotope", "isotopes": "modern_day_chondritic", "ref_time_s": 1.4e17}}
    merged = merge_configs(base, {"radiogenics": {"model": "fixed", "fixed_heat_production_w_kg": 1.0e-11}})
    assert merged["radiogenics"] == {"model": "fixed", "fixed_heat_production_w_kg": 1.0e-11}


def test_same_model_still_merges_key_by_key():
    base = {"radiogenics": {"model": "isotope", "isotopes": "modern_day_chondritic"}}
    merged = merge_configs(base, {"radiogenics": {"model": "isotope", "ref_time_s": 0.0}})
    assert merged["radiogenics"]["isotopes"] == "modern_day_chondritic"


@pytest.mark.parametrize("value", ("false", 0, 1.0))
def test_switches_must_be_booleans(value):
    with pytest.raises(ValueError, match="must be true or false"):
        validate_layer_config("mantle", {"class": "physics", "radius_fraction": 1.0, "is_solid": value})


def test_numpy_values_are_written_as_numbers():
    config = {
        "a": np.float64(4.2e8),
        "b": [np.int64(3), np.int64(2)],
        "c": np.array([1.0, 2.0]),
        "d": {"e": np.bool_(True)}}
    loaded = toml.loads(toml.dumps(plain_config(config)))
    assert loaded == {"a": 4.2e8, "b": [3, 2], "c": [1.0, 2.0], "d": {"e": True}}
