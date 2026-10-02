"""A spec model's parameter has two spellings, its argument name and its config key. A keyword overrides a config
value under either spelling, and one table giving both spellings is refused rather than resolved by key order."""
import pytest

from TidalPy.Material import Material, Phase
from TidalPy.Utilities.classes import canonical_parameter_keys
from TidalPy.Viscosity import ConstantViscosity, make_viscosity


@pytest.mark.parametrize("config_key, keyword", [
    ("reference_viscosity_pas", "reference_viscosity"),
    ("reference_viscosity", "reference_viscosity_pas"),
])
def test_a_keyword_overrides_the_config_under_either_spelling(config_key, keyword):
    model = ConstantViscosity(config={config_key: 1.0e20}, **{keyword: 1.0e10})
    assert model.get_parameter("reference_viscosity") == 1.0e10


def test_composites_take_keywords_over_their_config():
    assert Phase(config={"thermal_conductivity_w_mk": 3.0}, thermal_conductivity=5.0).thermal_conductivity == 5.0
    assert Material(config={"solid": {}, "latent_heat_j_kg": 4.0e5}, latent_heat=1.0).latent_heat == 1.0


def test_one_table_giving_both_spellings_is_refused():
    with pytest.raises((TypeError, ValueError), match="same parameter|twice"):
        make_viscosity("constant", {"reference_viscosity_pas": 1.0e20, "reference_viscosity": 1.0e10})
    with pytest.raises(TypeError, match="same parameter"):
        canonical_parameter_keys(Material, {"latent_heat": 1.0, "latent_heat_j_kg": 2.0})


def test_canonical_keys_leave_other_keys_alone():
    table = canonical_parameter_keys(Phase, {"thermal_conductivity": 3.0, "eos": {"model": "constant"}})
    assert table == {"thermal_conductivity_w_mk": 3.0, "eos": {"model": "constant"}}
