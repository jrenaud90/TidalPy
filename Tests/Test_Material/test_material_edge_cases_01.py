"""Edge cases of the material laws: the modified polytrope at zero pressure, a NaN solid viscosity through melt
weakening, a rheology on the liquid phase of a two-phase material, and Material.replace by either spelling."""
import math

import pytest

from TidalPy.Material import Material, Phase, load_material
from TidalPy.Material.laws import make_eos
from TidalPy.PartialMelt.melting import make_melt_weakening

REFERENCE_DENSITY = 3000.0
POLYTROPE_COEFFICIENT = 1.0e-6


@pytest.mark.parametrize("exponent, expected", [
    (0.5, 0.0), (1.0, REFERENCE_DENSITY / POLYTROPE_COEFFICIENT), (1.5, math.inf)], ids=["below_1", "at_1", "above_1"])
def test_the_modified_polytrope_bulk_modulus_at_zero_pressure_is_its_limit(exponent, expected):
    """K = rho P^(1 - n) / (n c) falls to 0, holds at rho0 / c, or grows without bound as P goes to 0."""
    law = make_eos("modified_polytrope", {"reference_density_kg_m3": REFERENCE_DENSITY,
                                          "polytrope_coefficient": POLYTROPE_COEFFICIENT,
                                          "polytrope_exponent": exponent})
    assert law.calc_eos(0.0)["bulk_modulus"] == expected
    if math.isfinite(expected) and expected > 0.0:
        # Continuous from above.
        assert law.calc_eos(1.0)["bulk_modulus"] == pytest.approx(expected, rel=1.0e-3)


def test_a_nan_solid_viscosity_stays_nan_through_weakening():
    weakening = make_melt_weakening("henning")
    shear, viscosity = weakening.calc_weakening(
        temperature=1650.0, solidus=1600.0, liquidus=2000.0, solid_shear=5.0e10, solid_viscosity=math.nan,
        liquid_shear=0.0, liquid_viscosity=0.1)
    assert math.isnan(viscosity)
    assert shear > 0.0


def _two_phase_config():
    config = load_material("peridotite").get_config_dict()
    for key in ("schema_version", "description", "category"):
        config.pop(key, None)
    return config


def test_a_rheology_on_the_liquid_of_a_two_phase_material_is_refused():
    config = _two_phase_config()
    config["liquid"]["shear_rheology"] = {"model": "maxwell"}
    with pytest.raises(ValueError, match="liquid phase has a 'shear_rheology'"):
        load_material(config)


def test_a_liquid_only_material_keeps_its_rheology():
    liquid = Phase(eos={"model": "constant", "reference_density_kg_m3": 1000.0},
                   shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e-3},
                   shear_rheology={"model": "maxwell"})
    assert Material(liquid=liquid).liquid.shear_rheology is not None


@pytest.mark.parametrize("key", ["latent_heat", "latent_heat_j_kg"])
def test_replace_takes_either_spelling_of_a_parameter(key):
    changed = load_material("peridotite").replace(**{key: 2.0e5})
    assert changed.get_parameter("latent_heat_j_kg") == 2.0e5


def test_a_parameter_given_twice_names_both_spellings():
    with pytest.raises(TypeError, match="'latent_heat' and 'latent_heat_j_kg'"):
        load_material("peridotite").replace(latent_heat=1.0, latent_heat_j_kg=2.0)
