"""Phases and materials (TidalPy.Material): construction, validation, the evaluator under each physics switch, and
round trips. The laws inside are tested in their own files; here the composition is."""

import math

import numpy as np
import pytest

from TidalPy.Material import Material, Phase, make_eos, make_material, make_phase
from TidalPy.Material.laws import ConstantEOS, MurnaghanEOS
from TidalPy.PartialMelt.melting import HenningMeltWeakening
from TidalPy.Viscosity import ConstantViscosity

_T_SOL, _T_LIQ = 1600.0, 2000.0


def _solid():
    return Phase(
        eos={"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 1.3e11,
             "thermal_expansion_1_k": 3.0e-5},
        shear_modulus={"model": "constant", "shear_modulus_pa": 7.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e21},
        bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e22},
        thermal_conductivity=4.0,
        heat_capacity=1200.0)


def _liquid():
    return Phase(
        eos={"model": "constant", "reference_density_kg_m3": 2800.0, "bulk_modulus_pa": 2.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 0.2},
        thermal_conductivity=2.0,
        heat_capacity=1500.0)


def _material(**slots):
    components = dict(
        solid=_solid(), liquid=_liquid(),
        solidus={"model": "constant", "temperature_k": _T_SOL},
        liquidus={"model": "constant", "temperature_k": _T_LIQ},
        latent_heat=4.0e5)
    components.update(slots)
    return Material(**components)


_MELTING = dict(use_melting=True)


# =====================================================================================================================
# Construction and validation
# =====================================================================================================================
def test_components_read_back_as_their_classes():
    material = _material(weakening="henning")
    assert isinstance(material.solid, Phase)
    assert isinstance(material.liquid.eos, ConstantEOS)
    assert isinstance(material.weakening, HenningMeltWeakening)
    assert isinstance(material.solid.shear_viscosity, ConstantViscosity)
    assert material.bulk_modulus_mixing is None
    assert material.can_melt


def test_defaults_build():
    assert isinstance(Phase().eos, ConstantEOS)
    assert not Material().can_melt
    assert not Material().is_liquid_only
    assert Phase().shear_modulus is None


@pytest.mark.parametrize("bad, message", [
    (dict(solidus=None), "needs a 'solidus'"),
    (dict(liquid=Phase(eos="murnaghan")), "needs a 'shear_viscosity'"),
])
def test_material_validation(bad, message):
    with pytest.raises(ValueError, match=message):
        _material(**bad)


def test_melting_laws_need_both_phases():
    with pytest.raises(ValueError, match="not both a solid"):
        Material(_solid(), weakening="henning")
    with pytest.raises(ValueError, match="not both a solid"):
        Material(liquid=_liquid(), solidus={"model": "constant", "temperature_k": _T_SOL})


def test_a_liquid_only_material_is_liquid_everywhere():
    material = Material(liquid=_liquid())
    assert material.is_liquid_only and not material.can_melt
    assert material.solid is None
    for switches in ({}, _MELTING, dict(use_melting=True, use_melt_density=True)):
        state = material.calc_state(1.0e9, 300.0, **switches)
        assert state["phase"] == "liquid"
        assert state["melt_fraction"] == 1.0
        assert state["density"] == 2800.0
        assert state["shear_modulus"] == 0.0
        assert state["shear_viscosity"] == 0.2
    assert material.get_config_dict().keys() == {"liquid", "latent_heat_j_kg"}


def test_wrong_family_in_a_slot():
    with pytest.raises(ValueError, match="takes a equation of state model"):
        Phase(eos=ConstantViscosity())


def test_misspelled_slot_names_the_closest():
    with pytest.raises(ValueError, match="did you mean 'shear_viscosity'"):
        make_phase({"shear_viscosty": {"model": "constant"}})


# =====================================================================================================================
# The evaluator
# =====================================================================================================================
def test_solid_material_is_its_solid_phase():
    state = _material().calc_state(1.0e9, 1800.0)
    assert state["phase"] == "solid"
    assert state["melt_fraction"] == 0.0
    assert state["density"] == 3300.0
    assert state["shear_modulus"] == 7.0e10
    assert math.isnan(state["solidus"])


def test_thermal_expansion_switch():
    material = _material()
    assert material.calc_state(0.0, 1300.0)["density"] == 3300.0
    assert material.calc_state(0.0, 1300.0, use_thermal_expansion=True)["density"] == pytest.approx(
        3300.0 * math.exp(-3.0e-5 * 1000.0), rel=1e-14)


@pytest.mark.parametrize("temperature", [1500.0, 1700.0, 1900.0, 2100.0])
def test_melt_fraction_and_phase(temperature):
    state = _material().calc_state(0.0, temperature, **_MELTING)
    phi = min(max((temperature - _T_SOL) / (_T_LIQ - _T_SOL), 0.0), 1.0)
    assert state["melt_fraction"] == pytest.approx(phi)
    assert state["phase"] == ("solid" if phi == 0.0 else "liquid" if phi == 1.0 else "partial")
    assert state["solidus"] == _T_SOL and state["liquidus"] == _T_LIQ


def test_no_weakening_keeps_the_solid_until_fully_molten():
    material = _material()
    partial = material.calc_state(0.0, 1900.0, **_MELTING)
    assert partial["shear_modulus"] == 7.0e10
    assert partial["shear_viscosity"] == 1.0e21
    assert partial["bulk_modulus"] == 1.3e11
    molten = material.calc_state(0.0, 2100.0, **_MELTING)
    assert molten["shear_modulus"] == 0.0
    assert molten["shear_viscosity"] == 0.2
    assert molten["bulk_modulus"] == 2.0e10
    assert math.isnan(molten["bulk_viscosity"])


def test_density_mixes_only_with_its_switch():
    material = _material()
    assert material.calc_state(0.0, 1800.0, **_MELTING)["density"] == 3300.0
    assert material.calc_state(0.0, 1800.0, use_melting=True, use_melt_density=True)["density"] == pytest.approx(
        0.5 * 3300.0 + 0.5 * 2800.0)


def test_thermal_properties_mix_and_latent_heat():
    state = _material().calc_state(0.0, 1800.0, **_MELTING)
    assert state["thermal_conductivity"] == pytest.approx(3.0)
    assert state["heat_capacity"] == pytest.approx(0.5 * 1200.0 + 0.5 * 1500.0 + 4.0e5 / (_T_LIQ - _T_SOL))


def test_a_single_melting_temperature_is_a_step():
    material = _material(liquidus={"model": "constant", "temperature_k": _T_SOL})
    assert material.calc_state(0.0, _T_SOL - 1.0, **_MELTING)["phase"] == "solid"
    above = material.calc_state(0.0, _T_SOL + 1.0, **_MELTING)
    assert above["phase"] == "liquid"
    assert above["heat_capacity"] == pytest.approx(1500.0)


def test_pressure_melting_switch():
    material = _material(solidus={"model": "simon_glatzel", "temperature_k": _T_SOL, "simon_a_pa": 1.0e9,
                                  "simon_c": 5.0})
    assert material.calc_state(5.0e9, 1700.0, **_MELTING)["solidus"] == _T_SOL
    assert material.calc_state(5.0e9, 1700.0, use_melting=True, use_pressure_melting=True)["solidus"] == \
        pytest.approx(_T_SOL * 6.0 ** 0.2)


def test_mixing_laws():
    material = _material(bulk_modulus_mixing="hashin_shtrikman",
                         bulk_viscosity_mixing={"model": "compaction", "coefficient": 1.0, "exponent": 1.0})
    state = material.calc_state(0.0, 1800.0, **_MELTING)
    k_s, k_l, phi = 1.3e11, 2.0e10, 0.5
    mu = state["shear_modulus"]
    hashin_shtrikman = k_s + phi / (1.0 / (k_l - k_s) + (1.0 - phi) / (k_s + 4.0 / 3.0 * mu))
    assert state["bulk_modulus"] == pytest.approx(hashin_shtrikman)
    assert state["bulk_viscosity"] == pytest.approx(1.0 / (1.0 / 1.0e22 + phi / state["shear_viscosity"]))


def test_without_a_temperature_there_is_no_melt_state():
    state = _material().calc_state(0.0, math.nan, **_MELTING)
    assert math.isnan(state["melt_fraction"])
    assert state["phase"] == "solid"


def test_temperature_dependent_thermal_properties():
    phase = Phase(thermal_conductivity=2.0, conductivity_temperature_exponent=-1.0, heat_capacity=2000.0,
                  heat_capacity_temperature_exponent=1.0, thermal_reference_temperature=250.0)
    state = phase.calc_state(0.0, 125.0)
    assert state["thermal_conductivity"] == pytest.approx(4.0)
    assert state["heat_capacity"] == pytest.approx(1000.0)


def test_vectorized_evaluation():
    temperatures = np.linspace(1500.0, 2100.0, 13)
    state = _material().calc_state(np.array([0.0, 1.0e9])[:, None], temperatures[None, :], **_MELTING)
    assert state["density"].shape == (2, 13)
    assert state["phase"].shape == (2, 13)
    assert state["melt_fraction"][0, -1] == 1.0


# =====================================================================================================================
# Round trips and functional updates
# =====================================================================================================================
def test_config_and_binary_round_trips(tmp_path):
    material = _material(weakening="spohn", bulk_modulus_mixing="hs")
    config = material.get_config_dict()
    assert set(config) == {"solid", "liquid", "melting", "latent_heat_j_kg"}
    assert set(config["melting"]) == {"solidus", "liquidus", "weakening", "bulk_modulus_mixing"}
    assert make_material(config).get_config_dict() == config
    path = str(tmp_path / "material.tpyb")
    material.save_binary(path)
    loaded = Material()
    loaded.load_binary(path)
    assert loaded.get_config_dict() == config


def test_replace_and_with_parameters_leave_the_original():
    material = _material()
    swapped = material.replace(liquid=Phase(eos=MurnaghanEOS(), shear_viscosity=ConstantViscosity(1.0)))
    assert isinstance(swapped.liquid.eos, MurnaghanEOS)
    assert isinstance(material.liquid.eos, ConstantEOS)
    assert swapped.latent_heat == 4.0e5
    hotter = material.with_parameters(latent_heat=1.0e5)
    assert hotter.latent_heat == 1.0e5 and material.latent_heat == 4.0e5
    assert make_eos("constant", {}).model_name == "constant"
