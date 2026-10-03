"""Viscosity laws: the reference law's anchor pressure, the radius-tabulated law, and the constructor conventions."""

import math

import numpy as np
import pytest
from scipy.constants import R as _R   # molar gas constant [J/mol/K], the same source as the C++ config

from TidalPy.Viscosity import (
    ArrheniusViscosity,
    CompositeViscosity,
    ConstantViscosity,
    InterpolatedViscosity,
    ReferenceViscosity,
    make_viscosity,
)


@pytest.mark.parametrize("temperature, pressure", [(1600.0, 2.0e10), (1400.0, 5.0e10), (2000.0, 0.0)])
def test_reference_law_anchored_at_a_pressure(temperature, pressure):
    """eta = eta_ref exp((E + P V) / (R T) - (E + P_ref V) / (R T_ref)), and eta_ref at (T_ref, P_ref)."""
    law = ReferenceViscosity(
        reference_viscosity=1.0e21,
        reference_temperature=1600.0,
        molar_activation_energy=3.0e5,
        molar_activation_volume=4.0e-6,
        reference_pressure=2.0e10)
    expected = 1.0e21 * math.exp(
        (3.0e5 + pressure * 4.0e-6) / (_R * temperature) - (3.0e5 + 2.0e10 * 4.0e-6) / (_R * 1600.0))
    assert law.calc_viscosity(temperature, pressure) == pytest.approx(expected, rel=1e-12)
    assert law.calc_viscosity(1600.0, 2.0e10) == pytest.approx(1.0e21, rel=1e-13)


def test_reference_temperature_must_be_positive():
    with pytest.raises(ValueError, match="reference_temperature_k"):
        make_viscosity("reference", {"reference_temperature_k": 0.0})


def test_interpolated_law_reads_the_radius():
    law = InterpolatedViscosity(radius=[0.0, 1.0e6], viscosity=[1.0e20, 3.0e20])
    assert law.calc_viscosity(1500.0, 0.0, 5.0e5) == pytest.approx(2.0e20)
    assert np.allclose(law.calc_viscosity(1500.0, 0.0, np.array([0.0, 2.0e6])), [1.0e20, 3.0e20])
    assert math.isnan(law.calc_viscosity(1500.0))


def test_positional_and_keyword_arguments():
    by_position = ReferenceViscosity(1.0e21, 1500.0)
    assert by_position.reference_viscosity == 1.0e21
    assert by_position.reference_temperature == 1500.0
    assert ReferenceViscosity(reference_viscosity_pas=2.0e21).reference_viscosity == 2.0e21
    with pytest.raises(TypeError, match="multiple values"):
        ReferenceViscosity(1.0e21, reference_viscosity=2.0e21)


def test_composite_adds_mechanisms_in_parallel():
    """1 / eta = sum of 1 / eta_i; a rigid mechanism adds nothing."""
    diffusion = ArrheniusViscosity(arrhenius_coeff=1.0e5, molar_activation_energy=6.0e4)
    dislocation = ArrheniusViscosity(arrhenius_coeff=1.0e-12, molar_activation_energy=1.8e5)
    composite = CompositeViscosity(diffusion, dislocation, ConstantViscosity(math.inf))
    for temperature in (200.0, 250.0, 270.0):
        expected = 1.0 / (1.0 / diffusion.calc_viscosity(temperature) + 1.0 / dislocation.calc_viscosity(temperature))
        assert composite.calc_viscosity(temperature) == pytest.approx(expected, rel=1e-14)
    assert [type(model) for model in composite.mechanisms] == [ArrheniusViscosity, ArrheniusViscosity, ConstantViscosity]


def test_composite_config_and_binary_round_trip(tmp_path):
    composite = make_viscosity("composite", {"mechanisms": [
        {"model": "constant", "reference_viscosity_pas": 1.0e20},
        {"model": "reference", "reference_viscosity_pas": 1.0e21}]})
    config = composite.get_config_dict()
    assert [table["model"] for table in config["mechanisms"]] == ["constant", "reference"]
    assert make_viscosity("parallel", config).get_config_dict() == config
    path = str(tmp_path / "composite.tpyb")
    composite.save_binary(path)
    loaded = CompositeViscosity()
    loaded.load_binary(path)
    assert loaded.get_config_dict() == config


def test_composite_rejects_a_model_of_another_family():
    from TidalPy.Rheology import Maxwell
    with pytest.raises(TypeError, match="viscosity model or a config table"):
        CompositeViscosity(Maxwell())
