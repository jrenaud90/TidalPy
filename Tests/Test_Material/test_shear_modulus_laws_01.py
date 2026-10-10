"""Shear-modulus laws (TidalPy.Material.laws) against closed forms; the generic spec behavior is tested in
Tests/Test_Utilities/Test_Classes/test_spec_models_01.py."""

import math

import numpy as np
import pytest

from TidalPy.Material.laws import ConstantShearModulus, LinearShearModulus, make_shear_modulus


def test_constant():
    assert ConstantShearModulus(6.0e10).calc_shear_modulus(1.0e10, 1500.0) == 6.0e10


@pytest.mark.parametrize("pressure, temperature", [(0.0, 300.0), (2.0e10, 1500.0), (1.0e11, 3000.0)])
def test_linear(pressure, temperature):
    law = LinearShearModulus(
        shear_modulus=6.0e10,
        pressure_derivative=1.5,
        temperature_derivative=-1.0e7,
        reference_pressure=1.0e9,
        reference_temperature=300.0)
    expected = 6.0e10 + 1.5 * (pressure - 1.0e9) - 1.0e7 * (temperature - 300.0)
    assert law.calc_shear_modulus(pressure, temperature) == pytest.approx(expected, rel=1e-15)


def test_linear_without_temperature_drops_the_temperature_term():
    law = make_shear_modulus("linear", {"shear_modulus_pa": 6.0e10, "temperature_derivative_pa_k": -1.0e7})
    assert law.calc_shear_modulus(0.0) == 6.0e10
    assert law.calc_shear_modulus(0.0, math.nan) == 6.0e10


def test_interpolated():
    law = make_shear_modulus("interp", {"radius_m": [0.0, 2.0e6], "shear_modulus_pa": [0.0, 8.0e10]})
    radii = np.array([0.0, 5.0e5, 2.0e6, 3.0e6])
    assert np.allclose(law.calc_shear_modulus(0.0, radius=radii), [0.0, 2.0e10, 8.0e10, 8.0e10])
