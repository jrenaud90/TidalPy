"""Stellar luminosity model hierarchy (``TidalPy.stellar_x``).

Covers the factory (names, aliases, unknown-name ValueError, rich subclass), the three models
(Fixed, MassToLuminosity, PowerLaw), the shared Stefan-Boltzmann temperature conversions, scalar +
vectorized ``calc_luminosity``, the direct convenience functions, the config dict, the binary
round-trip, and the isinstance chain.
"""
import math
import os
import tempfile

import numpy as np
import pytest
import scipy.constants

from TidalPy.stellar_x import (
    LuminosityBase,
    FixedLuminosity,
    MassToLuminosity,
    PowerLawLuminosity,
    make_luminosity,
    fixed,
    mass_to_luminosity,
    power_law,
)

# Solar anchors used by the C++ models (TidalPyConstants::d_MASS_SOLAR / d_LUMINOSITY_SOLAR).
MASS_SOLAR = 1.988435e30
LUM_SOLAR = 3.848e26
SIGMA = scipy.constants.Stefan_Boltzmann


def _reference_mass_luminosity(mass_kg):
    """Independent Python re-implementation of the piecewise main-sequence L(M) relation."""
    ratio = mass_kg / MASS_SOLAR
    if ratio < 0.2:
        return LUM_SOLAR * 0.23 * ratio ** 2.3
    if ratio < 0.85:
        exponent = (-141.7 * ratio ** 4 + 232.4 * ratio ** 3
                    - 129.1 * ratio ** 2 + 33.29 * ratio + 0.215)
        return LUM_SOLAR * ratio ** exponent
    if ratio < 2.0:
        return LUM_SOLAR * ratio ** 4
    if ratio < 20.0:
        return LUM_SOLAR * 1.4 * ratio ** 3.5
    return LUM_SOLAR * 3.2e4 * ratio


# =====================================================================================================================
# Factory
# =====================================================================================================================
@pytest.mark.parametrize("name, cls", [
    ("fixed", FixedLuminosity),
    ("constant", FixedLuminosity),
    ("mass_to_luminosity", MassToLuminosity),
    ("cuntz_wang", MassToLuminosity),
    ("cw", MassToLuminosity),
    ("CW", MassToLuminosity),          # case-insensitive
    ("power_law", PowerLawLuminosity),
    ("powerlaw", PowerLawLuminosity),
])
def test_factory_names_and_aliases(name, cls):
    model = make_luminosity(name)
    assert isinstance(model, cls)
    assert isinstance(model, LuminosityBase)


def test_factory_unknown_name_raises():
    with pytest.raises(ValueError):
        make_luminosity("not_a_model")


def test_base_is_abstract():
    with pytest.raises(TypeError):
        LuminosityBase()


# =====================================================================================================================
# MassToLuminosity - piecewise relation
# =====================================================================================================================
@pytest.mark.parametrize("ratio", [0.05, 0.1, 0.19, 0.3, 0.5, 0.84, 0.9, 1.0, 1.5, 1.99, 5.0, 19.0, 25.0, 50.0])
def test_mass_to_luminosity_matches_reference(ratio):
    model = MassToLuminosity()
    mass = ratio * MASS_SOLAR
    assert math.isclose(model.calc_luminosity(mass), _reference_mass_luminosity(mass), rel_tol=1e-12)


def test_mass_to_luminosity_solar_is_solar():
    """At exactly one solar mass the (M/Msun)^4 branch gives L = Lsun."""
    model = MassToLuminosity()
    assert math.isclose(model.calc_luminosity(MASS_SOLAR), LUM_SOLAR, rel_tol=1e-12)


def test_mass_to_luminosity_monotonic():
    model = MassToLuminosity()
    masses = np.linspace(0.05, 50.0, 200) * MASS_SOLAR
    lums = model.calc_luminosity(masses)
    assert np.all(np.diff(lums) > 0.0)


def test_mass_to_luminosity_nonpositive_mass_is_nan():
    model = MassToLuminosity()
    assert np.isnan(model.calc_luminosity(0.0))
    assert np.isnan(model.calc_luminosity(-1.0e30))


# =====================================================================================================================
# FixedLuminosity
# =====================================================================================================================
def test_fixed_luminosity_ignores_mass():
    model = make_luminosity("fixed", {"luminosity_w": 1.5e26})
    assert model.luminosity == 1.5e26
    assert model.calc_luminosity(1.0e30) == 1.5e26
    assert model.calc_luminosity(9.9e30) == 1.5e26


# =====================================================================================================================
# PowerLawLuminosity
# =====================================================================================================================
def test_power_law_defaults_and_values():
    model = PowerLawLuminosity()  # coeff 1, exponent 3.5
    assert model.coeff == 1.0
    assert model.exponent == 3.5
    mass = 3.0 * MASS_SOLAR
    assert math.isclose(model.calc_luminosity(mass), LUM_SOLAR * 3.0 ** 3.5, rel_tol=1e-12)


def test_power_law_custom_params():
    model = make_luminosity("power_law", {"power_law_coeff": 0.23, "power_law_exponent": 2.3})
    mass = 0.1 * MASS_SOLAR
    assert math.isclose(model.calc_luminosity(mass), LUM_SOLAR * 0.23 * 0.1 ** 2.3, rel_tol=1e-12)


# =====================================================================================================================
# Stefan-Boltzmann conversions (shared on the base)
# =====================================================================================================================
def test_stefan_boltzmann_absolute():
    model = MassToLuminosity()
    radius, temperature = 6.957e8, 5772.0
    expected = 4.0 * math.pi * radius ** 2 * SIGMA * temperature ** 4
    assert math.isclose(model.calc_luminosity_from_temperature(temperature, radius), expected, rel_tol=1e-12)


def test_stefan_boltzmann_round_trip():
    model = MassToLuminosity()
    radius, temperature = 6.957e8, 5772.0
    lum = model.calc_luminosity_from_temperature(temperature, radius)
    assert math.isclose(model.calc_temperature_from_luminosity(lum, radius), temperature, rel_tol=1e-12)


def test_effective_temperature_from_mass():
    """calc_effective_temperature == temperature_from_luminosity(calc_luminosity(mass))."""
    model = MassToLuminosity()
    radius = 6.957e8
    mass = 1.2 * MASS_SOLAR
    lum = model.calc_luminosity(mass)
    assert math.isclose(model.calc_effective_temperature(mass, radius),
                        model.calc_temperature_from_luminosity(lum, radius), rel_tol=1e-14)


@pytest.mark.parametrize("temperature, radius", [(-1.0, 6.957e8), (5772.0, -1.0), (0.0, 6.957e8)])
def test_stefan_boltzmann_bad_inputs_nan(temperature, radius):
    model = MassToLuminosity()
    assert np.isnan(model.calc_luminosity_from_temperature(temperature, radius))


# =====================================================================================================================
# Vectorization + direct convenience functions
# =====================================================================================================================
def test_calc_luminosity_vectorized_shape():
    model = MassToLuminosity()
    masses = np.array([0.1, 0.5, 1.0, 5.0]).reshape(2, 2) * MASS_SOLAR
    out = model.calc_luminosity(masses)
    assert isinstance(out, np.ndarray)
    assert out.shape == (2, 2)
    for mass, lum in zip(masses.ravel(), out.ravel()):
        assert math.isclose(lum, _reference_mass_luminosity(mass), rel_tol=1e-12)


def test_direct_functions_match_models():
    mass = 0.5 * MASS_SOLAR
    assert math.isclose(mass_to_luminosity(mass), MassToLuminosity().calc_luminosity(mass), rel_tol=1e-14)
    assert fixed(mass, 2.0e26) == 2.0e26
    assert math.isclose(power_law(mass, 1.0, 4.0), LUM_SOLAR * 0.5 ** 4, rel_tol=1e-12)


def test_direct_function_broadcast():
    masses = np.array([0.5, 1.0, 2.0]) * MASS_SOLAR
    out = mass_to_luminosity(masses)
    assert isinstance(out, np.ndarray)
    assert out.shape == (3,)


# =====================================================================================================================
# Config dict + binary round-trip
# =====================================================================================================================
def test_config_dict():
    fixed_model = make_luminosity("fixed", {"luminosity_w": 3.0e26})
    assert fixed_model.get_config_dict()["model_name"] == "fixed"
    assert fixed_model.get_config_dict()["luminosity_w"] == 3.0e26

    power = make_luminosity("power_law", {"power_law_coeff": 1.4, "power_law_exponent": 3.5})
    d = power.get_config_dict()
    assert d["model_name"] == "power_law"
    assert d["coeff"] == 1.4
    assert d["exponent"] == 3.5


@pytest.mark.parametrize("model, config", [
    ("fixed", {"luminosity_w": 3.0e26}),
    ("mass_to_luminosity", {}),
    ("power_law", {"power_law_coeff": 1.4, "power_law_exponent": 3.2}),
])
def test_binary_round_trip(model, config):
    original = make_luminosity(model, config)
    mass = 0.7 * MASS_SOLAR
    before = original.calc_luminosity(mass)

    path = os.path.join(tempfile.gettempdir(), f"tpy_lum_{model}.tpyb")
    original.save_binary(path)
    try:
        restored = make_luminosity(model, config)
        restored.load_binary(path)
        after = restored.calc_luminosity(mass)
        if np.isnan(before):
            assert np.isnan(after)
        else:
            assert math.isclose(after, before, rel_tol=1e-14)
        assert restored.get_config_dict() == original.get_config_dict()
    finally:
        if os.path.isfile(path):
            os.remove(path)
