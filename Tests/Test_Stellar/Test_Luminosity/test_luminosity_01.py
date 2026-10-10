"""Stellar luminosity models (fixed, mass-to-luminosity, power law): physics, Stefan-Boltzmann, factory, and I/O."""
import math

import numpy as np
import pytest
import scipy.constants

from TidalPy.Stellar import (
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
LUM_SOLAR = 3.828e26
SIGMA = scipy.constants.Stefan_Boltzmann
RADIUS_SOLAR = 6.957e8


def _reference_mass_luminosity(mass):
    """Independent re-implementation of the piecewise main-sequence L(M) relation."""
    ratio = mass / MASS_SOLAR
    if ratio < 0.2:
        return LUM_SOLAR * 0.23 * ratio ** 2.3
    if ratio < 0.85:
        exponent = (-141.7 * ratio ** 4 + 232.4 * ratio ** 3
                    - 129.1 * ratio ** 2 + 33.29 * ratio + 0.215)
        return LUM_SOLAR * ratio ** exponent
    if ratio < 2.0:
        return LUM_SOLAR * ratio ** 4
    # The linear branch takes over where it meets 1.4 M^3.5 (about 55 Msun), keeping continuity.
    if ratio < 55.0:
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
    ("CW", MassToLuminosity),
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
# MassToLuminosity
# =====================================================================================================================
@pytest.mark.parametrize("boundary_ratio, tolerance", [(0.2, 0.30), (0.85, 0.02), (2.0, 0.02), (55.0, 0.05)])
def test_mass_luminosity_branch_continuity(boundary_ratio, tolerance):
    """Adjacent branches of the reference and the C++ L(M) agree at their boundaries."""
    # A one-part-per-million offset selects each branch, so this measures the jump, not the local slope.
    epsilon = 1.0e-6
    below_mass = (1.0 - epsilon) * boundary_ratio * MASS_SOLAR
    above_mass = (1.0 + epsilon) * boundary_ratio * MASS_SOLAR
    for luminosity in (_reference_mass_luminosity, MassToLuminosity().calc_luminosity):
        below = luminosity(below_mass)
        assert abs(luminosity(above_mass) - below) / below < tolerance


@pytest.mark.parametrize("ratio", [0.05, 0.1, 0.19, 0.3, 0.5, 0.84, 0.9, 1.0, 1.5, 1.99, 5.0, 19.0, 25.0, 50.0])
def test_mass_to_luminosity_matches_reference(ratio):
    mass = ratio * MASS_SOLAR
    assert math.isclose(MassToLuminosity().calc_luminosity(mass), _reference_mass_luminosity(mass), rel_tol=1e-12)


def test_mass_to_luminosity_solar_is_solar():
    assert math.isclose(MassToLuminosity().calc_luminosity(MASS_SOLAR), LUM_SOLAR, rel_tol=1e-12)


def test_mass_to_luminosity_monotonic():
    lums = MassToLuminosity().calc_luminosity(np.linspace(0.05, 50.0, 200) * MASS_SOLAR)
    assert np.all(np.diff(lums) > 0.0)


@pytest.mark.parametrize("mass", [0.0, -1.0e30])
def test_mass_to_luminosity_nonpositive_mass_is_nan(mass):
    assert np.isnan(MassToLuminosity().calc_luminosity(mass))


# =====================================================================================================================
# FixedLuminosity and PowerLawLuminosity
# =====================================================================================================================
def test_fixed_luminosity_ignores_mass():
    model = make_luminosity("fixed", {"luminosity_w": 1.5e26})
    assert model.luminosity == 1.5e26
    assert model.calc_luminosity(1.0e30) == 1.5e26
    assert model.calc_luminosity(9.9e30) == 1.5e26


def test_power_law_defaults_and_values():
    model = PowerLawLuminosity()
    assert model.coeff == 1.0
    assert model.exponent == 3.5
    assert math.isclose(model.calc_luminosity(3.0 * MASS_SOLAR), LUM_SOLAR * 3.0 ** 3.5, rel_tol=1e-12)


def test_power_law_custom_params():
    model = make_luminosity("power_law", {"power_law_coeff": 0.23, "power_law_exponent": 2.3})
    assert math.isclose(model.calc_luminosity(0.1 * MASS_SOLAR), LUM_SOLAR * 0.23 * 0.1 ** 2.3, rel_tol=1e-12)


# =====================================================================================================================
# Stefan-Boltzmann conversions (shared on the base)
# =====================================================================================================================
def test_stefan_boltzmann_absolute():
    expected = 4.0 * math.pi * RADIUS_SOLAR ** 2 * SIGMA * 5772.0 ** 4
    luminosity = MassToLuminosity().calc_luminosity_from_temperature(5772.0, RADIUS_SOLAR)
    assert math.isclose(luminosity, expected, rel_tol=1e-12)


def test_stefan_boltzmann_round_trip():
    model = MassToLuminosity()
    luminosity = model.calc_luminosity_from_temperature(5772.0, RADIUS_SOLAR)
    assert math.isclose(model.calc_temperature_from_luminosity(luminosity, RADIUS_SOLAR), 5772.0, rel_tol=1e-12)


def test_effective_temperature_from_mass():
    """calc_effective_temperature is temperature_from_luminosity of calc_luminosity."""
    model = MassToLuminosity()
    mass = 1.2 * MASS_SOLAR
    luminosity = model.calc_luminosity(mass)
    assert math.isclose(model.calc_effective_temperature(mass, RADIUS_SOLAR),
                        model.calc_temperature_from_luminosity(luminosity, RADIUS_SOLAR), rel_tol=1e-14)


@pytest.mark.parametrize("temperature, radius", [(-1.0, RADIUS_SOLAR), (5772.0, -1.0), (0.0, RADIUS_SOLAR)])
def test_stefan_boltzmann_bad_inputs_nan(temperature, radius):
    assert np.isnan(MassToLuminosity().calc_luminosity_from_temperature(temperature, radius))


# =====================================================================================================================
# Vectorization and direct convenience functions
# =====================================================================================================================
def test_calc_luminosity_vectorized_shape():
    masses = np.array([0.1, 0.5, 1.0, 5.0]).reshape(2, 2) * MASS_SOLAR
    out = MassToLuminosity().calc_luminosity(masses)
    assert isinstance(out, np.ndarray)
    assert out.shape == (2, 2)
    for mass, luminosity in zip(masses.ravel(), out.ravel()):
        assert math.isclose(luminosity, _reference_mass_luminosity(mass), rel_tol=1e-12)


def test_direct_functions_match_models():
    mass = 0.5 * MASS_SOLAR
    assert math.isclose(mass_to_luminosity(mass), MassToLuminosity().calc_luminosity(mass), rel_tol=1e-14)
    assert fixed(mass, 2.0e26) == 2.0e26
    assert math.isclose(power_law(mass, 1.0, 4.0), LUM_SOLAR * 0.5 ** 4, rel_tol=1e-12)


def test_direct_function_broadcast():
    out = mass_to_luminosity(np.array([0.5, 1.0, 2.0]) * MASS_SOLAR)
    assert isinstance(out, np.ndarray)
    assert out.shape == (3,)


def test_read_only_mass_is_accepted():
    masses = np.array([0.5, 1.0, 2.0]) * MASS_SOLAR
    masses.setflags(write=False)
    expected = np.array([MassToLuminosity().calc_luminosity(float(mass)) for mass in masses])
    assert np.array_equal(MassToLuminosity().calc_luminosity(masses), expected)
    assert np.array_equal(mass_to_luminosity(masses), expected)


# =====================================================================================================================
# Config dict and binary round trip
# =====================================================================================================================
@pytest.mark.parametrize("model, config", [
    ("fixed", {"luminosity_w": 3.0e26}),
    ("power_law", {"power_law_coeff": 1.4, "power_law_exponent": 3.5}),
])
def test_config_dict(model, config):
    emitted = make_luminosity(model, config).get_config_dict()
    assert emitted["model"] == model
    for key, value in config.items():
        assert emitted[key] == value


@pytest.mark.parametrize("model, config", [
    ("fixed", {"luminosity_w": 3.0e26}),
    ("mass_to_luminosity", {}),
    ("power_law", {"power_law_coeff": 1.4, "power_law_exponent": 3.2}),
])
def test_binary_round_trip(model, config, tmp_path):
    original = make_luminosity(model, config)
    mass = 0.7 * MASS_SOLAR
    before = original.calc_luminosity(mass)
    path = str(tmp_path / f"tpy_lum_{model}.tpyb")
    original.save_binary(path)
    restored = make_luminosity(model, config)
    restored.load_binary(path)
    after = restored.calc_luminosity(mass)
    if np.isnan(before):
        assert np.isnan(after)
    else:
        assert math.isclose(after, before, rel_tol=1e-14)
    assert restored.get_config_dict() == original.get_config_dict()
