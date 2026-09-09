"""Star-world attached luminosity model (``StarWorld`` + ``LuminosityBase``).

A star can hold a luminosity model and derive its luminosity (and effective temperature) from its own
mass and radius. These tests check the attach/transfer-of-ownership, the mass-derived luminosity and
temperature, the update mutator, and the no-model error path.
"""
import math

import numpy as np
import pytest

from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.stellar_x import make_luminosity, MassToLuminosity, FixedLuminosity

MASS_SOLAR = 1.988435e30
RADIUS_SOLAR = 6.957e8
LUM_SOLAR = 3.848e26


def test_no_model_raises():
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR)
    assert star.luminosity_model_set is False
    with pytest.raises(RuntimeError):
        star.calc_luminosity_from_mass()
    with pytest.raises(RuntimeError):
        star.calc_effective_temperature_from_mass()
    with pytest.raises(RuntimeError):
        star.update_luminosity_from_mass()


def test_attach_transfers_ownership():
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR)
    model = MassToLuminosity()
    star.set_luminosity_model(model)
    assert star.luminosity_model_set is True
    # The wrapper is now an empty shell; reusing it must raise.
    with pytest.raises(ValueError):
        star.set_luminosity_model(model)


def test_luminosity_from_mass_matches_model():
    """The star feeds its own mass to the model."""
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR)
    star.set_luminosity_model(MassToLuminosity())
    # A solar-mass star sits in the (M/Msun)^4 branch, so L = Lsun.
    assert math.isclose(star.calc_luminosity_from_mass(), LUM_SOLAR, rel_tol=1e-12)


def test_effective_temperature_from_mass():
    """T_eff derived from mass round-trips through Stefan-Boltzmann on the star's radius."""
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR)
    star.set_luminosity_model(MassToLuminosity())
    lum = star.calc_luminosity_from_mass()
    expected_t = star.calc_temperature_from_luminosity(lum)
    assert math.isclose(star.calc_effective_temperature_from_mass(), expected_t, rel_tol=1e-14)
    # A solar-mass, solar-radius blackbody with L=Lsun is ~5772 K.
    assert math.isclose(star.calc_effective_temperature_from_mass(), 5772.0, rel_tol=2e-3)


def test_update_luminosity_from_mass():
    """The update mutator writes the mass-derived L and T back onto the star's scalar fields."""
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR, effective_temperature_k=3000.0, luminosity_w=1.0)
    star.set_luminosity_model(MassToLuminosity())
    star.update_luminosity_from_mass()
    assert math.isclose(star.luminosity, LUM_SOLAR, rel_tol=1e-12)
    assert math.isclose(star.effective_temperature, star.calc_temperature_from_luminosity(LUM_SOLAR),
                        rel_tol=1e-14)


def test_low_mass_star_uses_low_mass_branch():
    """A TRAPPIST-1-like 0.0898 Msun star lands in the 0.23*(M/Msun)^2.3 branch."""
    ratio = 0.0898
    star = StarWorld("trap1", 0.1192 * RADIUS_SOLAR, ratio * MASS_SOLAR)
    star.set_luminosity_model(MassToLuminosity())
    expected = LUM_SOLAR * 0.23 * ratio ** 2.3
    assert math.isclose(star.calc_luminosity_from_mass(), expected, rel_tol=1e-12)


def test_fixed_model_on_star_ignores_mass():
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR)
    star.set_luminosity_model(make_luminosity("fixed", {"luminosity_w": 2.5e26}))
    assert star.calc_luminosity_from_mass() == 2.5e26


def test_factory_model_attaches():
    star = StarWorld("s", RADIUS_SOLAR, 2.0 * MASS_SOLAR)
    star.set_luminosity_model(make_luminosity("power_law", {"power_law_coeff": 1.4, "power_law_exponent": 3.5}))
    expected = LUM_SOLAR * 1.4 * 2.0 ** 3.5
    assert math.isclose(star.calc_luminosity_from_mass(), expected, rel_tol=1e-12)
