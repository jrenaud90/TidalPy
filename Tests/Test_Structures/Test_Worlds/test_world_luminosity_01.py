"""StarWorld luminosity models: ownership transfer, mass-derived luminosity and temperature, and the no-model error."""
import math

import pytest

from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Stellar import make_luminosity, MassToLuminosity

MASS_SOLAR = 1.988435e30
RADIUS_SOLAR = 6.957e8
LUM_SOLAR = 3.828e26


def test_no_model_raises():
    """Mass-derived calls raise RuntimeError when no luminosity model is set."""
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR)
    assert star.luminosity_model_set is False
    with pytest.raises(RuntimeError):
        star.calc_luminosity_from_mass()
    with pytest.raises(RuntimeError):
        star.calc_effective_temperature_from_mass()
    with pytest.raises(RuntimeError):
        star.update_luminosity_from_mass()


def test_attach_transfers_ownership():
    """Attaching a model consumes its wrapper, so reattaching it raises."""
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR)
    model = MassToLuminosity()
    star.set_luminosity_model(model)
    assert star.luminosity_model_set is True
    with pytest.raises(ValueError):
        star.set_luminosity_model(model)


@pytest.mark.parametrize(
    "radius, mass, make_model, expected, rel_tol",
    [
        # A solar-mass star sits in the (M/Msun)^4 branch, so L = Lsun.
        pytest.param(RADIUS_SOLAR, MASS_SOLAR, MassToLuminosity, LUM_SOLAR, 1e-12, id="solar_mass"),
        # A TRAPPIST-1-like star lands in the 0.23 (M/Msun)^2.3 branch.
        pytest.param(
            0.1192 * RADIUS_SOLAR,
            0.0898 * MASS_SOLAR,
            MassToLuminosity,
            LUM_SOLAR * 0.23 * 0.0898 ** 2.3,
            1e-12,
            id="low_mass",
        ),
        # A fixed model ignores the mass; rel_tol 0 demands exact equality.
        pytest.param(
            RADIUS_SOLAR,
            MASS_SOLAR,
            lambda: make_luminosity("fixed", {"luminosity_w": 2.5e26}),
            2.5e26,
            0.0,
            id="fixed",
        ),
        pytest.param(
            RADIUS_SOLAR,
            2.0 * MASS_SOLAR,
            lambda: make_luminosity("power_law", {"power_law_coeff": 1.4, "power_law_exponent": 3.5}),
            LUM_SOLAR * 1.4 * 2.0 ** 3.5,
            1e-12,
            id="power_law_factory",
        ),
    ],
)
def test_luminosity_from_mass(
        radius,
        mass,
        make_model,
        expected,
        rel_tol,
):
    """The star feeds its own mass to the attached luminosity model."""
    star = StarWorld("s", radius, mass)
    star.set_luminosity_model(make_model())
    assert math.isclose(star.calc_luminosity_from_mass(), expected, rel_tol=rel_tol)


def test_effective_temperature_from_mass():
    """T_eff from mass round-trips through Stefan-Boltzmann on the star's radius."""
    star = StarWorld("s", RADIUS_SOLAR, MASS_SOLAR)
    star.set_luminosity_model(MassToLuminosity())
    lum = star.calc_luminosity_from_mass()
    expected_t = star.calc_temperature_from_luminosity(lum)
    assert math.isclose(star.calc_effective_temperature_from_mass(), expected_t, rel_tol=1e-14)
    assert math.isclose(star.calc_effective_temperature_from_mass(), 5772.0, rel_tol=1e-3)


def test_update_luminosity_from_mass():
    """update_luminosity_from_mass writes the mass-derived L and T onto the star."""
    star = StarWorld(
        "s",
        RADIUS_SOLAR,
        MASS_SOLAR,
        effective_temperature=3000.0,
        luminosity=1.0,
    )
    star.set_luminosity_model(MassToLuminosity())
    star.update_luminosity_from_mass()
    assert math.isclose(star.luminosity, LUM_SOLAR, rel_tol=1e-12)
    assert math.isclose(star.effective_temperature, star.calc_temperature_from_luminosity(LUM_SOLAR),
                        rel_tol=1e-14)
