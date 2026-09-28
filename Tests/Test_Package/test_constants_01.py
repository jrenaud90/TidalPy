"""Tests for constants added or corrected in TidalPy.constants."""
import math

from TidalPy.constants import G, au, k_boltzmann, luminosity_solar, mass_jupiter, radius_jupiter, year, yr


def test_year_is_julian_year():
    """The year constant is the Julian year in seconds, with a matching yr alias."""
    assert math.isclose(year, 31_557_600.0, rel_tol=1e-12)
    assert yr == year


def test_third_party_constants_populated():
    """The SciPy-sourced constants are populated."""
    assert math.isclose(G, 6.6743e-11, rel_tol=1e-3)
    assert math.isclose(au, 1.495978707e11, rel_tol=1e-6)
    assert math.isclose(k_boltzmann, 1.380649e-23, rel_tol=1e-9)


def test_reference_body_constants():
    """The solar luminosity and Jupiter radius are IAU nominal values; Jupiter's mass is the NASA fact sheet value."""
    assert math.isclose(luminosity_solar, 3.828e26, rel_tol=1e-12)
    assert math.isclose(mass_jupiter, 1.898125e27, rel_tol=1e-12)
    assert math.isclose(radius_jupiter, 6.9911e7, rel_tol=1e-12)
