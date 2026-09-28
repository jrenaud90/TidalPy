"""Tests for the dynamically loaded third-party constants in TidalPy.constants."""
import math

import pytest

from TidalPy.constants import (
    G, au, k, k_boltzman, k_boltzmann, luminosity_solar, radius_jupiter, seconds_per_myr, year, yr)


@pytest.mark.parametrize(
    "value, expected, rel_tol",
    [
        (year, 31_557_600.0, 1e-12),
        (G, 6.6743e-11, 1e-3),
        (au, 1.495978707e11, 1e-6),
        (luminosity_solar, 3.828e26, 1e-12),
        (radius_jupiter, 6.9911e7, 1e-12),
        (k_boltzmann, 1.380649e-23, 1e-9),
    ],
    ids=["year_julian", "G", "au", "luminosity_solar_iau2015", "radius_jupiter_iau2015", "k_boltzmann"])
def test_constant_value(value, expected, rel_tol):
    """Each SciPy-sourced or reference-body constant is populated with its expected value."""
    assert math.isclose(value, expected, rel_tol=rel_tol)


def test_constant_aliases():
    """Aliases (yr, k, k_boltzman) equal their primary constants."""
    assert yr == year
    assert k == k_boltzmann
    assert k_boltzman == k_boltzmann


def test_seconds_per_myr_is_exact_julian():
    """The compile-time mega-year equals one million of SciPy's Julian years, exactly."""
    assert seconds_per_myr == 3.15576e13
    assert math.isclose(seconds_per_myr, 1.0e6 * year, rel_tol=1e-15)
