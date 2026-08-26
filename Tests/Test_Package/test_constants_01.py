"""Tests for the dynamically loaded third-party constants in TidalPy.constants."""
import math

from TidalPy.constants import G, au, year, yr


def test_year_is_julian_year():
    """The year constant is the Julian year in seconds, with a matching yr alias."""
    assert math.isclose(year, 31_557_600.0, rel_tol=1e-12)
    assert yr == year


def test_third_party_constants_populated():
    """The SciPy-sourced constants are populated (not NaN) after initialization."""
    assert math.isclose(G, 6.6743e-11, rel_tol=1e-3)
    assert math.isclose(au, 1.495978707e11, rel_tol=1e-6)
