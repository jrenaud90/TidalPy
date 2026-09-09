"""Tests for ``TidalPy.Utilities_x.conversions`` (unit and orbital-element conversions)."""

from math import isclose, pi

import pytest

from TidalPy.constants import G, au
from TidalPy.exceptions import BadValueError
from TidalPy.Utilities_x.conversions import (
    Au2m,
    days2rads,
    m2Au,
    myr2sec,
    orbital_motion2semi_a,
    rads2days,
    sec2myr,
    semi_a2orbital_motion,
)


def test_distance_round_trip():
    assert isclose(Au2m(1.0), 149597870700.0)
    assert isclose(m2Au(Au2m(2.5)), 2.5)
    assert isclose(Au2m(1.0), au, rel_tol=1e-9)


def test_frequency_period_round_trip():
    one_day_freq = 2.0 * pi / 86400.0
    assert isclose(rads2days(one_day_freq), 1.0)
    assert isclose(days2rads(rads2days(0.5 * one_day_freq)), 0.5 * one_day_freq)


def test_time_round_trip():
    assert isclose(sec2myr(myr2sec(123.0)), 123.0)


def test_kepler_round_trip():
    host_mass = 1.989e30
    target_mass = 5.972e24
    semi_major_axis = 1.4959787e11
    orbital_motion = semi_a2orbital_motion(semi_major_axis, host_mass, target_mass)
    expected = (G * (host_mass + target_mass) / semi_major_axis ** 3) ** 0.5
    assert isclose(orbital_motion, expected, rel_tol=1e-12)
    assert isclose(orbital_motion2semi_a(orbital_motion, host_mass, target_mass),
                   semi_major_axis, rel_tol=1e-12)


def test_kepler_matches_classic_implementation():
    """The ported conversions equal the classic numpy implementation while both backends exist."""
    from TidalPy.utilities.conversions import (
        orbital_motion2semi_a as classic_n2a,
        semi_a2orbital_motion as classic_a2n,
    )
    host_mass = 6.4e23
    orbital_motion = 2.0 * pi / (2.0 * 86400.0)
    assert isclose(orbital_motion2semi_a(orbital_motion, host_mass),
                   float(classic_n2a(orbital_motion, host_mass)), rel_tol=1e-12)
    semi_major_axis = 4.0e8
    assert isclose(semi_a2orbital_motion(semi_major_axis, host_mass),
                   float(classic_a2n(semi_major_axis, host_mass)), rel_tol=1e-12)


def test_bad_masses_raise():
    with pytest.raises(BadValueError):
        orbital_motion2semi_a(1.0e-5, 0.0)
    with pytest.raises(BadValueError):
        semi_a2orbital_motion(1.0e8, 1.0e24, target_mass=-1.0)
