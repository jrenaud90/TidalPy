"""Tests for ``TidalPy.Utilities.conversions`` (unit and orbital element conversions).

Classic Kepler results were frozen from the 0.8.0 snapshot 8b8e0b12 into ``frozen/test_conversions_01.npz``.
"""
from math import isclose, pi
from pathlib import Path

import numpy as np
import pytest

from TidalPy.constants import G, au, seconds_per_myr, year
from TidalPy.Utilities.conversions import (
    Au2m,
    days2rads,
    m2Au,
    myr2sec,
    orbital_motion2semi_a,
    rads2days,
    sec2myr,
    semi_a2orbital_motion,
)

FROZEN_PATH = Path(__file__).parent / "frozen" / "test_conversions_01.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_REFERENCE = {key: float(frozen_file[key]) for key in frozen_file.files}


def test_distance_round_trip():
    """Au and meter conversions agree with the au constant and round trip."""
    assert isclose(Au2m(1.0), 149597870700.0)
    assert isclose(m2Au(Au2m(2.5)), 2.5)
    assert isclose(Au2m(1.0), au, rel_tol=1e-9)


def test_frequency_period_round_trip():
    """Frequency and period (days) conversions round trip."""
    one_day_freq = 2.0 * pi / 86400.0
    assert isclose(rads2days(one_day_freq), 1.0)
    assert isclose(days2rads(rads2days(0.5 * one_day_freq)), 0.5 * one_day_freq)


def test_time_round_trip():
    """Seconds and mega-year conversions round trip."""
    assert isclose(sec2myr(myr2sec(123.0)), 123.0)


def test_myr_is_julian():
    """A mega-year is one million Julian years."""
    assert isclose(myr2sec(1.0), 1.0e6 * year, rel_tol=1e-15)
    assert myr2sec(1.0) == seconds_per_myr
    assert isclose(sec2myr(3.15576e13), 1.0, rel_tol=1e-15)


def test_kepler_round_trip():
    """Kepler's third law conversions match the analytic value and round trip."""
    host_mass = 1.989e30
    target_mass = 5.972e24
    semi_major_axis = 1.4959787e11
    orbital_motion = semi_a2orbital_motion(semi_major_axis, host_mass, target_mass)
    expected = (G * (host_mass + target_mass) / semi_major_axis ** 3) ** 0.5
    assert isclose(orbital_motion, expected, rel_tol=1e-12)
    assert isclose(orbital_motion2semi_a(orbital_motion, host_mass, target_mass),
                   semi_major_axis, rel_tol=1e-12)


def test_kepler_matches_classic_implementation():
    """Kepler conversions equal the frozen classic results."""
    host_mass = 6.4e23
    orbital_motion = 2.0 * pi / (2.0 * 86400.0)
    assert isclose(orbital_motion2semi_a(orbital_motion, host_mass),
                   CLASSIC_REFERENCE["classic_orbital_motion2semi_a"], rel_tol=1e-12)
    semi_major_axis = 4.0e8
    assert isclose(semi_a2orbital_motion(semi_major_axis, host_mass),
                   CLASSIC_REFERENCE["classic_semi_a2orbital_motion"], rel_tol=1e-12)


@pytest.mark.parametrize(
    "convert",
    [lambda: orbital_motion2semi_a(1.0e-5, 0.0), lambda: semi_a2orbital_motion(1.0e8, 1.0e24, target_mass=-1.0)],
    ids=["zero_host_mass", "negative_target_mass"])
def test_bad_masses_raise(convert):
    """Invalid masses raise ValueError."""
    with pytest.raises(ValueError):
        convert()
