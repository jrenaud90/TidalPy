"""Truncation accuracy limits by degree: a solve through degree L takes the tightest limit of degrees 2 to L."""
import numpy as np
import pytest

from TidalPy.Tides.classes.collapse import collapse_global_tides
from TidalPy.Tides.eccentricity import (
    ECCENTRICITY_TRUNCATIONS, eccentricity_accuracy_limit, recommend_eccentricity_truncation)
from TidalPy.Tides.obliquity import obliquity_accuracy_limit, recommend_obliquity_truncation

DEGREES = (2, 3, 4, 5, 6, 7, 8, 9, 10)
TOLERANCES = (1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1)
OBLIQUITY_LEVELS = (2, 4)
SPIN_RATIOS = (0.5, 1.0, 2.3)

_BODY = dict(planet_radius=1.8215e6, orbital_frequency=4.11e-5, semi_major_axis=4.217e8, host_mass=1.898e27,
             G_to_use=6.674e-11)
_TIDES = (("cpl", {"fixed_k": [0.3] * 9, "fixed_q": [100.0] * 9}),
          ("ctl", {"fixed_k": [0.3] * 9, "fixed_dt_s": [100.0] * 9}))


def _single_degree_heating(
        degree_l,
        spin_ratio,
        tide_model,
        tide_config,
        eccentricity,
        obliquity,
        eccentricity_truncation,
        obliquity_truncation):
    """Heating of one degree alone."""
    return collapse_global_tides(
        **_BODY,
        spin_frequency=spin_ratio * _BODY["orbital_frequency"],
        eccentricity=eccentricity,
        obliquity=obliquity,
        tide_model=tide_model,
        tide_config=tide_config,
        min_degree_l=degree_l,
        max_degree_l=degree_l,
        eccentricity_truncation=eccentricity_truncation,
        eccentricity_exact_tolerance=1.0e-10,
        obliquity_truncation=obliquity_truncation)["tidal_heating"]


def _worst_eccentricity_error(level, degree_l, eccentricity):
    """Largest relative heating error of one degree at a level against the exact functions."""
    worst = 0.0
    for tide_model, tide_config in _TIDES:
        for spin_ratio in SPIN_RATIOS:
            args = (degree_l, spin_ratio, tide_model, tide_config, eccentricity, 0.0)
            truncated = _single_degree_heating(*args, level, 0)
            exact = _single_degree_heating(*args, "exact", 0)
            worst = max(worst, abs(truncated / exact - 1.0))
    return worst


def _worst_obliquity_error(level, degree_l, obliquity):
    """Largest relative heating error of one degree at a level against the general functions, at e = 0."""
    worst = 0.0
    for tide_model, tide_config in _TIDES:
        for spin_ratio in SPIN_RATIOS:
            args = (degree_l, spin_ratio, tide_model, tide_config, 0.0, obliquity, 2)
            worst = max(worst, abs(_single_degree_heating(*args, level) / _single_degree_heating(*args, "gen") - 1.0))
    return worst


@pytest.mark.parametrize("level, degree_two, degree_three", [
    (2, 0.075, 0.06),
    (10, 0.405, 0.36),
    (20, 0.595, 0.535),
    (50, 0.78, 0.75),
])
def test_eccentricity_limits_through_degree_three_are_unchanged(level, degree_two, degree_three):
    """The 10% limits of a solve through degree 2 or 3 are the measured degree-2 and degree-3 values."""
    assert eccentricity_accuracy_limit(level, 1.0e-1, 2) == degree_two
    assert eccentricity_accuracy_limit(level, 1.0e-1, 3) == degree_three


@pytest.mark.parametrize("level, degree_two, degree_three", [(2, 0.465, 0.31), (4, 0.83, 0.565)])
def test_obliquity_limits_through_degree_three_are_unchanged(level, degree_two, degree_three):
    """The 10% limits of a solve through degree 2 or 3 are the measured degree-2 and degree-3 values."""
    assert obliquity_accuracy_limit(level, 1.0e-1, 2) == degree_two
    assert obliquity_accuracy_limit(level, 1.0e-1, 3) == degree_three


@pytest.mark.parametrize("tolerance", TOLERANCES)
@pytest.mark.parametrize("level", ECCENTRICITY_TRUNCATIONS)
def test_eccentricity_limits_never_grow_with_degree(level, tolerance):
    """Including a higher degree never raises the limit; a degree past the table uses degree 10's."""
    limits = [eccentricity_accuracy_limit(level, tolerance, degree_l) for degree_l in DEGREES]
    assert all(higher <= lower for lower, higher in zip(limits, limits[1:]))
    assert eccentricity_accuracy_limit(level, tolerance, 11) == limits[-1]


@pytest.mark.parametrize("tolerance", TOLERANCES)
@pytest.mark.parametrize("level", OBLIQUITY_LEVELS)
def test_obliquity_limits_never_grow_with_degree(level, tolerance):
    """Including a higher degree never raises the limit; a degree past the table uses degree 10's."""
    limits = [obliquity_accuracy_limit(level, tolerance, degree_l) for degree_l in DEGREES]
    assert all(higher <= lower for lower, higher in zip(limits, limits[1:]))
    assert obliquity_accuracy_limit(level, tolerance, 11) == limits[-1]


def test_high_degrees_tighten_the_limits():
    """Degrees 4 to 10 lose accuracy sooner than degree 3, so a solve through degree 10 has tighter limits."""
    for level in ECCENTRICITY_TRUNCATIONS:
        assert eccentricity_accuracy_limit(level, 1.0e-1, 10) < eccentricity_accuracy_limit(level, 1.0e-1, 3)
    for level in OBLIQUITY_LEVELS:
        assert obliquity_accuracy_limit(level, 1.0e-1, 10) < obliquity_accuracy_limit(level, 1.0e-1, 3)


@pytest.mark.parametrize("tolerance", (1.0e-2, 1.0e-1))
@pytest.mark.parametrize("level", ECCENTRICITY_TRUNCATIONS)
def test_eccentricity_limit_through_degree_ten_holds_at_every_degree(level, tolerance):
    """At a level's limit for a solve through degree 10, every degree's heating is within the tolerance."""
    eccentricity = eccentricity_accuracy_limit(level, tolerance, 10)
    assert eccentricity > 0.0
    for degree_l in DEGREES:
        assert _worst_eccentricity_error(level, degree_l, eccentricity) <= tolerance, degree_l


@pytest.mark.parametrize("tolerance", (1.0e-2, 1.0e-1))
@pytest.mark.parametrize("level", OBLIQUITY_LEVELS)
def test_obliquity_limit_through_degree_ten_holds_at_every_degree(level, tolerance):
    """At a level's limit for a solve through degree 10, every degree's heating is within the tolerance."""
    obliquity = obliquity_accuracy_limit(level, tolerance, 10)
    assert obliquity > 0.0
    for degree_l in DEGREES:
        assert _worst_obliquity_error(level, degree_l, obliquity) <= tolerance, degree_l


@pytest.mark.parametrize("eccentricity", (0.05, 0.2, 0.36))
def test_recommended_eccentricity_level_holds_at_every_degree(eccentricity):
    """The level recommended for a solve through degree 10 holds 10% at every degree (level 10 at e = 0.36 did not)."""
    level = recommend_eccentricity_truncation(eccentricity, tolerance=1.0e-1, max_degree_l=10)
    for degree_l in DEGREES:
        assert _worst_eccentricity_error(level, degree_l, eccentricity) <= 1.0e-1, degree_l


@pytest.mark.parametrize("obliquity", (0.05, 0.15))
def test_recommended_obliquity_level_holds_at_every_degree(obliquity):
    """The level recommended for a solve through degree 10 holds 10% at every degree."""
    level = recommend_obliquity_truncation(obliquity, tolerance=1.0e-1, max_degree_l=10)
    assert level in OBLIQUITY_LEVELS
    for degree_l in DEGREES:
        assert _worst_obliquity_error(level, degree_l, obliquity) <= 1.0e-1, degree_l
