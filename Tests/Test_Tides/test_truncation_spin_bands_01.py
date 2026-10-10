"""Truncation accuracy limits by spin band.

The limits are measured on their own for |spin rate / mean motion| in [0, 1.5], [1.5, 5], and [5, 30], each band
including its edges; a ratio past 30, NaN, or None takes the limits over every measured spin rate.
"""
import math

import pytest

from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.classes.collapse import collapse_global_tides
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Tides.eccentricity import (
    ECCENTRICITY_TRUNCATIONS, eccentricity_accuracy_limit, recommend_eccentricity_truncation)
from TidalPy.Tides.obliquity import obliquity_accuracy_limit, recommend_obliquity_truncation

DEGREES = (2, 3, 4, 5, 6, 7, 8, 9, 10)
TOLERANCES = (1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1)
OBLIQUITY_LEVELS = (2, 4)
# One spin ratio inside each band, and the ratios each band is checked at (both edges and inside).
BAND_RATIOS = (1.0, 3.0, 10.0)
BAND_SPINS = ((0.0, 0.75, 1.5), (1.5, 3.7, 5.0), (5.0, 13.0, 30.0))

_BODY = dict(planet_radius=1.8215e6, orbital_frequency=4.11e-5, semi_major_axis=4.217e8, host_mass=1.898e27,
             G_to_use=6.674e-11)
_TIDES = (("cpl", {"fixed_k": [0.3] * 9, "fixed_q": [100.0] * 9}),
          ("ctl", {"fixed_k": [0.3] * 9, "fixed_dt_s": [100.0] * 9}))
_WARNING_TEXT = "can misstate the tides by 10% or more"


def _single_degree_heating(degree_l, spin_ratio, tide_model, tide_config, eccentricity, eccentricity_truncation):
    """Heating of one degree alone at zero obliquity."""
    return collapse_global_tides(
        **_BODY,
        spin_frequency=spin_ratio * _BODY["orbital_frequency"],
        eccentricity=eccentricity,
        obliquity=0.0,
        tide_model=tide_model,
        tide_config=tide_config,
        min_degree_l=degree_l,
        max_degree_l=degree_l,
        eccentricity_truncation=eccentricity_truncation,
        eccentricity_exact_tolerance=1.0e-10,
        obliquity_truncation=0)["tidal_heating"]


@pytest.mark.parametrize("level, spin_ratio, expected", [
    (50, 1.0, 0.805),
    (50, -1.0, 0.805),
    (50, 3.0, 0.8),
    (50, 10.0, 0.57),
    (50, None, 0.57),
    (50, 100.0, 0.57),
    (50, math.nan, 0.57),
    (10, 1.0, 0.4),
    (10, 10.0, 0.52),
    (10, None, 0.395),
])
def test_eccentricity_limits_are_the_measured_band_values(level, spin_ratio, expected):
    """The degree-2 10% limits of each band; past 30, NaN, and None take the any-spin table."""
    assert eccentricity_accuracy_limit(level, 1.0e-1, 2, spin_ratio=spin_ratio) == expected


@pytest.mark.parametrize("level, spin_ratio, expected", [(2, 1.0, 0.35), (2, 3.0, 0.58), (2, 10.0, 0.565),
                                                         (2, None, 0.35), (4, 3.0, 0.795)])
def test_obliquity_limits_are_the_measured_band_values(level, spin_ratio, expected):
    """Synchronous rotation sets the near-synchronous band's limits, which are those for any spin rate."""
    assert obliquity_accuracy_limit(level, 1.0e-1, 2, spin_ratio=spin_ratio) == expected


@pytest.mark.parametrize("spin_ratio, band_ratio", [
    (1.5, 1.0), (1.5000001, 3.0), (5.0, 3.0), (5.0000001, 10.0), (30.0, 10.0), (30.0000001, None), (0.0, 1.0)])
def test_band_edges(spin_ratio, band_ratio):
    """A band holds its upper edge; a ratio past it takes the next band (level 6 and 10 limits differ in each)."""
    for level in (6, 10):
        assert (eccentricity_accuracy_limit(level, 1.0e-1, 2, spin_ratio=spin_ratio)
                == eccentricity_accuracy_limit(level, 1.0e-1, 2, spin_ratio=band_ratio))


@pytest.mark.parametrize("spin_ratio", BAND_RATIOS)
def test_a_band_never_promises_less_than_any_spin_rate(spin_ratio):
    """Every measured spin rate includes the band's, so a band's limit is never below the any-spin one."""
    for tolerance in TOLERANCES:
        for degree_l in DEGREES:
            for level in ECCENTRICITY_TRUNCATIONS:
                assert (eccentricity_accuracy_limit(level, tolerance, degree_l, spin_ratio=spin_ratio)
                        >= eccentricity_accuracy_limit(level, tolerance, degree_l))
            for level in OBLIQUITY_LEVELS:
                assert (obliquity_accuracy_limit(level, tolerance, degree_l, spin_ratio=spin_ratio)
                        >= obliquity_accuracy_limit(level, tolerance, degree_l))


@pytest.mark.parametrize("band_i", range(len(BAND_RATIOS)))
@pytest.mark.parametrize("level", ECCENTRICITY_TRUNCATIONS)
def test_band_limit_holds_across_the_band(level, band_i):
    """At a level's 10% limit for a band and a solve through degree 10, every degree's constant-phase-lag and
    constant-time-lag heating is within 10% at the band's edges and inside it."""
    eccentricity = eccentricity_accuracy_limit(level, 1.0e-1, 10, spin_ratio=BAND_RATIOS[band_i])
    assert eccentricity > 0.0
    for degree_l in DEGREES:
        for tide_model, tide_config in _TIDES:
            for spin_ratio in BAND_SPINS[band_i]:
                args = (degree_l, spin_ratio, tide_model, tide_config, eccentricity)
                error = abs(_single_degree_heating(*args, level) / _single_degree_heating(*args, "exact") - 1.0)
                assert error <= 1.0e-1, (degree_l, tide_model, spin_ratio)


def test_recommendations_follow_the_spin_band():
    """Near synchronous rotation level 50 holds 10% to e = 0.805; for any spin rate e = 0.7 needs 'exact'."""
    assert recommend_eccentricity_truncation(0.7, tolerance=1.0e-1, spin_ratio=1.0) == 50
    assert recommend_eccentricity_truncation(0.7, tolerance=1.0e-1) == "exact"
    # A fast rotator holds the low levels further: level 8 to e = 0.475, where a synchronous body needs level 20.
    assert recommend_eccentricity_truncation(0.45, tolerance=1.0e-1, spin_ratio=10.0) == 8
    assert recommend_eccentricity_truncation(0.45, tolerance=1.0e-1, spin_ratio=1.0) == 20
    # Past synchronous rotation level 2 holds 10% to 0.565 rad; near it, level 4 is needed past 0.35 rad.
    assert recommend_obliquity_truncation(0.5, tolerance=1.0e-1, spin_ratio=10.0) == 2
    assert recommend_obliquity_truncation(0.5, tolerance=1.0e-1) == 4


@pytest.mark.parametrize("spin_ratio, warns", [(1.0, True), (10.0, False)])
def test_a_world_warns_by_its_spin_band(spdlog_text, spin_ratio, warns):
    """At e = 0.45, level 10 is past its near-synchronous range (0.4) but inside its fast-rotator range (0.52)."""
    world = StarWorld(f"spin_{spin_ratio:g}", 7.0e7, 1.9e27)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [1.0e4]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=10, obliquity_truncation=0)
    world.calc_tides(orbital_frequency=2.0e-5, spin_frequency=spin_ratio * 2.0e-5, eccentricity=0.45, obliquity=0.0,
                     semi_major_axis=1.0e10, host_mass=2.0e30)
    assert (_WARNING_TEXT in spdlog_text()) == warns
