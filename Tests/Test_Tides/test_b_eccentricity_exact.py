"""The exact eccentricity functions (Hansen coefficients from the exact Kepler orbit) and the truncation helper.

``eccentricity_trunc_lvl = "exact"`` takes G_lpq from a transform of the exact orbit, keeping the modes whose q^2-weighted
squares leave a tail below the exact tolerance; ``recommend_eccentricity_truncation`` picks the lowest tabulated level
that holds a heating tolerance at an eccentricity (see Documentation/Tides_x/eccentricity.md).
"""
from math import isclose

import numpy as np
import pytest

from TidalPy.Tides_x.eccentricity import (
    ECCENTRICITY_EXACT, ECCENTRICITY_TRUNCATIONS, eccentricity_accuracy_limit, eccentricity_func,
    eccentricity_squared_func, recommend_eccentricity_truncation, validate_eccentricity_exact_tolerance,
    validate_eccentricity_truncation)


def _hut_synchronous(eccentricity):
    """Hut (1981) constant-time-lag heating at synchronous rotation in units of the l = 2 mode sum below."""
    e2 = eccentricity * eccentricity
    n_a = (1 + 31 / 2 * e2 + 255 / 8 * e2**2 + 185 / 16 * e2**3 + 25 / 64 * e2**4) / (1 - e2)**7.5
    n_e = (1 + 15 / 2 * e2 + 45 / 8 * e2**2 + 5 / 16 * e2**3) / (1 - e2)**6
    omega_e = (1 + 3 * e2 + 3 / 8 * e2**2) / (1 - e2)**4.5
    return 3.0 * (n_a - 2.0 * n_e + omega_e)


def _synchronous_mode_sum(by_lpq):
    weight = {0: 0.75, 1: 0.25}
    return sum(weight[p] * value * q**2 for (l, p, q), value in by_lpq if p in weight)


def test_exact_is_a_truncation():
    assert validate_eccentricity_truncation("exact") == ECCENTRICITY_EXACT
    assert validate_eccentricity_truncation("EXACT") == ECCENTRICITY_EXACT
    assert ECCENTRICITY_EXACT not in ECCENTRICITY_TRUNCATIONS
    assert validate_eccentricity_exact_tolerance(1.0e-6) == 1.0e-6
    for bad in (0.0, 1.0, -1.0e-3):
        with pytest.raises(ValueError):
            validate_eccentricity_exact_tolerance(bad)
    with pytest.raises(ValueError):
        eccentricity_func(0.3, 2, "exact", exact_tolerance=2.0)


@pytest.mark.parametrize("degree_l", (2, 3, 4, 7, 10))
def test_exact_matches_the_tables_at_moderate_eccentricity(degree_l):
    """At e = 0.2 level 50 is exact to rounding, so the two agree on every mode both hold."""
    exact = dict(eccentricity_func(0.2, degree_l, "exact", exact_tolerance=1.0e-12)[0])
    table = dict(eccentricity_func(0.2, degree_l, 50)[0])
    shared = set(exact) & set(table)
    assert len(shared) > 0
    for key in shared:
        assert abs(exact[key] - table[key]) <= 1.0e-11 * max(1.0, abs(table[key])), key


def test_exact_k0_mode_is_the_closed_form():
    e = 0.7
    by_lpq = dict(eccentricity_func(e, 2, "exact")[0])
    assert isclose(by_lpq[(2, 1, 0)], (1.0 - e**2)**-1.5, rel_tol=1.0e-12)
    # The k = 0 modes with |m| >= l vanish identically and are not returned.
    assert (2, 0, -2) not in by_lpq and (2, 2, 2) not in by_lpq


@pytest.mark.parametrize("eccentricity", (0.1, 0.5, 0.7, 0.8, 0.9))
def test_exact_heating_matches_hut(eccentricity):
    """The tolerance bounds the relative error of the synchronous constant-time-lag heating."""
    tolerance = 1.0e-4
    squares = eccentricity_squared_func(eccentricity, 2, "exact", exact_tolerance=tolerance)[0]
    error = _synchronous_mode_sum(squares) / _hut_synchronous(eccentricity) - 1.0
    assert -tolerance <= error <= 1.0e-10


def test_exact_mode_range_grows_with_eccentricity_and_tolerance():
    def max_q(e, tolerance):
        return max(abs(q) for (l, p, q), _ in eccentricity_func(e, 2, "exact", exact_tolerance=tolerance)[0])
    assert max_q(0.3, 1.0e-4) < max_q(0.6, 1.0e-4) < max_q(0.9, 1.0e-4)
    assert max_q(0.6, 1.0e-2) < max_q(0.6, 1.0e-6)
    # A circular orbit keeps only q = 0.
    assert max_q(0.0, 1.0e-4) == 0


def test_recommend_eccentricity_truncation():
    assert recommend_eccentricity_truncation(0.0) == 2
    assert recommend_eccentricity_truncation(0.05) == 4
    assert recommend_eccentricity_truncation(0.3, tolerance=1.0e-2) == 10
    assert recommend_eccentricity_truncation(0.3, tolerance=1.0e-2, max_degree_l=3) == 20
    assert recommend_eccentricity_truncation(0.55, tolerance=1.0e-6) == 50
    assert recommend_eccentricity_truncation(0.6, tolerance=1.0e-6) == "exact"
    assert recommend_eccentricity_truncation(0.8) == "exact"
    assert recommend_eccentricity_truncation(0.3, tolerance=1.0e-12) == "exact"
    with pytest.raises(ValueError):
        recommend_eccentricity_truncation(1.0)


@pytest.mark.parametrize("level", ECCENTRICITY_TRUNCATIONS)
def test_recommended_levels_hold_their_tolerance(level):
    """At each level's own 1% limit the synchronous constant-time-lag heating is within 1%."""
    limit = eccentricity_accuracy_limit(level, 1.0e-2)
    assert 0.0 < limit < 0.8
    error = _synchronous_mode_sum(eccentricity_squared_func(limit, 2, level)[0]) / _hut_synchronous(limit) - 1.0
    assert abs(error) < 1.0e-2
    assert np.isnan(eccentricity_accuracy_limit("exact", 1.0e-2))
