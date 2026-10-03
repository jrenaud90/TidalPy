"""Eccentricity functions G_lpq(e): the tabulated levels, the exact functions, and the truncation helpers."""
from math import isclose

import numpy as np
import pytest

from TidalPy.Tides.eccentricity import (
    ECCENTRICITY_EXACT, ECCENTRICITY_TRUNCATIONS, eccentricity_accuracy_limit, eccentricity_func,
    eccentricity_squared_func, recommend_eccentricity_truncation, validate_eccentricity_exact_tolerance,
    validate_eccentricity_truncation)
from TidalPy.Utilities.lookups import IntMap3

DEGREES = (2, 3, 4, 5, 6, 7, 8, 9, 10)


def _hut_synchronous(eccentricity):
    """Hut (1981) synchronous constant-time-lag heating in the units of `_synchronous_mode_sum`."""
    e2 = eccentricity * eccentricity
    n_a = (1 + 31 / 2 * e2 + 255 / 8 * e2**2 + 185 / 16 * e2**3 + 25 / 64 * e2**4) / (1 - e2)**7.5
    n_e = (1 + 15 / 2 * e2 + 45 / 8 * e2**2 + 5 / 16 * e2**3) / (1 - e2)**6
    omega_e = (1 + 3 * e2 + 3 / 8 * e2**2) / (1 - e2)**4.5
    return 3.0 * (n_a - 2.0 * n_e + omega_e)


def _synchronous_mode_sum(squares_by_lpq):
    """Degree-2, zero-obliquity, synchronous constant-time-lag heating from G^2 (each mode's frequency is q n)."""
    weight = {0: 0.75, 1: 0.25}
    return sum(weight[p] * value * q**2 for (l, p, q), value in squares_by_lpq if p in weight)


def _vanishing(degree_l, max_q):
    """Count of k = 0 modes with |l - 2p| >= l, which vanish identically."""
    return sum(1 for p in (0, degree_l) if degree_l <= max_q)


def test_tabulated_levels():
    """The tabulated levels, and rejection of untabulated ones."""
    assert ECCENTRICITY_TRUNCATIONS == (2, 4, 6, 8, 10, 20, 50)
    assert validate_eccentricity_truncation("20") == 20
    for bad in (0, 1, 3, 7, 60):
        with pytest.raises(NotImplementedError):
            validate_eccentricity_truncation(bad)
    with pytest.raises(TypeError):
        validate_eccentricity_truncation(True)


def test_exact_is_a_truncation():
    """The "exact" truncation is accepted and its tolerance is validated."""
    assert validate_eccentricity_truncation("exact") == ECCENTRICITY_EXACT
    assert validate_eccentricity_truncation("EXACT") == ECCENTRICITY_EXACT
    assert ECCENTRICITY_EXACT not in ECCENTRICITY_TRUNCATIONS
    assert validate_eccentricity_exact_tolerance(1.0e-6) == 1.0e-6
    for bad in (0.0, 1.0, -1.0e-3):
        with pytest.raises(ValueError):
            validate_eccentricity_exact_tolerance(bad)
    with pytest.raises(ValueError):
        eccentricity_func(0.3, 2, "exact", exact_tolerance=2.0)


@pytest.mark.parametrize('degree_l', DEGREES)
def test_untabulated_level_raises(degree_l):
    """Level 0 is not tabulated."""
    with pytest.raises(NotImplementedError):
        eccentricity_func(0.5, degree_l, 0)


@pytest.mark.parametrize('degree_l', DEGREES)
@pytest.mark.parametrize('truncation', (2, 4, 6, 10))
def test_eccentricity_funcs(degree_l, truncation):
    """Both return maps hold int keys within the degree and level, and float values."""
    by_lpq, by_lp = eccentricity_func(0.5, degree_l, truncation)
    assert isinstance(by_lpq, IntMap3)
    assert isinstance(by_lp, dict)

    assert len(by_lpq) > 0
    for (l, p, q), value in by_lpq:
        assert isinstance(l, int)
        assert isinstance(p, int)
        assert isinstance(q, int)
        assert l == degree_l
        assert p <= l
        assert abs(q) <= truncation
        assert isinstance(value, float)

    for (l, p), by_q in by_lp.items():
        assert isinstance(l, int)
        assert isinstance(p, int)
        assert l == degree_l
        assert p <= l
        for (q,), value in by_q:
            assert isinstance(q, int)
            assert isinstance(value, float)


def test_degree_two_level_two_values():
    """Degree-2 level-2 counts and values (Kaula 1964; Cayley 1861)."""
    e = 0.5
    by_lpq, by_lp = eccentricity_func(e, 2, 2)
    # 3 (2 * 2 + 1) modes less G_20,-2 and G_22,2, which vanish identically.
    assert len(by_lpq) == 13
    assert isclose(by_lpq[(2, 0, 0)], 1.0 - 2.5 * e**2)
    assert isclose(by_lpq[(2, 0, -1)], -0.5 * e)
    assert isclose(by_lpq[(2, 0, 1)], 3.5 * e)
    assert isclose(by_lpq[(2, 0, 2)], 8.5 * e**2)
    assert isclose(by_lpq[(2, 1, 0)], (1.0 - e**2)**-1.5)
    assert len(by_lp) == 3
    assert len(by_lp[(2, 0)]) == 4
    assert len(by_lp[(2, 1)]) == 5
    assert len(by_lp[(2, 2)]) == 4
    assert isclose(by_lp[(2, 1)][(1,)], 1.5 * e)
    assert isclose(by_lp[(2, 1)][(-2,)], 2.25 * e**2)


def test_degree_two_level_four_values():
    """Degree-2 level-4 values through e^4."""
    e = 0.5
    by_lpq = eccentricity_func(e, 2, 4)[0]
    assert isclose(by_lpq[(2, 0, 0)], 1.0 - 2.5 * e**2 + 0.8125 * e**4)
    assert isclose(by_lpq[(2, 0, 2)], 8.5 * e**2 - (115.0 / 6.0) * e**4)
    assert isclose(by_lpq[(2, 0, 1)], 3.5 * e - 123.0 / 16.0 * e**3)
    assert isclose(by_lpq[(2, 1, 1)], 1.5 * e + 27.0 / 16.0 * e**3)


@pytest.mark.parametrize('degree_l', DEGREES)
@pytest.mark.parametrize('truncation', ECCENTRICITY_TRUNCATIONS)
def test_mode_count(degree_l, truncation):
    """At most (l + 1)(2N + 1) modes and (l + 1)(N + 1) squares; exact counts at degree 2."""
    by_lpq = dict(eccentricity_func(0.3, degree_l, truncation)[0])
    by_lpq_squared = dict(eccentricity_squared_func(0.3, degree_l, truncation)[0])
    assert len(by_lpq) <= (degree_l + 1) * (2 * truncation + 1)
    assert set(by_lpq_squared) <= set(by_lpq)
    assert all(2 * abs(q) <= truncation for (l, p, q) in by_lpq_squared)
    # Past degree 2 some modes with a zero leading coefficient drop out at low levels, so only degree 2 is exact.
    if degree_l == 2:
        assert len(by_lpq) == 3 * (2 * truncation + 1) - _vanishing(2, truncation)
        assert len(by_lpq_squared) == 3 * (truncation + 1) - _vanishing(2, truncation // 2)


def test_squared_spot_checks():
    """Level-2 squares are the leading-order squares; the k = 0 modes are exact."""
    e = 0.3
    by_lpq, _ = eccentricity_squared_func(e, 2, 2)
    assert isclose(by_lpq[(2, 0, 0)], 1.0 - 5.0 * e**2)
    assert isclose(by_lpq[(2, 0, 1)], 12.25 * e**2)
    assert isclose(by_lpq[(2, 0, -1)], 0.25 * e**2)
    assert isclose(by_lpq[(2, 1, 1)], 2.25 * e**2)
    assert isclose(by_lpq[(2, 1, 0)], (1.0 - e**2)**-3)
    assert (2, 0, 2) not in dict(by_lpq)   # its square starts at e^4


@pytest.mark.parametrize('degree_l', DEGREES)
@pytest.mark.parametrize('truncation', (2, 4))
def test_eccentricity_truncation_order(degree_l, truncation):
    """Level N carries G and its cut square through e^N: halving e shrinks the gap to level 50 by at least 2^(N+1)."""
    # Higher levels are not checked: their gaps at small e are below double precision.
    large_eccentricity = 0.01
    small_eccentricity = 0.005
    for function in (eccentricity_func, eccentricity_squared_func):
        gaps = dict()
        for eccentricity in (large_eccentricity, small_eccentricity):
            reference = dict(function(eccentricity, degree_l, 50)[0])
            truncated = dict(function(eccentricity, degree_l, truncation)[0])
            assert set(truncated) <= set(reference)
            gaps[eccentricity] = {key: abs(truncated.get(key, 0.0) - value) for key, value in reference.items()}
        num_checked = 0
        for key, gap_small in gaps[small_eccentricity].items():
            # Gaps at round-off carry no ratio information.
            if gap_small < 1.0e-15:
                continue
            num_checked += 1
            assert gaps[large_eccentricity][key] / gap_small > 0.9 * 2**(truncation + 1), (function.__name__, key)
        assert num_checked > 0


@pytest.mark.parametrize('eccentricity', (0.2, 0.4, 0.5, 0.6))
def test_levels_converge_to_hut(eccentricity):
    """Against Hut (1981) every level errs low, a higher level is never worse, and level 50 is within 1e-5."""
    exact = _hut_synchronous(eccentricity)
    errors = [_synchronous_mode_sum(eccentricity_squared_func(eccentricity, 2, level)[0]) / exact - 1.0
              for level in ECCENTRICITY_TRUNCATIONS]
    for error in errors:
        assert error <= 1.0e-12
    for lower, higher in zip(errors, errors[1:]):
        assert abs(higher) <= max(abs(lower), 1.0e-13)
    assert abs(errors[-1]) < 1.0e-5


@pytest.mark.parametrize("degree_l", (2, 3, 4, 7, 10))
def test_exact_matches_the_tables_at_moderate_eccentricity(degree_l):
    """At e = 0.2 the exact functions match level 50 on every shared mode."""
    exact = dict(eccentricity_func(0.2, degree_l, "exact", exact_tolerance=1.0e-12)[0])
    table = dict(eccentricity_func(0.2, degree_l, 50)[0])
    shared = set(exact) & set(table)
    assert len(shared) > 0
    for key in shared:
        assert abs(exact[key] - table[key]) <= 1.0e-11 * max(1.0, abs(table[key])), key


def test_exact_k0_mode_is_the_closed_form():
    """The exact G_210 is (1 - e^2)^(-3/2) and the vanishing k = 0 modes are omitted."""
    e = 0.7
    by_lpq = dict(eccentricity_func(e, 2, "exact")[0])
    assert isclose(by_lpq[(2, 1, 0)], (1.0 - e**2)**-1.5, rel_tol=1.0e-12)
    assert (2, 0, -2) not in by_lpq and (2, 2, 2) not in by_lpq


@pytest.mark.parametrize("eccentricity", (0.1, 0.5, 0.7, 0.8, 0.9))
def test_exact_heating_matches_hut(eccentricity):
    """The exact tolerance bounds the relative error of the Hut (1981) heating."""
    tolerance = 1.0e-4
    squares = eccentricity_squared_func(eccentricity, 2, "exact", exact_tolerance=tolerance)[0]
    error = _synchronous_mode_sum(squares) / _hut_synchronous(eccentricity) - 1.0
    assert -tolerance <= error <= 1.0e-10


def test_exact_mode_range_grows_with_eccentricity_and_tolerance():
    """The exact functions keep more q modes at higher e and tighter tolerance."""
    def max_q(e, tolerance):
        return max(abs(q) for (l, p, q), _ in eccentricity_func(e, 2, "exact", exact_tolerance=tolerance)[0])
    assert max_q(0.3, 1.0e-4) < max_q(0.6, 1.0e-4) < max_q(0.9, 1.0e-4)
    assert max_q(0.6, 1.0e-2) < max_q(0.6, 1.0e-6)
    assert max_q(0.0, 1.0e-4) == 0


def test_recommend_eccentricity_truncation():
    """The recommended level for a few eccentricities and tolerances."""
    assert recommend_eccentricity_truncation(0.0) == 2
    assert recommend_eccentricity_truncation(0.05) == 4
    assert recommend_eccentricity_truncation(0.25, tolerance=1.0e-2) == 10
    assert recommend_eccentricity_truncation(0.3, tolerance=1.0e-2) == 20
    assert recommend_eccentricity_truncation(0.3, tolerance=1.0e-2, max_degree_l=3) == 20
    assert recommend_eccentricity_truncation(0.4, tolerance=1.0e-6) == 50
    assert recommend_eccentricity_truncation(0.45, tolerance=1.0e-6) == "exact"
    assert recommend_eccentricity_truncation(5.0e-4, tolerance=1.0e-4) == 2
    assert recommend_eccentricity_truncation(0.8) == "exact"
    assert recommend_eccentricity_truncation(0.3, tolerance=1.0e-12) == "exact"
    with pytest.raises(ValueError):
        recommend_eccentricity_truncation(1.0)


@pytest.mark.parametrize("level", ECCENTRICITY_TRUNCATIONS)
def test_recommended_levels_hold_their_tolerance(level):
    """At each level's own 1% limit the Hut (1981) heating is within 1%."""
    limit = eccentricity_accuracy_limit(level, 1.0e-2)
    assert 0.0 < limit < 0.8
    error = _synchronous_mode_sum(eccentricity_squared_func(limit, 2, level)[0]) / _hut_synchronous(limit) - 1.0
    assert abs(error) < 1.0e-2
    assert np.isnan(eccentricity_accuracy_limit("exact", 1.0e-2))
