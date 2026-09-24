"""TidalPy's tabulated eccentricity functions G_lpq(e).

Truncation level N keeps every product of two eccentricity functions through e^N: the unsquared functions hold every
mode with |q| <= N, each through e^N, and the cut squares every mode with |q| <= N / 2, each square through e^N (see
Documentation/Tides_x/eccentricity.md).
"""
from math import isclose

import pytest

from TidalPy.Tides_x.eccentricity import (
    ECCENTRICITY_TRUNCATIONS, eccentricity_func, eccentricity_squared_func, validate_eccentricity_truncation)
from TidalPy.Utilities_x.lookups import IntMap3

DEGREES = (2, 3, 4, 5, 6, 7, 8, 9, 10)


def _hut_synchronous(eccentricity):
    """Hut (1981) constant-time-lag heating at synchronous rotation, in the units of the l = 2 mode sum below,
    sum (3/4) G_20q^2 q^2 + (1/4) G_21q^2 q^2 (the ratio 3 follows from the small-e limits 10.5 e^2 and 3.5 e^2)."""
    e2 = eccentricity * eccentricity
    n_a = (1 + 31 / 2 * e2 + 255 / 8 * e2**2 + 185 / 16 * e2**3 + 25 / 64 * e2**4) / (1 - e2)**7.5
    n_e = (1 + 15 / 2 * e2 + 45 / 8 * e2**2 + 5 / 16 * e2**3) / (1 - e2)**6
    omega_e = (1 + 3 * e2 + 3 / 8 * e2**2) / (1 - e2)**4.5
    return 3.0 * (n_a - 2.0 * n_e + omega_e)


def _synchronous_mode_sum(eccentricity, truncation):
    """The l = 2, zero-obliquity, synchronous constant-time-lag heating from the cut squares: each mode's frequency
    is q n."""
    by_lpq, _ = eccentricity_squared_func(eccentricity, 2, truncation)
    weight = {0: 0.75, 1: 0.25}
    return sum(weight[p] * value * q**2 for (l, p, q), value in by_lpq if p in weight)


def _vanishing(degree_l, max_q):
    """Modes with k = l - 2p + q = 0 and |l - 2p| >= l vanish identically (X^{-(l+1), m}_0 = 0 for |m| >= l)."""
    return sum(1 for p in (0, degree_l) if degree_l <= max_q)


def test_tabulated_levels():
    assert ECCENTRICITY_TRUNCATIONS == (2, 4, 6, 8, 10, 20, 50)
    assert validate_eccentricity_truncation("20") == 20
    for bad in (0, 1, 3, 7, 60):
        with pytest.raises(NotImplementedError):
            validate_eccentricity_truncation(bad)
    with pytest.raises(TypeError):
        validate_eccentricity_truncation(True)


@pytest.mark.parametrize('degree_l', DEGREES)
@pytest.mark.parametrize('truncation', (2, 4, 6, 10, 0))
def test_eccentricity_funcs(degree_l, truncation):
    """Tests TidalPy's eccentricity functions for various degree_ls and various truncation levels."""

    eccentricity = 0.5

    if truncation == 0:
        # Truncation 0 is not tabulated, check that it raises an error.
        with pytest.raises(NotImplementedError):
            eccentricity_func(eccentricity, degree_l, truncation)
        return

    ecc_results_by_lpq, ecc_results_by_lp = eccentricity_func(eccentricity, degree_l, truncation)

    # Check return types
    assert isinstance(ecc_results_by_lpq, IntMap3)
    assert isinstance(ecc_results_by_lp, dict)

    assert len(ecc_results_by_lpq) > 0
    for (l, p, q), ecc_result in ecc_results_by_lpq:
        assert isinstance(l, int)
        assert isinstance(p, int)
        assert isinstance(q, int)
        assert l == degree_l
        assert p <= l
        assert abs(q) <= truncation
        assert isinstance(ecc_result, float)

    for (l, p), ecc_results_by_q in ecc_results_by_lp.items():
        assert isinstance(l, int)
        assert isinstance(p, int)
        assert l == degree_l
        assert p <= l
        for (q,), ecc_result in ecc_results_by_q:
            assert isinstance(q, int)
            assert isinstance(ecc_result, float)

    # Spot checks (Kaula 1964; Cayley 1861).
    e = eccentricity
    if degree_l == 2 and truncation == 2:
        # 3 (2 * 2 + 1) modes less G_20,-2 and G_22,2, which vanish identically; each G through e^2.
        assert len(ecc_results_by_lpq) == 13
        assert isclose(ecc_results_by_lpq[(2, 0, 0)], 1.0 - 2.5 * e**2)
        assert isclose(ecc_results_by_lpq[(2, 0, -1)], -0.5 * e)
        assert isclose(ecc_results_by_lpq[(2, 0, 1)], 3.5 * e)
        assert isclose(ecc_results_by_lpq[(2, 0, 2)], 8.5 * e**2)
        assert isclose(ecc_results_by_lpq[(2, 1, 0)], (1.0 - e**2)**-1.5)
        assert len(ecc_results_by_lp) == 3
        assert len(ecc_results_by_lp[(2, 0)]) == 4
        assert len(ecc_results_by_lp[(2, 1)]) == 5
        assert len(ecc_results_by_lp[(2, 2)]) == 4
        assert isclose(ecc_results_by_lp[(2, 1)][(1,)], 1.5 * e)
        assert isclose(ecc_results_by_lp[(2, 1)][(-2,)], 2.25 * e**2)
    elif degree_l == 2 and truncation == 4:
        # G_200 and G_202 through e^4.
        assert isclose(ecc_results_by_lpq[(2, 0, 0)], 1.0 - 2.5 * e**2 + 0.8125 * e**4)
        assert isclose(ecc_results_by_lpq[(2, 0, 2)], 8.5 * e**2 - (115.0 / 6.0) * e**4)
        assert isclose(ecc_results_by_lpq[(2, 0, 1)], 3.5 * e - 123.0 / 16.0 * e**3)
        assert isclose(ecc_results_by_lpq[(2, 1, 1)], 1.5 * e + 27.0 / 16.0 * e**3)


@pytest.mark.parametrize('degree_l', DEGREES)
@pytest.mark.parametrize('truncation', ECCENTRICITY_TRUNCATIONS)
def test_mode_count(degree_l, truncation):
    """The unsquared functions hold at most (l + 1)(2N + 1) modes and the cut squares at most (l + 1)(N + 1): a mode
    is left out when it vanishes identically or, at a low level, when its series starts past e^N (the leading
    coefficient of a few odd-degree modes is zero). At l = 2 only the two k = 0 modes with |m| = l vanish."""
    by_lpq = dict(eccentricity_func(0.3, degree_l, truncation)[0])
    by_lpq_squared = dict(eccentricity_squared_func(0.3, degree_l, truncation)[0])
    assert len(by_lpq) <= (degree_l + 1) * (2 * truncation + 1)
    assert set(by_lpq_squared) <= set(by_lpq)
    assert all(2 * abs(q) <= truncation for (l, p, q) in by_lpq_squared)
    if degree_l == 2:
        assert len(by_lpq) == 3 * (2 * truncation + 1) - _vanishing(2, truncation)
        assert len(by_lpq_squared) == 3 * (truncation + 1) - _vanishing(2, truncation // 2)


def test_squared_spot_checks():
    """At level 2 the cut squares are the leading-order squares; the k = 0 modes are exact."""
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
    """Level N carries each G_lpq and each cut square through e^N: halving e shrinks the gap to the level-50 table by
    at least 2^(N + 1).

    Higher levels are not checked this way because their gaps at small e are below double precision.
    """
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
            # Skip gaps at round-off, where the ratio carries no information.
            if gap_small < 1.0e-15:
                continue
            num_checked += 1
            assert gaps[large_eccentricity][key] / gap_small > 0.9 * 2**(truncation + 1), (function.__name__, key)
        assert num_checked > 0


@pytest.mark.parametrize('eccentricity', (0.2, 0.4, 0.5, 0.6))
def test_levels_converge_to_hut(eccentricity):
    """Against the exact constant-time-lag heating of Hut (1981), every level errs low and a higher level is never
    worse; level 50 stays within 1e-5 up to e = 0.6."""
    exact = _hut_synchronous(eccentricity)
    errors = [_synchronous_mode_sum(eccentricity, level) / exact - 1.0 for level in ECCENTRICITY_TRUNCATIONS]
    for error in errors:
        assert error <= 1.0e-12
    for lower, higher in zip(errors, errors[1:]):
        assert abs(higher) <= max(abs(lower), 1.0e-13)
    assert abs(errors[-1]) < 1.0e-5
