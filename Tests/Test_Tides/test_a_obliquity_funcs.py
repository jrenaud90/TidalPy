"""Obliquity functions F_lmp(I): each truncation level, their cut squares, and the truncation helpers."""
from math import isclose
import warnings

import pytest
import numpy as np

from TidalPy.Tides.obliquity import (
    OBLIQUITY_GENERAL, OBLIQUITY_TRUNCATIONS, obliquity_accuracy_limit, obliquity_func, obliquity_squared_func,
    obliquity_truncation_name, promote_obliquity_truncation, recommend_obliquity_truncation,
    validate_obliquity_truncation)
from TidalPy.Tides.classes.collapse import collapse_global_tides
from TidalPy.Utilities.lookups import IntMap3

DEGREES = (2, 3, 4, 5, 6, 7, 8, 9, 10)
OBLIQUITY = 0.5
_HALF_SIN = np.sin(OBLIQUITY / 2)
_HALF_COS = np.cos(OBLIQUITY / 2)

# Degree-2 values at I = 0.5: total (l, m, p) count, (l, m, p) values, (l, m) count, (l, m) sizes, (l, m) values.
_DEGREE_TWO_SPOT_CHECKS = {
    'off': (
        2,
        {(2, 0, 1): -0.5, (2, 2, 0): 3.0},
        2,
        {(2, 0): 1, (2, 2): 1},
        {(2, 0, 1): -0.5, (2, 2, 0): 3.0}),
    2: (
        7,
        {(2, 0, 1): 0.75 * OBLIQUITY**2 - 0.5, (2, 1, 0): 1.5 * OBLIQUITY},
        3,
        {(2, 1): 2, (2, 2): 2},
        {(2, 1, 0): 1.5 * OBLIQUITY, (2, 1, 1): -1.5 * OBLIQUITY, (2, 2, 0): 3.0 - 1.5 * OBLIQUITY**2}),
    4: (
        9,
        {(2, 0, 1): -0.5 + 0.75 * OBLIQUITY**2 - 0.25 * OBLIQUITY**4,
         (2, 2, 0): 3.0 - 1.5 * OBLIQUITY**2 + 0.3125 * OBLIQUITY**4,
         (2, 2, 2): 0.1875 * OBLIQUITY**4},
        None,
        {},
        {}),
    'gen': (
        9,
        {(2, 0, 1): -_HALF_SIN**4 + _HALF_SIN**2 + 0.5 * np.sin(OBLIQUITY)**2 - 0.5,
         (2, 1, 0): 3.0 * _HALF_SIN * _HALF_COS**3},
        3,
        {(2, 1): 3, (2, 2): 3},
        {(2, 2, 0): 3.0 * _HALF_COS**4, (2, 2, 2): 3.0 * _HALF_SIN**4}),
}


@pytest.mark.parametrize('degree_l, truncation', [(1, truncation) for truncation in ('gen', 2, 4, 'off', 1, 3, 10)]
                         + [(degree_l, truncation) for degree_l in DEGREES for truncation in (1, 3, 10)])
def test_obliquity_func_unsupported_raises(degree_l, truncation):
    """Degree 1 and untabulated levels (the old 1 and 10 among them) raise."""
    with pytest.raises(NotImplementedError):
        obliquity_func(OBLIQUITY, degree_l, truncation)


@pytest.mark.parametrize('degree_l', DEGREES)
@pytest.mark.parametrize('truncation', ('gen', 2, 4, 'off'))
def test_obliquity_func_structure(degree_l, truncation):
    """Both return maps hold int keys within the degree and float values."""
    by_lmp, by_lm = obliquity_func(OBLIQUITY, degree_l, truncation)
    assert isinstance(by_lmp, IntMap3)
    assert isinstance(by_lm, dict)

    assert len(by_lmp) > 0
    for (l, m, p), value in by_lmp:
        assert isinstance(l, int)
        assert isinstance(m, int)
        assert isinstance(p, int)
        assert l == degree_l
        assert m <= l
        assert p <= l
        assert isinstance(value, float)

    for (l, m), by_p in by_lm.items():
        assert isinstance(l, int)
        assert isinstance(m, int)
        assert l == degree_l
        assert m <= l
        for (p,), value in by_p:
            assert isinstance(p, int)
            assert isinstance(value, float)
            assert p <= degree_l


@pytest.mark.parametrize('truncation', tuple(_DEGREE_TWO_SPOT_CHECKS))
def test_obliquity_func_degree_two_values(truncation):
    """Degree-2 mode counts and values at each level."""
    num_lmp, lmp_values, num_lm, lm_sizes, lm_values = _DEGREE_TWO_SPOT_CHECKS[truncation]
    by_lmp, by_lm = obliquity_func(OBLIQUITY, 2, truncation)
    assert len(by_lmp) == num_lmp
    for key, expected in lmp_values.items():
        assert isclose(by_lmp[key], expected), key
    if num_lm is not None:
        assert len(by_lm) == num_lm
    for key, size in lm_sizes.items():
        assert len(by_lm[key]) == size, key
    for (l, m, p), expected in lm_values.items():
        assert isclose(by_lm[(l, m)][(p,)], expected), (l, m, p)


@pytest.mark.parametrize('degree_l', DEGREES)
@pytest.mark.parametrize('truncation', (2, 4))
def test_obliquity_truncation_order(degree_l, truncation):
    """Level N carries F through I^N: halving I shrinks the gap to the general form by at least 2^(N+1)."""
    large_obliquity = 0.02
    small_obliquity = 0.01
    gaps = dict()
    for obliquity in (large_obliquity, small_obliquity):
        exact = dict(obliquity_func(obliquity, degree_l, 'gen')[0])
        truncated = dict(obliquity_func(obliquity, degree_l, truncation)[0])
        assert set(truncated) <= set(exact)
        # A mode the truncation drops counts as zero.
        gaps[obliquity] = {key: abs(truncated.get(key, 0.0) - value) for key, value in exact.items()}

    num_checked = 0
    for key, gap_small in gaps[small_obliquity].items():
        # Gaps at round-off carry no ratio information.
        if gap_small < 1.0e-12 * max(1.0, abs(dict(obliquity_func(small_obliquity, degree_l, 'gen')[0])[key])):
            continue
        num_checked += 1
        assert gaps[large_obliquity][key] / gap_small > 0.9 * 2**(truncation + 1), key
    assert num_checked > 0


@pytest.mark.parametrize('degree_l', DEGREES)
@pytest.mark.parametrize('truncation', (0, 2, 4))
def test_obliquity_squares_are_cut_at_the_level(degree_l, truncation):
    """F^2 cut at I^N keeps only modes starting at or below I^(N/2) and converges at 2^(N+1); general squares are
    plain squares."""
    squares = {obliquity: dict(obliquity_squared_func(obliquity, degree_l, truncation)[0])
               for obliquity in (0.02, 0.01)}
    for (l, m, p) in squares[0.01]:
        assert 2 * abs(l - m - 2 * p) <= truncation
    exact = {obliquity: {key: value**2 for key, value in obliquity_func(obliquity, degree_l, 'gen')[0]}
             for obliquity in (0.02, 0.01)}
    num_checked = 0
    for key, value in exact[0.01].items():
        gap_small = abs(squares[0.01].get(key, 0.0) - value)
        if gap_small < 1.0e-12 * max(1.0, value):
            continue
        num_checked += 1
        gap_large = abs(squares[0.02].get(key, 0.0) - exact[0.02][key])
        assert gap_large / gap_small > 0.9 * 2**(truncation + 1), key
    assert num_checked > 0
    general = dict(obliquity_squared_func(0.4, degree_l, 'gen')[0])
    for key, value in obliquity_func(0.4, degree_l, 'gen')[0]:
        assert isclose(general[key], value**2, rel_tol=1.0e-14)


def test_level_two_squares_at_degree_two():
    """Level-2 degree-2 squares hold the expected modes and values."""
    obliquity = 0.2
    squares = dict(obliquity_squared_func(obliquity, 2, 2)[0])
    assert set(squares) == {(2, 0, 1), (2, 1, 0), (2, 1, 1), (2, 2, 0)}
    assert isclose(squares[(2, 0, 1)], 0.25 - 0.75 * obliquity**2)
    assert isclose(squares[(2, 1, 0)], 2.25 * obliquity**2)
    assert isclose(squares[(2, 2, 0)], 9.0 - 9.0 * obliquity**2)


def test_truncation_names_and_promotion():
    """Truncation validation, naming, and promotion of old configured levels."""
    assert OBLIQUITY_TRUNCATIONS == (0, 2, 4)
    assert validate_obliquity_truncation("off") == 0
    assert validate_obliquity_truncation("gen") == OBLIQUITY_GENERAL
    assert validate_obliquity_truncation("General") == OBLIQUITY_GENERAL
    assert validate_obliquity_truncation("4") == 4
    assert validate_obliquity_truncation(2.0) == 2
    assert validate_obliquity_truncation(OBLIQUITY_GENERAL) == OBLIQUITY_GENERAL
    for bad in (1, 3, 10, "1"):
        with pytest.raises(NotImplementedError):
            validate_obliquity_truncation(bad)
    for bad in (True, 2.5):
        with pytest.raises(TypeError):
            validate_obliquity_truncation(bad)
    assert obliquity_truncation_name(OBLIQUITY_GENERAL) == "gen"
    assert obliquity_truncation_name(2) == 2
    # The old general code 10 becomes the general functions.
    for level, promoted in ((1, 2), (3, 4), (10, OBLIQUITY_GENERAL)):
        with pytest.warns(UserWarning):
            assert promote_obliquity_truncation(level, warned_levels=set(), warn=True) == promoted
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert promote_obliquity_truncation("gen", warned_levels=set(), warn=True) == OBLIQUITY_GENERAL
        assert promote_obliquity_truncation(4, warned_levels=set(), warn=True) == 4
    with pytest.raises(ValueError):
        promote_obliquity_truncation(-2, warned_levels=set(), warn=False)


def test_recommend_obliquity_truncation():
    """The recommended level and accuracy limit for a few obliquities and tolerances."""
    assert recommend_obliquity_truncation(0.0) == 0
    assert recommend_obliquity_truncation(0.1) == 2
    assert recommend_obliquity_truncation(-0.1) == 2
    assert recommend_obliquity_truncation(0.1, max_degree_l=3) == 4
    assert recommend_obliquity_truncation(0.3) == 4
    assert recommend_obliquity_truncation(0.6) == "gen"
    assert recommend_obliquity_truncation(0.3, tolerance=1.0e-12) == "gen"
    assert recommend_obliquity_truncation(0.0, tolerance=1.0e-12) == 0
    with pytest.raises(ValueError):
        recommend_obliquity_truncation(0.1, tolerance=0.0)
    assert obliquity_accuracy_limit("off") == 0.0
    assert np.isnan(obliquity_accuracy_limit("gen"))


_BODY = dict(planet_radius=1.8215e6, orbital_frequency=4.11e-5, semi_major_axis=4.217e8, host_mass=1.898e27,
             G_to_use=6.674e-11)


@pytest.mark.parametrize("level", (2, 4))
@pytest.mark.parametrize("tide_model, config", [("cpl", {"fixed_k": [0.3], "fixed_q": [100.0]}),
                                                ("ctl", {"fixed_k": [0.3], "fixed_dt_s": [100.0]})])
@pytest.mark.parametrize("spin_ratio", (0.5, 1.0, 2.3))
def test_levels_hold_their_tolerance(level, tide_model, config, spin_ratio):
    """At each level's own 1% limit the collapsed 1D heating is within 1% of the general functions'."""
    limit = obliquity_accuracy_limit(level, 1.0e-2)
    assert 0.0 < limit < 0.5

    def heating(truncation):
        return collapse_global_tides(
            **_BODY,
            spin_frequency=spin_ratio * _BODY["orbital_frequency"],
            eccentricity=0.0,
            obliquity=limit,
            tide_model=tide_model,
            tide_config=config,
            eccentricity_truncation=2,
            obliquity_truncation=truncation)["tidal_heating"]

    assert abs(heating(level) / heating("gen") - 1.0) < 1.0e-2
