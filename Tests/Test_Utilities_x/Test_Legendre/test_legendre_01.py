"""Tests for the associated-Legendre utilities (TidalPy.Utilities_x.legendre).

Both the precomputed closed-form tables (``legendre``, degrees l = 2..10) and the generic xsf-backed
evaluator (``legendre_generic``, any degree) must reproduce, for every supported (l, m):

* the value P_lm(cos theta) and its first/second colatitude derivatives vs ``scipy.special``
  (the same unnormalized, Condon-Shortley-phase, branch_cut=2 convention);
* each other (table vs generic);
* the hand-coded l = 2 reference forms exactly.
"""
import numpy as np
import pytest

from scipy.special import assoc_legendre_p

from TidalPy.Utilities_x.legendre import legendre, legendre_generic


# Colatitudes away from the exact poles (the m >= 1 derivatives have removable 1/sin theta structure
# at theta = 0, pi, as in the legacy hand-coded forms).
_COLATS = np.linspace(0.15, np.pi - 0.15, 13)
_LM = [(l, m) for l in range(2, 11) for m in range(0, l + 1)]


def _scipy_triple(l, m, colat):
    """(P, dP/dtheta, d2P/dtheta2) from scipy via the chain rule x = cos(theta)."""
    x = np.cos(colat)
    s = np.sin(colat)
    p, dpdx, d2pdx2 = assoc_legendre_p(l, m, x, branch_cut=2, diff_n=2)
    return (float(p),
            float(-s * dpdx),
            float(s * s * d2pdx2 - x * dpdx))


@pytest.mark.parametrize("l,m", _LM)
def test_table_matches_scipy(l, m):
    for colat in _COLATS:
        got = legendre(l, m, float(colat))
        ref = _scipy_triple(l, m, colat)
        assert np.allclose(got, ref, rtol=1e-11, atol=1e-11), \
            f"table l={l} m={m} colat={colat}: {got} vs scipy {ref}"


@pytest.mark.parametrize("l,m", _LM)
def test_generic_matches_scipy(l, m):
    for colat in _COLATS:
        got = legendre_generic(l, m, float(colat))
        ref = _scipy_triple(l, m, colat)
        assert np.allclose(got, ref, rtol=1e-11, atol=1e-11), \
            f"generic l={l} m={m} colat={colat}: {got} vs scipy {ref}"


@pytest.mark.parametrize("l,m", _LM)
def test_table_matches_generic(l, m):
    for colat in _COLATS:
        assert np.allclose(legendre(l, m, float(colat)),
                           legendre_generic(l, m, float(colat)),
                           rtol=1e-11, atol=1e-11)


def test_l2_hardcoded_reference():
    """The l = 2 table reproduces the hand-coded forms used by the tidal potential kernel."""
    for colat in _COLATS:
        c = np.cos(colat)
        s = np.sin(colat)
        # (P, dP/dtheta, d2P/dtheta2) for m = 0, 1, 2.
        ref = {
            0: (0.5 * (3 * c * c - 1.0), -3.0 * c * s, 3.0 * (s * s - c * c)),
            1: (-3.0 * c * s, 3.0 * (s * s - c * c), 12.0 * c * s),
            2: (3.0 * s * s, 6.0 * c * s, 6.0 * (c * c - s * s)),
        }
        for m, expected in ref.items():
            assert np.allclose(legendre(2, m, float(colat)), expected, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("func", [legendre, legendre_generic])
def test_out_of_range_order_is_nan(func):
    # m > l and m < 0 are undefined -> NaN triple.
    assert all(np.isnan(func(2, 3, 0.7)))
    assert all(np.isnan(func(3, -1, 0.7)))


def test_generic_supports_high_degree():
    """The generic evaluator works beyond the table range (l = 11) and matches scipy."""
    for colat in _COLATS:
        got = legendre_generic(11, 4, float(colat))
        ref = _scipy_triple(11, 4, colat)
        assert np.allclose(got, ref, rtol=1e-10, atol=1e-10)
    # The table does not cover l = 11.
    assert all(np.isnan(legendre(11, 4, 0.7)))
