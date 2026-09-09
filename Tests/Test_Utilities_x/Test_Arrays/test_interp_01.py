"""
Tests for TidalPy.Utilities_x.arrays.interp — 1-D linear interpolation.

Cross-checks the C++-backed ``interp`` against ``numpy.interp`` for scalar and
array queries, endpoint clamping, shape preservation, and small domains.

Requires the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import numpy as np
import pytest


def _import_interp():
    from TidalPy.Utilities_x.arrays.interp import interp
    return interp


_XP = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
_FP = np.array([0.0, 10.0, 5.0, 7.0, 1.0, 9.0])


def test_scalar_matches_numpy():
    interp = _import_interp()
    for x in (0.5, 1.0, 2.5, 3.9, 4.999):
        assert interp(x, _XP, _FP) == pytest.approx(float(np.interp(x, _XP, _FP)))


def test_scalar_returns_float():
    interp = _import_interp()
    assert isinstance(interp(2.5, _XP, _FP), float)


def test_endpoint_clamping():
    interp = _import_interp()
    # Below/above the domain clamps to the endpoint values (numpy.interp default).
    assert interp(-10.0, _XP, _FP) == pytest.approx(_FP[0])
    assert interp(100.0, _XP, _FP) == pytest.approx(_FP[-1])
    # Exact node values.
    assert interp(0.0, _XP, _FP) == pytest.approx(_FP[0])
    assert interp(5.0, _XP, _FP) == pytest.approx(_FP[-1])
    assert interp(3.0, _XP, _FP) == pytest.approx(_FP[3])


def test_array_matches_numpy():
    interp = _import_interp()
    x = np.linspace(-1.0, 6.0, 50)
    got = interp(x, _XP, _FP)
    assert isinstance(got, np.ndarray)
    np.testing.assert_allclose(got, np.interp(x, _XP, _FP), rtol=1e-12, atol=0.0)


def test_shape_preserved():
    interp = _import_interp()
    x = np.array([[0.5, 1.5], [2.5, 3.5]])
    got = interp(x, _XP, _FP)
    assert got.shape == x.shape
    np.testing.assert_allclose(got, np.interp(x.ravel(), _XP, _FP).reshape(x.shape), rtol=1e-12)


def test_two_point_domain():
    interp = _import_interp()
    xp = np.array([0.0, 10.0])
    fp = np.array([100.0, 200.0])
    assert interp(5.0, xp, fp) == pytest.approx(150.0)
    assert interp(-1.0, xp, fp) == pytest.approx(100.0)
    assert interp(11.0, xp, fp) == pytest.approx(200.0)


def test_single_point_domain():
    interp = _import_interp()
    assert interp(3.0, np.array([1.0]), np.array([42.0])) == pytest.approx(42.0)


def test_length_mismatch_raises():
    interp = _import_interp()
    with pytest.raises(ValueError):
        interp(1.0, np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0]))


def test_empty_domain_raises():
    interp = _import_interp()
    with pytest.raises(ValueError):
        interp(1.0, np.array([]), np.array([]))


def test_accepts_python_lists():
    interp = _import_interp()
    assert interp(0.5, [0.0, 1.0], [0.0, 2.0]) == pytest.approx(1.0)
