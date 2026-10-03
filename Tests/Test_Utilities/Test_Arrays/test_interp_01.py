"""Tests for ``TidalPy.Utilities.arrays.interp`` (1D linear interpolation) against ``numpy.interp``."""
import numpy as np
import pytest

from TidalPy.Utilities.arrays.interp import interp

_XP = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
_FP = np.array([0.0, 10.0, 5.0, 7.0, 1.0, 9.0])


@pytest.mark.parametrize("x", [0.5, 1.0, 2.5, 3.9, 4.999])
def test_scalar_matches_numpy(x):
    """A scalar query matches numpy.interp."""
    assert interp(x, _XP, _FP) == pytest.approx(float(np.interp(x, _XP, _FP)))


def test_scalar_returns_float():
    """A scalar query returns a Python float."""
    assert isinstance(interp(2.5, _XP, _FP), float)


@pytest.mark.parametrize(
    "x, expected",
    [(-10.0, _FP[0]), (100.0, _FP[-1]), (0.0, _FP[0]), (5.0, _FP[-1]), (3.0, _FP[3])],
    ids=["below", "above", "first_node", "last_node", "interior_node"])
def test_endpoint_clamping(x, expected):
    """Queries outside the domain clamp to the endpoints; node queries return the node value."""
    assert interp(x, _XP, _FP) == pytest.approx(expected)


def test_array_matches_numpy():
    """An array query returns an ndarray matching numpy.interp."""
    x = np.linspace(-1.0, 6.0, 50)
    got = interp(x, _XP, _FP)
    assert isinstance(got, np.ndarray)
    np.testing.assert_allclose(got, np.interp(x, _XP, _FP), rtol=1e-12, atol=0.0)


def test_shape_preserved():
    """A 2D query keeps its shape."""
    x = np.array([[0.5, 1.5], [2.5, 3.5]])
    got = interp(x, _XP, _FP)
    assert got.shape == x.shape
    np.testing.assert_allclose(got, np.interp(x.ravel(), _XP, _FP).reshape(x.shape), rtol=1e-12)


@pytest.mark.parametrize(
    "x, xp, fp, expected",
    [
        (5.0, [0.0, 10.0], [100.0, 200.0], 150.0),
        (-1.0, [0.0, 10.0], [100.0, 200.0], 100.0),
        (11.0, [0.0, 10.0], [100.0, 200.0], 200.0),
        (3.0, [1.0], [42.0], 42.0),
    ],
    ids=["two_point_inside", "two_point_below", "two_point_above", "single_point"])
def test_small_domain(x, xp, fp, expected):
    """Two point and single point domains interpolate or clamp correctly."""
    assert interp(x, np.array(xp), np.array(fp)) == pytest.approx(expected)


@pytest.mark.parametrize(
    "xp, fp",
    [([0.0, 1.0, 2.0], [0.0, 1.0]), ([], [])],
    ids=["length_mismatch", "empty"])
def test_bad_domain_raises(xp, fp):
    """A length mismatch or an empty domain raises ValueError."""
    with pytest.raises(ValueError):
        interp(1.0, np.array(xp), np.array(fp))


def test_accepts_python_lists():
    """Python lists are accepted for the domain."""
    assert interp(0.5, [0.0, 1.0], [0.0, 2.0]) == pytest.approx(1.0)


def test_accepts_read_only_arrays():
    """Read-only inputs are accepted for array and scalar queries."""
    xp = _XP.copy()
    fp = _FP.copy()
    x = np.array([0.5, 2.5, 4.5])
    for array in (xp, fp, x):
        array.setflags(write=False)
    np.testing.assert_allclose(interp(x, xp, fp), np.interp(x, xp, fp))
    assert interp(1.5, xp, fp) == pytest.approx(float(np.interp(1.5, xp, fp)))


def test_complex_values_keep_their_imaginary_part():
    """Complex samples interpolate both parts, as numpy.interp does, rather than dropping the imaginary one."""
    fp = _FP + 1j * (2.0 * _FP + 1.0)
    x = np.array([0.5, 2.5, 4.5, -1.0, 10.0])
    np.testing.assert_allclose(interp(x, _XP, fp), np.interp(x, _XP, fp))
    assert interp(1.5, _XP, fp) == pytest.approx(complex(np.interp(1.5, _XP, fp)))
    assert isinstance(interp(1.5, _XP, fp), complex)
