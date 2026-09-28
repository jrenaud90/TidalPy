"""The angular Gram table ``angular_gram`` against an independent high-precision quadrature."""
import math

import numpy as np
import pytest

sympy = pytest.importorskip("sympy")
mpmath = pytest.importorskip("mpmath")

from TidalPy.Tides.multilayer.stress_strain import angular_gram


def _basis(degree_l, order_m):
    """The six angular basis functions in x = cos(theta)."""
    x = sympy.symbols("x", real=True)
    sin_t = sympy.sqrt(1 - x ** 2)
    cot_t = x / sin_t
    P = sympy.assoc_legendre(degree_l, order_m, x)
    dPdx = sympy.diff(P, x)
    dP = -sin_t * dPdx                                   # dP/dtheta
    d2P = (1 - x ** 2) * sympy.diff(P, x, 2) - x * dPdx  # d2P/dtheta2
    funcs = [P, dP, d2P, P / sin_t,
             -order_m ** 2 * P / sin_t ** 2 + cot_t * dP,
             (dP - cot_t * P) / sin_t]
    return x, [sympy.simplify(f) for f in funcs]


def _reference_gram(degree_l, order_m):
    mpmath.mp.dps = 40
    x, funcs = _basis(degree_l, order_m)
    lam = [sympy.lambdify(x, f, "mpmath") for f in funcs]
    gram = np.zeros((6, 6))
    for i in range(6):
        for j in range(i, 6):
            if order_m == 0 and (i in (3, 5) or j in (3, 5)):
                continue  # f4 and f6 diverge for m = 0 but are multiplied by i m = 0 (tabulated as 0)
            value = float(mpmath.quad(lambda xx, i=i, j=j: lam[i](xx) * lam[j](xx), [-1, 0, 1]))
            gram[i, j] = gram[j, i] = value
    return gram


@pytest.mark.parametrize("degree_l", range(2, 11))
def test_table_matches_reference(degree_l):
    """Each order's table is symmetric, matches the reference, and zeroes the unused m = 0 rows."""
    for order_m in range(0, degree_l + 1):
        table = angular_gram(degree_l, order_m)
        assert table.shape == (6, 6)
        np.testing.assert_allclose(table, table.T, rtol=0, atol=0)
        reference = _reference_gram(degree_l, order_m)
        scale = np.maximum(np.abs(reference), 1.0)
        np.testing.assert_allclose(table, reference, rtol=1e-10, atol=1e-10 * scale.max())
        if order_m == 0:
            assert np.all(table[3, :] == 0.0) and np.all(table[:, 3] == 0.0)
            assert np.all(table[5, :] == 0.0) and np.all(table[:, 5] == 0.0)


def test_low_order_closed_forms():
    """Two entries with simple closed forms."""
    assert math.isclose(angular_gram(2, 1)[0, 3], 9.0 * math.pi / 8.0, rel_tol=1e-12)
    assert math.isclose(angular_gram(2, 2)[4, 4], 1032.0 / 5.0, rel_tol=1e-12)


@pytest.mark.parametrize("degree_l, order_m", [(11, 0), (2, 3)], ids=["degree_above_10", "order_above_degree"])
def test_out_of_range_raises(degree_l, order_m):
    with pytest.raises(ValueError):
        angular_gram(degree_l, order_m)
