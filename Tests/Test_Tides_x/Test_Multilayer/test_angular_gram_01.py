"""Angular Gram table (``TidalPy.Tides_x.multilayer.stress_strain.angular_gram``).

The table stores ``G_ij(l, m) = int_0^pi f_i(theta) f_j(theta) sin(theta) dtheta`` for the bounded
6-function angular basis used by the analytic colatitude collapse. This test recomputes each entry
independently (sympy Legendre forms + high-precision tanh-sinh quadrature in ``x = cos theta``) and
checks the compiled table matches, plus the structural properties (symmetry, the ``m = 0`` zeros).
"""
import math

import numpy as np
import pytest

sympy = pytest.importorskip("sympy")
mpmath = pytest.importorskip("mpmath")

from TidalPy.Tides_x.multilayer.stress_strain import angular_gram


def _basis(degree_l, order_m):
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
                continue  # f4/f6 diverge for m=0 but are multiplied by i*m=0 (unused, tabulated as 0)
            value = float(mpmath.quad(lambda xx, i=i, j=j: lam[i](xx) * lam[j](xx), [-1, 0, 1]))
            gram[i, j] = gram[j, i] = value
    return gram


@pytest.mark.parametrize("degree_l", range(2, 11))
def test_table_matches_reference(degree_l):
    for order_m in range(0, degree_l + 1):
        table = angular_gram(degree_l, order_m)
        assert table.shape == (6, 6)
        # symmetric
        np.testing.assert_allclose(table, table.T, rtol=0, atol=0)
        reference = _reference_gram(degree_l, order_m)
        scale = np.maximum(np.abs(reference), 1.0)
        np.testing.assert_allclose(table, reference, rtol=1e-10, atol=1e-10 * scale.max())
        if order_m == 0:
            # f4 (index 3) and f6 (index 5) entries are unused and stored as zero.
            assert np.all(table[3, :] == 0.0) and np.all(table[:, 3] == 0.0)
            assert np.all(table[5, :] == 0.0) and np.all(table[:, 5] == 0.0)


def test_low_order_closed_forms():
    """A couple of entries have simple closed forms; check the table hits them."""
    g21 = angular_gram(2, 1)
    assert math.isclose(g21[0, 3], 9.0 * math.pi / 8.0, rel_tol=1e-12)  # int P f4 sin dtheta = 9 pi / 8
    g22 = angular_gram(2, 2)
    assert math.isclose(g22[4, 4], 1032.0 / 5.0, rel_tol=1e-12)         # int f5^2 sin dtheta = 206.4


def test_out_of_range_raises():
    with pytest.raises(ValueError):
        angular_gram(11, 0)   # l > 10
    with pytest.raises(ValueError):
        angular_gram(2, 3)    # m > l
