""" Tests for TidalPy's custom RK integrator.
These integrators are comparable to SciPy's solve_ivp integration scheme but with numba tolerance.
"""
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp

import TidalPy
from TidalPy.utilities.integration.rk_integrator import rk_integrate
from TidalPy.utilities.performance import njit

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

TEVAL_RTOL = 0.8


def test_rk23_coeffs():
    """ Test that the RK23 coeffs are as expected. """

    from TidalPy.utilities.integration.rk_integrator import RK23_A as A
    from TidalPy.utilities.integration.rk_integrator import RK23_B as B
    from TidalPy.utilities.integration.rk_integrator import RK23_C as C

    assert_allclose(np.sum(B), 1, rtol=1e-15)
    assert_allclose(np.sum(A, axis=1), C, rtol=1e-14)

def test_rk45_coeffs():
    """ Test that the RK45 coeffs are as expected. """

    from TidalPy.utilities.integration.rk_integrator import RK45_A as A
    from TidalPy.utilities.integration.rk_integrator import RK45_B as B
    from TidalPy.utilities.integration.rk_integrator import RK45_C as C

    assert_allclose(np.sum(B), 1, rtol=1e-15)
    assert_allclose(np.sum(A, axis=1), C, rtol=1e-14)


def test_rk23_Yis1D_noTeval():
    """ Test the RK order 3 integrator for a 1D diffeq with no explicit Teval. """

    @njit(cacheable=False)
    def diffeq(t, y):
        y0 = y[0]
        dy0 = 0.5 * y0 - np.sin(2. * t)
        return np.asarray([dy0])

    initial_y = np.asarray([1.])
    t_span = (0., 100.)
    ts, ys, success, message = \
        rk_integrate(
            diffeq, t_span, initial_y, rk_method=0,
            rtol=1.0e-3, atol=1.0e-6
            )

    assert success
    assert type(message) == str
    assert type(ts) == np.ndarray
    assert ts.dtype in (float, np.float64)
    assert type(ys) == np.ndarray
    assert ys.dtype in (float, np.float64)
    np.testing.assert_almost_equal(ts[0], t_span[0])
    np.testing.assert_almost_equal(ts[-1], t_span[-1])
    assert ys.shape[0] == 1
    assert ys.shape[1] == ts.shape[0]

    # Compare to scipy's solution
    solution = solve_ivp(diffeq, t_span, initial_y, method='RK23', t_eval=None, rtol=1.0e-3, atol=1.0e-6)
    scipy_y = solution.y
    scipy_t = solution.t

    assert scipy_t.shape == ts.shape
    assert scipy_y.shape == ys.shape
    assert np.allclose(ys, scipy_y, rtol=0.1)


def test_rk45_Yis1D_noTeval():
    """ Test the RK order 5 integrator for a 1D diffeq with no explicit Teval. """

    @njit(cacheable=False)
    def diffeq(t, y):
        y0 = y[0]
        dy0 = 0.5 * y0 - np.sin(2. * t)
        return np.asarray([dy0])
    initial_y = np.asarray([1.])
    t_span = (0., 100.)
    ts, ys, success, message = \
        rk_integrate(
            diffeq, t_span, initial_y, rk_method=1,
            rtol=1.0e-5, atol=1.0e-6,
            )

    assert success
    assert type(message) == str
    assert type(ts) == np.ndarray
    assert ts.dtype in (float, np.float64)
    assert type(ys) == np.ndarray
    assert ys.dtype in (float, np.float64)
    np.testing.assert_almost_equal(ts[0], t_span[0])
    np.testing.assert_almost_equal(ts[-1], t_span[-1])
    assert ys.shape[0] == 1
    assert ys.shape[1] == ts.shape[0]

    # Compare to scipy's solution
    solution = solve_ivp(diffeq, t_span, initial_y, method='RK45', t_eval=None, rtol=1.0e-5, atol=1.0e-6)
    scipy_y = solution.y
    scipy_t = solution.t

    assert scipy_t.shape == ts.shape
    assert scipy_y.shape == ys.shape

    assert np.allclose(ys, scipy_y, rtol=0.1)


def test_rk23_Yis2D_noTeval():
    """ Test the RK order 3 integrator for a 2D diffeq with no explicit Teval. """

    @njit(cacheable=False)
    def diffeq(t, y):
        y0 = y[0]
        y1 = y[1]
        y2 = y[2]
        dy0 = 0.5 * y0 - np.sin(2. * t)
        dy1 = 0.5 * y1 - np.sin(2. * t) - 0.005 * y0
        dy2 = 0.5 * y2 - np.sin(2. * t) + 0.005 * y1
        return np.asarray([dy0, dy1, dy2])

    initial_y = np.asarray([1., (1. - 0.005), (1. - (1. - 0.005))])
    t_span = (0., 100.)
    ts, ys, success, message = \
        rk_integrate(
            diffeq, t_span, initial_y, rk_method=0,
            rtol=1.0e-3, atol=1.0e-6,
            )

    assert success
    assert type(message) == str
    assert type(ts) == np.ndarray
    assert ts.dtype in (float, np.float64)
    assert type(ys) == np.ndarray
    assert ys.dtype in (float, np.float64)
    np.testing.assert_almost_equal(ts[0], t_span[0])
    np.testing.assert_almost_equal(ts[-1], t_span[-1])
    assert ys.shape[0] == 3
    assert ys.shape[1] == ts.shape[0]

    # Compare to scipy's solution
    solution = solve_ivp(diffeq, t_span, initial_y, method='RK23', t_eval=None, rtol=1.0e-3, atol=1.0e-6)
    scipy_y = solution.y
    scipy_t = solution.t

    assert scipy_t.shape == ts.shape
    assert scipy_y.shape == ys.shape
    assert np.allclose(ys, scipy_y, rtol=0.1)


def test_rk45_Yis2D_noTeval():
    """ Test the RK order 5 integrator for a 2D diffeq with no explicit Teval. """

    @njit(cacheable=False)
    def diffeq(t, y):
        y0 = y[0]
        y1 = y[1]
        y2 = y[2]
        dy0 = 0.5 * y0 - np.sin(2. * t)
        dy1 = 0.5 * y1 - np.sin(2. * t) - 0.005 * y0
        dy2 = 0.5 * y2 - np.sin(2. * t) + 0.005 * y1
        return np.asarray([dy0, dy1, dy2])

    initial_y = np.asarray([1., (1. - 0.005), (1. - (1. - 0.005))])
    t_span = (0., 100.)
    ts, ys, success, message = \
        rk_integrate(
            diffeq, t_span, initial_y, rk_method=1,
            rtol=1.0e-3, atol=1.0e-6
            )

    assert success
    assert type(message) == str
    assert type(ts) == np.ndarray
    assert ts.dtype in (float, np.float64)
    assert type(ys) == np.ndarray
    assert ys.dtype in (float, np.float64)
    np.testing.assert_almost_equal(ts[0], t_span[0])
    np.testing.assert_almost_equal(ts[-1], t_span[-1])
    assert ys.shape[0] == 3
    assert ys.shape[1] == ts.shape[0]

    # Compare to scipy's solution
    solution = solve_ivp(diffeq, t_span, initial_y, method='RK45', t_eval=None, rtol=1.0e-3, atol=1.0e-6)
    scipy_y = solution.y
    scipy_t = solution.t

    assert scipy_t.shape == ts.shape
    assert scipy_y.shape == ys.shape
    assert np.allclose(ys, scipy_y, rtol=0.1)


def test_rk45_Yis2D_noTeval_complex():
    """ Test the RK order 5 integrator for a 2D diffeq with no explicit Teval.
    y is in the complex domain """

    @njit(cacheable=False)
    def diffeq(t, y):
        y0 = y[0]
        y1 = y[1]
        y2 = y[2]
        dy0 = 0.5 * y0 - np.sin(2. * t)
        dy1 = 0.5 * y1 - np.sin(2. * t) - 0.005 * np.imag(y0)
        dy2 = 0.5 * y2 - np.sin(2. * t) + 0.005 * np.imag(y1)
        return np.asarray([dy0, dy1, dy2])

    initial_y = np.asarray([1. + 1.j, (1. - 0.005) * (1. + 1.j), (1. - (1. - 0.005)) * (1. + 1.j)])
    assert initial_y.dtype == np.complex128

    t_span = (0., 100.)
    ts, ys, success, message = \
        rk_integrate(
            diffeq, t_span, initial_y, rk_method=1,
            rtol=1.0e-3, atol=1.0e-6,
            )

    assert success
    assert type(message) == str
    assert type(ts) == np.ndarray
    assert ts.dtype in (float, np.float64)
    assert type(ys) == np.ndarray
    assert ys.dtype in (complex, np.complex128)
    np.testing.assert_almost_equal(ts[0], t_span[0])
    np.testing.assert_almost_equal(ts[-1], t_span[-1])
    assert ys.shape[0] == 3
    assert ys.shape[1] == ts.shape[0]

    # Compare to scipy's solution
    solution = solve_ivp(diffeq, t_span, initial_y, method='RK45', t_eval=None, rtol=1.0e-3, atol=1.0e-6)
    scipy_y = solution.y
    scipy_t = solution.t

    assert scipy_t.shape == ts.shape
    assert scipy_y.shape == ys.shape
    assert np.allclose(ys, scipy_y, rtol=0.1)


def test_rk45_Yis2D_withTeval_smallerN():
    """ Test the RK order 5 integrator for a 2D diffeq with Teval. Where Teval is smaller than the solution domain """

    @njit(cacheable=False)
    def diffeq(t, z):
        a, b, c, d = 1.5, 1., 3., 1.
        x = z[0]
        y = z[1]
        return np.asarray([a * x - b * x * y, -c * y + d * x * y])

    initial_y = np.asarray([10., 5.])
    t_span = (0., 15.)
    t_eval = np.linspace(0., 15., 25)
    ts, ys, success, message = \
        rk_integrate(
            diffeq, t_span, initial_y, rk_method=1, t_eval=t_eval,
            rtol=1.0e-3, atol=1.0e-6,
            )

    assert success
    assert type(message) == str
    assert type(ts) == np.ndarray
    assert ts.dtype in (float, np.float64)
    assert type(ys) == np.ndarray
    assert ys.dtype in (float, np.float64)
    np.testing.assert_almost_equal(ts[0], t_span[0])
    np.testing.assert_almost_equal(ts[-1], t_span[-1])
    assert ys.shape[0] == 2
    assert ys.shape[1] == ts.shape[0]

    # Compare to scipy's solution
    solution = solve_ivp(diffeq, t_span, initial_y, method='RK45', t_eval=t_eval, rtol=1.0e-3, atol=1.0e-6)
    scipy_y = solution.y
    scipy_t = solution.t

    assert scipy_t.shape == ts.shape
    assert scipy_y.shape == ys.shape


def test_rk45_Yis2D_withTeval_largerN():
    """ Test the RK order 5 integrator for a 2D diffeq with Teval. Where Teval is larger than the solution domain """

    @njit(cacheable=False)
    def diffeq(t, z):
        a, b, c, d = 1.5, 1., 3., 1.
        x = z[0]
        y = z[1]
        return np.asarray([a * x - b * x * y, -c * y + d * x * y])

    initial_y = np.asarray([10., 5.])
    t_span = (0., 15.)
    t_eval = np.linspace(0., 15., 80)
    ts, ys, success, message = \
        rk_integrate(
            diffeq, t_span, initial_y, rk_method=1, t_eval=t_eval,
            rtol=1.0e-3, atol=1.0e-6,
            )

    assert success
    assert type(message) == str
    assert type(ts) == np.ndarray
    assert ts.dtype in (float, np.float64)
    assert type(ys) == np.ndarray
    assert ys.dtype in (float, np.float64)
    np.testing.assert_almost_equal(ts[0], t_span[0])
    np.testing.assert_almost_equal(ts[-1], t_span[-1])
    assert ys.shape[0] == 2
    assert ys.shape[1] == ts.shape[0]

    # Compare to scipy's solution
    solution = solve_ivp(diffeq, t_span, initial_y, method='RK45', t_eval=t_eval, rtol=1.0e-3, atol=1.0e-6)
    scipy_y = solution.y
    scipy_t = solution.t

    assert scipy_t.shape == ts.shape
    assert scipy_y.shape == ys.shape
