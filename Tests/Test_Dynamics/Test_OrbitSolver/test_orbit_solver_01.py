"""The orbital rate engine ``OrbitSolver``: da/dt, de/dt, dn/dt, and degenerate inputs."""
import math

import numpy as np
import pytest

from TidalPy.Dynamics import OrbitSolver

_N = 2.0e-5
_A = 1.0e9
_E = 0.1
_MT = 6.0e24
_MH = 1.9e27
_DUDM = 3.0e-9
_DUDW = -1.2e-9
# (n, a, e, target mass, host mass), the leading arguments of every rate.
_STATE = (_N, _A, _E, _MT, _MH)


def _dR(dU):
    """The disturbing-function derivative with the reduced-mass conversion."""
    return -((_MT + _MH) / _MT) * dU


def test_da_dt_formula():
    """da/dt = 2 / (n a) dR/dM."""
    expected = 2.0 / (_N * _A) * _dR(_DUDM)
    assert math.isclose(OrbitSolver().calc_da_dt(*_STATE, _DUDM), expected, rel_tol=1e-14)


def test_de_dt_formula():
    """de/dt = sqrt(1 - e^2) / (n a^2 e) (sqrt(1 - e^2) dR/dM - dR/dw)."""
    e_term = math.sqrt(1.0 - _E ** 2)
    expected = (e_term / (_N * _A ** 2 * _E)) * (e_term * _dR(_DUDM) - _dR(_DUDW))
    assert math.isclose(OrbitSolver().calc_de_dt(*_STATE, _DUDM, _DUDW), expected, rel_tol=1e-14)


def test_dn_dt_from_kepler():
    """dn/dt = -(3/2) (n / a) da/dt."""
    orbit = OrbitSolver()
    da_dt = orbit.calc_da_dt(*_STATE, _DUDM)
    assert math.isclose(orbit.calc_dn_dt(_N, _A, da_dt), -1.5 * (_N / _A) * da_dt, rel_tol=1e-14)


def test_circular_orbit_de_dt_is_zero():
    circular_state = (_N, _A, 0.0, _MT, _MH)
    assert OrbitSolver().calc_de_dt(*circular_state, _DUDM, _DUDW) == 0.0


def test_calc_derivatives_matches_individual():
    """calc_derivatives returns the individual rates."""
    orbit = OrbitSolver()
    out = orbit.calc_derivatives(*_STATE, _DUDM, _DUDW)
    assert math.isclose(out["da_dt"], orbit.calc_da_dt(*_STATE, _DUDM), rel_tol=1e-14)
    assert math.isclose(out["de_dt"], orbit.calc_de_dt(*_STATE, _DUDM, _DUDW), rel_tol=1e-14)
    assert math.isclose(out["dn_dt"], orbit.calc_dn_dt(_N, _A, out["da_dt"]), rel_tol=1e-14)


@pytest.mark.parametrize("state", [(0.0, _A, _E, _MT, _MH), (_N, _A, _E, 0.0, _MH)],
                         ids=["zero_mean_motion", "zero_target_mass"])
def test_degenerate_returns_nan(state):
    assert np.isnan(OrbitSolver().calc_da_dt(*state, _DUDM))
