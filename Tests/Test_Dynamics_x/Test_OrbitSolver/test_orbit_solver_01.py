"""Orbital rate engine (``TidalPy.dynamics_x.OrbitSolver``).

Checks the da/dt, de/dt, dn/dt formulas, the reduced-mass conversion, and the circular-orbit guard.
This is the low-level rate engine; the System class attaches it and pulls the orbital
state / tidal-potential derivatives from its worlds. End-to-end energy conservation with the world
spin is validated in ``Tests/Test_Structures_x/Test_Worlds/test_world_spin_01.py``.
"""
import math

import numpy as np
import pytest

from TidalPy.dynamics_x import OrbitSolver


_N = 2.0e-5
_A = 1.0e9
_E = 0.1
_MT = 6.0e24
_MH = 1.9e27
_DUDM = 3.0e-9
_DUDW = -1.2e-9


def _dR(dU):
    return -((_MT + _MH) / _MT) * dU


def test_da_dt_formula():
    orbit = OrbitSolver()
    expected = 2.0 / (_N * _A) * _dR(_DUDM)
    assert math.isclose(orbit.calc_da_dt(_N, _A, _E, _MT, _MH, _DUDM), expected, rel_tol=1e-14)


def test_de_dt_formula():
    orbit = OrbitSolver()
    e_term = math.sqrt(1.0 - _E ** 2)
    denom = _N * _A ** 2 * _E
    expected = (e_term / denom) * (e_term * _dR(_DUDM) - _dR(_DUDW))
    assert math.isclose(orbit.calc_de_dt(_N, _A, _E, _MT, _MH, _DUDM, _DUDW), expected, rel_tol=1e-14)


def test_dn_dt_from_kepler():
    orbit = OrbitSolver()
    da_dt = orbit.calc_da_dt(_N, _A, _E, _MT, _MH, _DUDM)
    assert math.isclose(orbit.calc_dn_dt(_N, _A, da_dt), -1.5 * (_N / _A) * da_dt, rel_tol=1e-14)


def test_circular_orbit_de_dt_is_zero():
    orbit = OrbitSolver()
    assert orbit.calc_de_dt(_N, _A, 0.0, _MT, _MH, _DUDM, _DUDW) == 0.0


def test_calc_derivatives_matches_individual():
    orbit = OrbitSolver()
    out = orbit.calc_derivatives(_N, _A, _E, _MT, _MH, _DUDM, _DUDW)
    assert math.isclose(out["da_dt"], orbit.calc_da_dt(_N, _A, _E, _MT, _MH, _DUDM), rel_tol=1e-14)
    assert math.isclose(out["de_dt"], orbit.calc_de_dt(_N, _A, _E, _MT, _MH, _DUDM, _DUDW), rel_tol=1e-14)
    assert math.isclose(out["dn_dt"], orbit.calc_dn_dt(_N, _A, out["da_dt"]), rel_tol=1e-14)


def test_degenerate_returns_nan():
    orbit = OrbitSolver()
    assert np.isnan(orbit.calc_da_dt(0.0, _A, _E, _MT, _MH, _DUDM))   # n = 0
    assert np.isnan(orbit.calc_da_dt(_N, _A, _E, 0.0, _MH, _DUDM))    # target mass = 0
