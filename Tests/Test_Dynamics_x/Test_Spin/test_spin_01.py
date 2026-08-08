"""Spin-dynamics calculator (``TidalPy.dynamics_x.Spin``).

Checks the moment-of-inertia formula (solid sphere, shell, structure factor), the synchronous spin,
and the tidal spin-rate change ``dspin/dt = M_host dU_dO / I``.
"""
import math

import numpy as np
import pytest

from TidalPy.dynamics_x import Spin


def test_moment_of_inertia_solid_sphere():
    """Solid uniform sphere: I = (2/5) M R^2."""
    spin = Spin()
    mass, radius = 5.0e24, 6.0e6
    assert math.isclose(spin.calc_moment_of_inertia(mass, radius),
                        0.4 * mass * radius ** 2, rel_tol=1e-14)


def test_moment_of_inertia_shell():
    """Uniform shell: I = (2/5) M (R_o^5 - R_i^5)/(R_o^3 - R_i^3)."""
    spin = Spin()
    mass, r_o, r_i = 3.0e23, 2.0e6, 1.0e6
    expected = 0.4 * mass * (r_o ** 5 - r_i ** 5) / (r_o ** 3 - r_i ** 3)
    assert math.isclose(spin.calc_moment_of_inertia(mass, r_o, r_i), expected, rel_tol=1e-14)


def test_moment_of_inertia_factor_scales():
    mass, radius, factor = 1.0e24, 1.0e6, 0.33
    spin = Spin(moment_of_inertia_factor=factor)
    assert spin.moment_of_inertia_factor == factor
    assert math.isclose(spin.calc_moment_of_inertia(mass, radius),
                        factor * 0.4 * mass * radius ** 2, rel_tol=1e-14)


def test_degenerate_shell_is_nan():
    spin = Spin()
    assert np.isnan(spin.calc_moment_of_inertia(1.0e24, 1.0e6, 1.0e6))  # zero thickness


def test_synchronous_spin_equals_mean_motion():
    spin = Spin()
    n = 2.0 * np.pi / 86400.0
    assert spin.calc_synchronous_spin(n) == n


def test_dspin_dt_formula_and_sign():
    """dspin/dt = M_host dU_dO / I; scales linearly and flips sign with dU_dO."""
    spin = Spin()
    mass, radius = 5.0e24, 6.0e6
    moi = spin.calc_moment_of_inertia(mass, radius)
    host = 1.9e27
    dU_dO = 3.0e-8
    expected = host * dU_dO / moi
    assert math.isclose(spin.calc_dspin_dt(host, dU_dO, moi), expected, rel_tol=1e-14)
    # opposite sign of dU_dO -> opposite dspin/dt
    assert math.isclose(spin.calc_dspin_dt(host, -dU_dO, moi), -expected, rel_tol=1e-14)


def test_dspin_dt_zero_moi_is_nan():
    spin = Spin()
    assert np.isnan(spin.calc_dspin_dt(1.9e27, 3.0e-8, 0.0))
