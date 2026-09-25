"""Spin-dynamics calculator (``TidalPy.dynamics_x.Spin``).

Checks the moment of inertia from the conventional factor ``C / (M R^2)``, the factor validation, the synchronous
spin, and the tidal spin-rate change ``dspin/dt = M_host dU_dO / I``.
"""
import math

import numpy as np
import pytest

from TidalPy.dynamics_x import Spin


def test_default_factor_is_uniform_sphere():
    """The default factor is the uniform-sphere value, so I = (2/5) M R^2."""
    spin = Spin()
    mass, radius = 5.0e24, 6.0e6
    assert spin.moment_of_inertia_factor == 0.4
    assert math.isclose(spin.calc_moment_of_inertia(mass, radius), 0.4 * mass * radius ** 2, rel_tol=1e-14)


@pytest.mark.parametrize("factor", [0.25, 0.3307, 0.4, 2.0 / 3.0])
def test_moment_of_inertia_uses_conventional_factor(factor):
    """I = factor * M R^2 for any physical factor, including the thin-shell limit 2/3."""
    mass, radius = 1.0e24, 1.0e6
    spin = Spin(moment_of_inertia_factor=factor)
    assert spin.moment_of_inertia_factor == factor
    assert math.isclose(spin.calc_moment_of_inertia(mass, radius), factor * mass * radius ** 2, rel_tol=1e-14)


@pytest.mark.parametrize("factor", [0.0, -0.4, 0.7, 1.0, math.nan, math.inf])
def test_unphysical_factor_raises(factor):
    """A factor outside (0, 2/3] is rejected; 1.0 is the old ratio-to-a-uniform-sphere convention."""
    with pytest.raises(ValueError):
        Spin(moment_of_inertia_factor=factor)


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
