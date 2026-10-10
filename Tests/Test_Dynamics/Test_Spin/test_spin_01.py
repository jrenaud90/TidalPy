"""The spin-dynamics calculator ``Spin``: moment of inertia, synchronous spin, and dspin/dt."""
import math

import numpy as np
import pytest

from TidalPy.Dynamics import Spin


def test_default_factor_is_uniform_sphere():
    """The default factor is 0.4, so I = (2/5) M R^2."""
    spin = Spin()
    mass, radius = 5.0e24, 6.0e6
    assert spin.moment_of_inertia_factor == 0.4
    assert math.isclose(spin.calc_moment_of_inertia(mass, radius), 0.4 * mass * radius ** 2, rel_tol=1e-14)


@pytest.mark.parametrize("factor", [0.25, 0.3307, 0.4, 2.0 / 3.0])
def test_moment_of_inertia_uses_conventional_factor(factor):
    """I = factor M R^2 for any physical factor, up to the thin-shell 2/3."""
    mass, radius = 1.0e24, 1.0e6
    spin = Spin(moment_of_inertia_factor=factor)
    assert spin.moment_of_inertia_factor == factor
    assert math.isclose(spin.calc_moment_of_inertia(mass, radius), factor * mass * radius ** 2, rel_tol=1e-14)


# 1.0 is the old ratio-to-a-uniform-sphere convention, which must now be rejected.
@pytest.mark.parametrize("factor", [0.0, -0.4, 0.7, 1.0, math.nan, math.inf])
def test_unphysical_factor_raises(factor):
    """A factor outside (0, 2/3] raises."""
    with pytest.raises(ValueError):
        Spin(moment_of_inertia_factor=factor)


def test_synchronous_spin_equals_mean_motion():
    orbital_frequency = 2.0 * np.pi / 86400.0
    assert Spin().calc_synchronous_spin(orbital_frequency) == orbital_frequency


@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_dspin_dt_formula_and_sign(sign):
    """dspin/dt = M_host dU/dO / I, which flips sign with dU/dO."""
    spin = Spin()
    moment_of_inertia = spin.calc_moment_of_inertia(5.0e24, 6.0e6)
    host_mass = 1.9e27
    dU_dO = sign * 3.0e-8
    expected = host_mass * dU_dO / moment_of_inertia
    assert math.isclose(spin.calc_dspin_dt(host_mass, dU_dO, moment_of_inertia), expected, rel_tol=1e-14)


def test_dspin_dt_zero_moi_is_nan():
    assert np.isnan(Spin().calc_dspin_dt(1.9e27, 3.0e-8, 0.0))
