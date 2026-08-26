"""Tests for the point volumetric heating binding (TidalPy.Tides_x.multilayer.stress_strain).

``volumetric_heating(stress, strain)`` takes the 6 complex stress and strain tensor components
(3 diagonal then 3 off-diagonal) and returns |sum_k w_k Im(sigma_k conj(eps_k))| with w = 1 on the
diagonals and 2 on the off-diagonals.
"""
import math

import numpy as np
import pytest

from TidalPy.Tides_x.multilayer.stress_strain import volumetric_heating

_WEIGHTS = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])


def _reference(stress, strain):
    return abs(np.sum(_WEIGHTS * (stress.imag * strain.real - stress.real * strain.imag)))


def test_single_diagonal_component():
    """One diagonal component: h = Im(sigma) Re(eps) - Re(sigma) Im(eps) = 2*3 - 1*4 = 2."""
    stress = np.zeros(6, dtype=np.complex128)
    strain = np.zeros(6, dtype=np.complex128)
    stress[0] = 1.0 + 2.0j
    strain[0] = 3.0 + 4.0j
    assert math.isclose(volumetric_heating(stress, strain), 2.0, rel_tol=1e-12)


def test_off_diagonal_component_is_doubled():
    """The same values on an off-diagonal component count twice (symmetric tensor)."""
    stress = np.zeros(6, dtype=np.complex128)
    strain = np.zeros(6, dtype=np.complex128)
    stress[3] = 1.0 + 2.0j
    strain[3] = 3.0 + 4.0j
    assert math.isclose(volumetric_heating(stress, strain), 4.0, rel_tol=1e-12)


def test_in_phase_response_dissipates_nothing():
    """Real (in phase) stress and strain carry no phase lag, so the dissipation is zero."""
    rng = np.random.default_rng(7)
    stress = rng.normal(size=6).astype(np.complex128)
    strain = rng.normal(size=6).astype(np.complex128)
    assert volumetric_heating(stress, strain) == 0.0


@pytest.mark.parametrize('seed', (1, 42, 1234))
def test_matches_reference_formula(seed):
    """Random complex tensors match an independent numpy evaluation of the same definition."""
    rng = np.random.default_rng(seed)
    stress = (rng.normal(size=6) + 1j * rng.normal(size=6)).astype(np.complex128)
    strain = (rng.normal(size=6) + 1j * rng.normal(size=6)).astype(np.complex128)
    assert math.isclose(volumetric_heating(stress, strain), _reference(stress, strain), rel_tol=1e-12)
