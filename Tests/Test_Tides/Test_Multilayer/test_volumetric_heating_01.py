"""The point volumetric heating ``volumetric_heating``: (|frequency| / 2) sum_k w_k Im(sigma_k conj(eps_k)), signed."""
import math

import numpy as np
import pytest

from TidalPy.Tides.multilayer.stress_strain import volumetric_heating

# Off-diagonal components count twice in the symmetric tensor.
_WEIGHTS = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])
# At 2 rad s-1 the |frequency| / 2 factor is one, leaving the bilinear form.
_UNIT_FACTOR_FREQUENCY = 2.0


def _reference(stress, strain):
    return np.sum(_WEIGHTS * (stress.imag * strain.real - stress.real * strain.imag))


@pytest.mark.parametrize("component, expected", [(0, 2.0), (3, 4.0)], ids=["diagonal", "off_diagonal"])
def test_single_component(component, expected):
    """sigma = 1 + 2i with eps = 3 + 4i on one component gives 2 on the diagonal and 4 off it."""
    stress = np.zeros(6, dtype=np.complex128)
    strain = np.zeros(6, dtype=np.complex128)
    stress[component] = 1.0 + 2.0j
    strain[component] = 3.0 + 4.0j
    assert math.isclose(volumetric_heating(stress, strain, _UNIT_FACTOR_FREQUENCY), expected, rel_tol=1e-12)


def test_in_phase_response_dissipates_nothing():
    """Real stress and strain carry no phase lag, so nothing dissipates."""
    rng = np.random.default_rng(7)
    stress = rng.normal(size=6).astype(np.complex128)
    strain = rng.normal(size=6).astype(np.complex128)
    assert volumetric_heating(stress, strain, _UNIT_FACTOR_FREQUENCY) == 0.0


@pytest.mark.parametrize('seed', (1, 42, 1234))
def test_matches_reference_formula(seed):
    """Random complex tensors match an independent numpy evaluation."""
    rng = np.random.default_rng(seed)
    stress = (rng.normal(size=6) + 1j * rng.normal(size=6)).astype(np.complex128)
    strain = (rng.normal(size=6) + 1j * rng.normal(size=6)).astype(np.complex128)
    assert math.isclose(
        volumetric_heating(stress, strain, _UNIT_FACTOR_FREQUENCY), _reference(stress, strain), rel_tol=1e-12)


def test_the_frequency_factor_is_half_its_magnitude():
    """The factor is |frequency| / 2, so a negative frequency heats as its magnitude does."""
    rng = np.random.default_rng(3)
    stress = (rng.normal(size=6) + 1j * rng.normal(size=6)).astype(np.complex128)
    strain = (rng.normal(size=6) + 1j * rng.normal(size=6)).astype(np.complex128)
    form = _reference(stress, strain)
    assert math.isclose(volumetric_heating(stress, strain, 1.0e-5), 0.5e-5 * form, rel_tol=1e-12)
    assert volumetric_heating(stress, strain, -1.0e-5) == volumetric_heating(stress, strain, 1.0e-5)
    assert volumetric_heating(stress, strain, 0.0) == 0.0


def test_a_non_dissipative_response_reads_negative():
    """A response with the wrong phase (Im(mu) < 0) gives negative heating rather than its magnitude, as the world's
    totals do."""
    stress = np.zeros(6, dtype=np.complex128)
    strain = np.zeros(6, dtype=np.complex128)
    stress[0] = 1.0 - 2.0j
    strain[0] = 3.0
    assert math.isclose(volumetric_heating(stress, strain, _UNIT_FACTOR_FREQUENCY), -6.0, rel_tol=1e-12)
