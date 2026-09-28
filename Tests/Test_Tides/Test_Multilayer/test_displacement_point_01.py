"""The point-wise tidal displacement kernel ``displacement_point``."""
from pathlib import Path

import numpy as np
import pytest

from TidalPy.Tides.multilayer.stress_strain import displacement_point

# Classic `calculate_displacements` at 200 inputs drawn with default_rng(7), from the 0.8.0 snapshot 8b8e0b12.
FROZEN_PATH = Path(__file__).parent / "frozen" / "test_displacement_point_01.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    LEGACY_REFERENCE = {key: frozen_file[key] for key in frozen_file.files}


def test_displacement_point_closed_form():
    """u = (y1 U, y3 dU/dtheta, y3 dU/dphi / sin(theta))."""
    y = np.array([1.5 - 0.2j, 3.0, 0.7 + 0.1j, 2.0, 0.5, 0.1], dtype=np.complex128)
    pot6 = (10.0, -3.0, 4.0, 1.0, 2.0, 0.5)
    theta = 1.1
    u = displacement_point(y, pot6, theta)
    assert u.shape == (3,) and u.dtype == np.complex128
    assert u[0] == pytest.approx(y[0] * pot6[0])
    assert u[1] == pytest.approx(y[2] * pot6[1])
    assert u[2] == pytest.approx(y[2] * pot6[2] / np.sin(theta))


def test_displacement_point_matches_legacy():
    """The kernel matches the frozen classic displacements."""
    worst = 0.0
    for sample_i in range(LEGACY_REFERENCE["colatitudes"].size):
        y = LEGACY_REFERENCE["radial_functions"][sample_i]
        pot6 = tuple(LEGACY_REFERENCE["potential_terms"][sample_i])
        theta = float(LEGACY_REFERENCE["colatitudes"][sample_i])
        reference = LEGACY_REFERENCE["classic_displacements"][sample_i]
        u = displacement_point(np.ascontiguousarray(y), pot6, theta)
        worst = max(worst, float(np.max(np.abs(u - reference) / (np.abs(reference) + 1e-30))))
    assert worst < 1e-13


def test_displacement_point_short_y_and_poles():
    """Three radial functions suffice; the azimuthal term is zero or NaN at the pole."""
    y = np.array([2.0 + 0j, 0.0, 0.5 + 0j], dtype=np.complex128)
    zonal = (1.0, 0.5, 0.0, 0.0, 0.0, 0.0)
    u = displacement_point(y, zonal, 0.0)
    assert u[0] == 2.0 and u[1] == 0.25 and u[2] == 0.0
    sectoral = (1.0, 0.5, 0.3, 0.0, 0.0, 0.0)
    assert np.isnan(displacement_point(y, sectoral, 0.0)[2])  # 1/sin(theta) is singular for m != 0 at the pole
    assert displacement_point(y, sectoral, np.pi / 2)[2] == pytest.approx(0.15)
