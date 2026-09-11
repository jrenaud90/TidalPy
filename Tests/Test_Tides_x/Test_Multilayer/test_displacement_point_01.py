"""Tests for the point-wise tidal displacement kernel (TidalPy.Tides_x.multilayer.stress_strain.displacement_point).

NOTE (0.9.0): the comparison against the classic `calculate_displacements` loses its reference when the legacy
tree is removed; the closed-form checks stand on their own.
"""
import numpy as np
import pytest

from TidalPy.tides.multilayer.displacements import calculate_displacements
from TidalPy.Tides_x.multilayer.stress_strain import displacement_point


def _legacy(y, pot6, theta):
    """Classic displacement arrays for one radius and one grid point."""
    potential = np.array([[[pot6[0]]]], dtype=np.float64)
    d_theta = np.array([[[pot6[1]]]], dtype=np.float64)
    d_phi = np.array([[[pot6[2]]]], dtype=np.float64)
    radial, polar, azimuthal = calculate_displacements(
        potential, d_theta, d_phi, np.asarray(y, dtype=np.complex128).reshape(6, 1), np.array([[[theta]]]))
    return np.array([radial[0, 0, 0, 0], polar[0, 0, 0, 0], azimuthal[0, 0, 0, 0]])


def test_displacement_point_closed_form():
    y = np.array([1.5 - 0.2j, 3.0, 0.7 + 0.1j, 2.0, 0.5, 0.1], dtype=np.complex128)
    pot6 = (10.0, -3.0, 4.0, 1.0, 2.0, 0.5)
    theta = 1.1
    u = displacement_point(y, pot6, theta)
    assert u.shape == (3,) and u.dtype == np.complex128
    assert u[0] == pytest.approx(y[0] * pot6[0])
    assert u[1] == pytest.approx(y[2] * pot6[1])
    assert u[2] == pytest.approx(y[2] * pot6[2] / np.sin(theta))


def test_displacement_point_matches_legacy():
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(200):
        y = rng.standard_normal(6) + 1j * rng.standard_normal(6)
        pot6 = tuple(rng.standard_normal(6) * 10.0)
        theta = rng.uniform(0.1, np.pi - 0.1)
        reference = _legacy(y, pot6, theta)
        u = displacement_point(np.ascontiguousarray(y), pot6, theta)
        worst = max(worst, float(np.max(np.abs(u - reference) / (np.abs(reference) + 1e-30))))
    assert worst < 1e-13


def test_displacement_point_short_y_and_poles():
    y = np.array([2.0 + 0j, 0.0, 0.5 + 0j], dtype=np.complex128)   # y1..y3 only
    zonal = (1.0, 0.5, 0.0, 0.0, 0.0, 0.0)
    u = displacement_point(y, zonal, 0.0)
    assert u[0] == 2.0 and u[1] == 0.25 and u[2] == 0.0      # zonal potential: no azimuthal motion at the pole
    sectoral = (1.0, 0.5, 0.3, 0.0, 0.0, 0.0)
    assert np.isnan(displacement_point(y, sectoral, 0.0)[2])  # 1/sin(theta) singular for m != 0 at the pole
    assert displacement_point(y, sectoral, np.pi / 2)[2] == pytest.approx(0.15)
