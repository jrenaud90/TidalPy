"""On-demand 3D tidal kernel (Tides_x.multilayer) point evaluation - validated against the legacy.

The new on-demand 3D system evaluates the tidal potential, the strain/stress tensors, and the
volumetric heating at a single point (no materialized grids). These tests check that the compiled
C++ kernel reproduces the legacy ``calculate_strain_stress`` + ``calculate_volumetric_heating`` and
the legacy ``synchronous_low_e`` tidal potential to machine precision (solid compressible - the case
the legacy supports).
"""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.tides.potential.synchronous_low_e import tidal_potential as legacy_potential
from TidalPy.tides.multilayer.stress_strain import calculate_strain_stress
from TidalPy.tides.heating import calculate_volumetric_heating

from TidalPy.Tides_x.multilayer.stress_strain import strain_stress_heating_point
from TidalPy.Tides_x.potential.tidal_potential import SyncLowEPotential


def _legacy_strain_stress_heat(y, shear, bulk, radius, l, pot6, theta):
    tp, tpt, tpp, tpt2, tpp2, tptp = (np.array([[[v]]], dtype=np.float64) for v in pot6)
    strains, stresses = calculate_strain_stress(
        tp, tpt, tpp, tpt2, tpp2, tptp,
        np.asarray(y, dtype=np.complex128).reshape(6, 1),
        np.array([0.0]), np.array([theta]), np.array([0.0]),
        np.array([radius]),
        np.array([shear], dtype=np.complex128), np.array([bulk], dtype=np.complex128),
        1.0, l)
    e = strains[:, 0, 0, 0, 0]
    s = stresses[:, 0, 0, 0, 0]
    vh = calculate_volumetric_heating(s.reshape(6, 1, 1, 1, 1), e.reshape(6, 1, 1, 1, 1))[0, 0, 0, 0]
    return e, s, vh


def test_potential_sync_low_e_matches_legacy():
    sync_potential = SyncLowEPotential()
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(150):
        radius = rng.uniform(1e5, 6e6); lon = rng.uniform(0, 2 * np.pi)
        colat = rng.uniform(0.15, np.pi - 0.15); t = rng.uniform(0, 1e5)
        n = 2 * np.pi / (86400.0 * rng.uniform(0.5, 5.0)); ecc = rng.uniform(0.001, 0.2)
        hm = 1e26 * rng.uniform(0.5, 5.0); a = rng.uniform(3e8, 3e9)
        _, pots = sync_potential.calc_modes(n, 0.0, ecc, hm, a, radius, colat, lon, t)
        mine = np.asarray(pots[0])
        _, _, ptup = legacy_potential(radius, lon, colat, t, n, ecc, hm, a)
        ref = np.array(ptup['n'])
        worst = max(worst, float(np.max(np.abs(mine - ref) / (np.abs(ref) + 1e-30))))
    assert worst < 1e-12, f"potential differs from legacy: {worst:.3e}"


def test_strain_stress_heating_matches_legacy():
    rng = np.random.default_rng(99)
    worst_e = worst_s = worst_h = 0.0
    for _ in range(300):
        y = rng.standard_normal(6) + 1j * rng.standard_normal(6)
        shear = complex(5e10 + rng.standard_normal() * 1e9, rng.uniform(1e7, 1e9))
        bulk = complex(1.3e11 + rng.standard_normal() * 1e9, rng.uniform(1e6, 1e8))
        radius = rng.uniform(1e5, 6e6)
        pot6 = tuple(rng.standard_normal(6) * 10.0)
        theta = rng.uniform(0.15, np.pi - 0.15)

        e_ref, s_ref, h_ref = _legacy_strain_stress_heat(y, shear, bulk, radius, 2.0, pot6, theta)
        e_c, s_c, h_c = strain_stress_heating_point(
            np.ascontiguousarray(y, dtype=np.complex128), shear, bulk, radius, 2.0, True, False,
            pot6, theta)
        worst_e = max(worst_e, float(np.max(np.abs(e_c - e_ref) / (np.abs(e_ref) + 1e-30))))
        worst_s = max(worst_s, float(np.max(np.abs(s_c - s_ref) / (np.abs(s_ref) + 1e-30))))
        worst_h = max(worst_h, abs(h_c - h_ref) / (abs(h_ref) + 1e-30))
    assert worst_e < 1e-11, f"strain differs: {worst_e:.3e}"
    assert worst_s < 1e-11, f"stress differs: {worst_s:.3e}"
    assert worst_h < 1e-10, f"heating differs: {worst_h:.3e}"


def test_liquid_returns_nan():
    """The shear kernel is solid-only; a liquid (shear=0) point returns NaN strains."""
    y = np.ones(6, dtype=np.complex128)
    e, s, h = strain_stress_heating_point(y, 0.0 + 0j, 1.3e11 + 0j, 1e6, 2.0, False, False,
                                          (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0)
    assert np.all(np.isnan(e)) and np.isnan(h)
