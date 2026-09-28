"""Propagation matrix Love numbers against the closed form, the frozen classic matrix method, and shooting.

Homogeneous static incompressible sphere closed form (complex mu via the correspondence principle):
    k_l = (3 / (2 (l - 1))) / (1 + mu_eff),  h_l = ((2 l + 1) / (2 (l - 1))) / (1 + mu_eff),
    mu_eff = ((2 l^2 + 4 l + 3) / l) * mu / (rho g R).
"""
from pathlib import Path

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Rheology import Maxwell
from TidalPy.RadialSolver.solver import radial_solver

# Classic propagation-matrix Love numbers, shape (1, 3) complex: k, h, l.
FROZEN_PATH = Path(__file__).parent / "frozen" / "test_c_matrix_parity.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_MATRIX_LOVE = {key: frozen_file[key] for key in frozen_file.files}

# The large bulk modulus plays no role once is_incompressible is set.
frequency = 1.0 / (86400. * 1.5)
N = 100
radius_array = np.linspace(0.0, 6000.e3, N)
density = 5500.
density_array = density * np.ones_like(radius_array)
bulk_modulus_array = 1.0e14 * np.ones(N, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)
complex_shear_modulus_array = Maxwell().calc_complex_modulus_vectorize_modulus(shear_array, viscosity_array, frequency)
planet_radius = radius_array[-1]
upper_radius_by_layer = np.asarray((planet_radius,))
surface_gravity = 4.0 * np.pi * G * density * planet_radius / 3.0
complex_shear = complex_shear_modulus_array[0]


def analytic_love(degree_l):
    """Closed-form homogeneous incompressible static k and h."""
    mu_eff = (2.0 * degree_l**2 + 4.0 * degree_l + 3.0) / degree_l * complex_shear \
        / (density * surface_gravity * planet_radius)
    k_l = (3.0 / (2.0 * (degree_l - 1.0))) / (1.0 + mu_eff)
    h_l = (2.0 * degree_l + 1.0) / (2.0 * (degree_l - 1.0)) / (1.0 + mu_eff)
    return k_l, h_l


def _solve(degree_l, is_static, **kwargs):
    return radial_solver(
        radius_array,
        density_array,
        bulk_modulus_array,
        complex_shear_modulus_array,
        frequency,
        density,
        ('solid',),
        (is_static,),
        (True,),
        upper_radius_by_layer,
        degree_l=degree_l,
        solve_for=('tidal',),
        verbose=False,
        nondimensionalize=True,
        raise_on_fail=True,
        warnings=False,
        **kwargs)


def run_matrix(degree_l, core_model):
    return _solve(degree_l, True, core_model=core_model, love_method='propagation_matrix')


@pytest.mark.parametrize('core_model', (0, 1, 2, 3))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_matrix_love_matches_analytic_homogeneous(degree_l, core_model):
    """The matrix solver reproduces the closed-form k and h, and h / k = (2l + 1) / 3."""
    k_ref, h_ref = analytic_love(degree_l)
    out = run_matrix(degree_l, core_model)
    assert out.success
    k, h, _ = out.love[0]

    # Core model 0 is consistent with the homogeneous solution; the alternate core starts perturb it to about 1e-5.
    rtol = 1.0e-12 if core_model == 0 else 1.0e-4
    np.testing.assert_allclose(k, k_ref, rtol=rtol)
    np.testing.assert_allclose(h, h_ref, rtol=rtol)
    np.testing.assert_allclose(h / k, (2.0 * degree_l + 1.0) / 3.0, rtol=1.0e-8)


@pytest.mark.parametrize('core_model', (0, 1, 2, 3))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_matrix_old_vs_new_parity(degree_l, core_model):
    """k, h, and l match the frozen classic matrix method."""
    out_new = run_matrix(degree_l, core_model)
    old_love = CLASSIC_MATRIX_LOVE[f"degree_l_{degree_l}__core_model_{core_model}"]
    assert out_new.success
    np.testing.assert_allclose(out_new.love, old_love, rtol=1.0e-9, atol=1.0e-10)


@pytest.mark.parametrize('degree_l', (2, 3))
def test_matrix_vs_shooting(degree_l):
    """The static matrix and dynamic shooting methods agree to the inertial-term level (1e-3)."""
    out_matrix = run_matrix(degree_l, core_model=0)
    out_shooting = _solve(
        degree_l,
        False,
        use_kamata=True,
        integration_method='DOP853',
        integration_rtol=1.0e-10,
        integration_atol=1.0e-12,
        scale_rtols_bylayer_type=False,
        love_method='radial_solver')
    assert out_matrix.success and out_shooting.success
    np.testing.assert_allclose(out_shooting.love, out_matrix.love, rtol=1.0e-3)
