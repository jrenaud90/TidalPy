"""Numerical parity tests for the propagation matrix method.

The propagation matrix method applies to solid, static, incompressible layers, for which a homogeneous
sphere has the classic closed-form Love numbers (with the correspondence principle extending them to a
complex shear modulus):

    k_l = (3 / (2 (l - 1))) / (1 + mu_eff),    h_l = ((2 l + 1) / (2 (l - 1))) / (1 + mu_eff),
    mu_eff = ((2 l^2 + 4 l + 3) / l) * mu / (rho g R),    h_l / k_l = (2 l + 1) / 3 exactly.

These tests pin the `_x` matrix solver against that analytic solution (core model 0 reproduces it to
machine precision; the alternate core starting conditions 1..3 agree to ~1e-5), against the original
solver's matrix method, and against the `_x` shooting method (dynamic incompressible, which differs
from static by the inertial term, ~2e-5 relative here).
"""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.rheology.models import Maxwell
from TidalPy.RadialSolver import radial_solver as radial_solver_old
from TidalPy.RadialSolver_x.solver import radial_solver as radial_solver_new

# Homogeneous solid planet; the large real bulk modulus plays no role once is_incompressible is set.
frequency = 1.0 / (86400. * 1.5)
N = 100
radius_array = np.linspace(0.0, 6000.e3, N)
density = 5500.
density_array = density * np.ones_like(radius_array)
bulk_modulus_array = 1.0e14 * np.ones(N, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)
complex_shear_modulus_array = np.empty(N, dtype=np.complex128)
Maxwell().vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear_modulus_array)
planet_radius = radius_array[-1]
upper_radius_by_layer = np.asarray((planet_radius,))
surface_gravity = 4.0 * np.pi * G * density * planet_radius / 3.0
complex_shear = complex_shear_modulus_array[0]


def analytic_love(degree_l):
    """Closed-form homogeneous incompressible static k and h at this planet's parameters."""
    mu_eff = (2.0 * degree_l**2 + 4.0 * degree_l + 3.0) / degree_l * complex_shear \
        / (density * surface_gravity * planet_radius)
    k_l = (3.0 / (2.0 * (degree_l - 1.0))) / (1.0 + mu_eff)
    h_l = (2.0 * degree_l + 1.0) / (2.0 * (degree_l - 1.0)) / (1.0 + mu_eff)
    return k_l, h_l


def run_matrix(solver_func, degree_l, core_model):
    return solver_func(
        radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
        frequency, density,
        ('solid',), (True,), (True,), upper_radius_by_layer,
        degree_l=degree_l, solve_for=('tidal',), core_model=core_model,
        use_prop_matrix=True, verbose=False, nondimensionalize=True,
        raise_on_fail=True, warnings=False)


@pytest.mark.parametrize('core_model', (0, 1, 2, 3))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_matrix_love_matches_analytic_homogeneous(degree_l, core_model):
    """The matrix solver reproduces the closed-form homogeneous incompressible Love numbers.

    Core model 0 uses the solid starting condition consistent with the homogeneous solution and lands on
    it to machine precision; the alternate core models (1..3) perturb the start and agree to ~1e-5.
    """
    k_ref, h_ref = analytic_love(degree_l)
    out = run_matrix(radial_solver_new, degree_l, core_model)
    assert out.success
    k, h, shida = out.love[0]

    rtol = 1.0e-12 if core_model == 0 else 1.0e-4
    np.testing.assert_allclose(k, k_ref, rtol=rtol)
    np.testing.assert_allclose(h, h_ref, rtol=rtol)
    # The h/k ratio is exactly (2l + 1)/3 for a homogeneous incompressible sphere, for any modulus.
    np.testing.assert_allclose(h / k, (2.0 * degree_l + 1.0) / 3.0, rtol=1.0e-8)


@pytest.mark.parametrize('core_model', (0, 1, 2, 3))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_matrix_old_vs_new_parity(degree_l, core_model):
    """The `_x` matrix solver matches the original solver's matrix method for k, h, and l."""
    out_new = run_matrix(radial_solver_new, degree_l, core_model)
    out_old = run_matrix(radial_solver_old, degree_l, core_model)
    assert out_new.success and out_old.success
    np.testing.assert_allclose(out_new.love, out_old.love, rtol=1.0e-9, atol=1.0e-10)


@pytest.mark.parametrize('degree_l', (2, 3))
def test_matrix_vs_shooting(degree_l):
    """The matrix (static) and shooting (dynamic) methods agree to the inertial-term level (~1e-3)."""
    out_matrix = run_matrix(radial_solver_new, degree_l, core_model=0)
    out_shooting = radial_solver_new(
        radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
        frequency, density,
        ('solid',), (False,), (True,), upper_radius_by_layer,
        degree_l=degree_l, solve_for=('tidal',), use_kamata=True,
        integration_method='DOP853', integration_rtol=1.0e-10, integration_atol=1.0e-12,
        scale_rtols_bylayer_type=False, use_prop_matrix=False,
        verbose=False, nondimensionalize=True, raise_on_fail=True, warnings=False)
    assert out_matrix.success and out_shooting.success
    np.testing.assert_allclose(out_shooting.love, out_matrix.love, rtol=1.0e-3)
