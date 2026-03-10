"""Compare RadialSolver_x radial_solver output against original RadialSolver.

This is the most important comparison test: it runs the full radial solver from both
modules with identical inputs and compares Love numbers and radial solutions.
"""
import pytest
import numpy as np

from TidalPy.rheology.models import Maxwell
from TidalPy.utilities.spherical_helper import calculate_mass_gravity_arrays

from TidalPy.RadialSolver import radial_solver as radial_solver_old
from TidalPy.RadialSolver_x.solver import radial_solver as radial_solver_new

# ---------- Setup: 1-layer homogeneous solid planet ----------
frequency = 1.0 / (86400. * 0.2)  # 0.2-day period
N = 100
radius_array = np.linspace(0.0, 6000.e3, N)
density_array = 3500. * np.ones_like(radius_array)
bulk_modulus_array = 1.0e11 * np.ones(N, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)

max_rho = Maxwell()
complex_shear_modulus_array = np.empty(N, dtype=np.complex128)
max_rho.vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear_modulus_array)

volume_array, mass_array, gravity_array = calculate_mass_gravity_arrays(radius_array, density_array)
planet_bulk_density = np.sum(mass_array) / np.sum(volume_array)
upper_radius_by_layer = np.asarray((radius_array[-1],))


@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('solve_for', (('tidal',), ('loading',), ('tidal', 'loading')))
def test_compare_radial_solver_1layer_solid(is_static, is_incompressible, degree_l, use_kamata, solve_for):
    """Compare Love numbers and radial solution arrays for 1-layer solid planets."""

    layer_type_by_layer = ('solid',)
    is_static_by_layer = (is_static,)
    is_incompressible_by_layer = (is_incompressible,)

    common_kwargs = dict(
        degree_l=degree_l, solve_for=solve_for, use_kamata=use_kamata,
        integration_method='rk45', integration_rtol=1.0e-7, integration_atol=1.0e-10,
        scale_rtols_bylayer_type=False,
        max_num_steps=5_000_000, expected_size=250, max_step=0,
        verbose=False, nondimensionalize=True, starting_radius=0.1,
        raise_on_fail=True,
    )

    try:
        old_out = radial_solver_old(
            radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
            frequency, planet_bulk_density,
            layer_type_by_layer, is_static_by_layer, is_incompressible_by_layer,
            upper_radius_by_layer,
            **common_kwargs
        )
    except NotImplementedError:
        pytest.skip('Not implemented in original RadialSolver.')

    try:
        new_out = radial_solver_new(
            radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
            frequency, planet_bulk_density,
            layer_type_by_layer, is_static_by_layer, is_incompressible_by_layer,
            upper_radius_by_layer,
            **common_kwargs
        )
    except NotImplementedError:
        pytest.skip('Not implemented in RadialSolver_x.')

    # Both should succeed.
    assert old_out.success, f"Old solver failed: {old_out.message}"
    assert new_out.success, f"New solver failed: {new_out.message}"

    # Compare result arrays (radial solutions).
    assert old_out.result.shape == new_out.result.shape
    try:
        np.testing.assert_allclose(new_out.result, old_out.result, rtol=1e-5, atol=1e-20,
                               err_msg="Radial solution arrays differ.")
    except:
        import pdb; pdb.set_trace()

    # Compare Love numbers (both return ndarray of shape (num_solve_for, 3)).
    assert old_out.love.shape == new_out.love.shape
    np.testing.assert_allclose(new_out.love, old_out.love, rtol=1e-5,
                               err_msg="Love numbers differ.")
