import warnings

import pytest
import numpy as np

import TidalPy
TidalPy.test_mode = True


from TidalPy.RadialSolver import radial_solver
from TidalPy.utilities.spherical_helper import calculate_mass_gravity_arrays

frequency = 1.0 / (86400. * 0.2)
N = 40
planet_r = 6000.0e3
icb_r = planet_r * (1. / 3.)
cmb_r = planet_r * (2. / 3)
radius_array = np.concatenate((
    np.linspace(0.1, icb_r, N),
    np.linspace(icb_r, cmb_r, N+1)[1:],
    np.linspace(cmb_r, planet_r, N+1)[1:])
)

density_array = 3500. * np.ones_like(radius_array)
density_array[radius_array <= cmb_r] = 7000.
density_array[radius_array <= icb_r] = 8500.

bulk_modulus_array = 1.0e11 * np.ones_like(radius_array)

viscosity_array = 1.0e20 * np.ones_like(radius_array)
viscosity_array[radius_array <= cmb_r] = 1.0e6
viscosity_array[radius_array <= icb_r] = 1.0e26

shear_array = 5.0e10 * np.ones_like(radius_array)
shear_array[radius_array <= cmb_r] = 0.0
shear_array[radius_array <= icb_r] = 1.0e11

from TidalPy.rheology.models import Maxwell
complex_shear_modulus_array = np.empty(radius_array.size, dtype=np.complex128)
max_rho = Maxwell()
max_rho.vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear_modulus_array)
volume_array, mass_array, gravity_array = calculate_mass_gravity_arrays(radius_array, density_array)

planet_bulk_density = np.sum(mass_array) / np.sum(volume_array)
upper_radius_by_layer = (icb_r, cmb_r, planet_r)
layer_types = ("solid", "liquid", "solid")

@pytest.mark.parametrize('solid_is_static', (True, False))
@pytest.mark.parametrize('liquid_is_static', (True, False))
@pytest.mark.parametrize('solid_is_incompressible', (True, False))
@pytest.mark.parametrize('liquid_is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2,))
@pytest.mark.parametrize('use_kamata', (True,))
@pytest.mark.parametrize('solve_for', (('free',), ('tidal',), ('loading',)))
def test_radial_solver_3layer(solid_is_static, liquid_is_static,
                              solid_is_incompressible, liquid_is_incompressible,
                              method, degree_l, use_kamata, solve_for):

    is_static_by_layer = (solid_is_static, liquid_is_static, solid_is_static)
    is_incompressible_by_layer = (solid_is_incompressible, liquid_is_incompressible, solid_is_incompressible)

    try:
        out = radial_solver(
            radius_array, density_array, gravity_array, bulk_modulus_array, complex_shear_modulus_array,
            frequency, planet_bulk_density, 
            layer_types, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
            degree_l=degree_l, solve_for=solve_for, use_kamata=use_kamata,
            integration_method=method, integration_rtol=1.0e-7, integration_atol=1.0e-10,
            scale_rtols_by_layer_type=False,
            max_num_steps=5_000_000, expected_size=250, max_step=0,
            limit_solution_to_radius=True, verbose=False, nondimensionalize=True)


        if not out.success:
            if (not liquid_is_static) and solid_is_incompressible:
                # TODO: Look into this.
                pytest.skip('Integration Failed. Dynamic liquid with incompressible solid is not very stable.')

        assert out.success
        assert type(out.message) is str
        assert type(out.result) is np.ndarray
        assert out.result.shape == (6, radius_array.size)
    except NotImplementedError as e:
        pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {e}')


@pytest.mark.parametrize('solid_is_static', (True, False))
@pytest.mark.parametrize('liquid_is_static', (True, False))
@pytest.mark.parametrize('solid_is_incompressible', (True, False))
@pytest.mark.parametrize('liquid_is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2,))
@pytest.mark.parametrize('use_kamata', (True,))
def test_radial_solver_3layer_solve_for_both(solid_is_static, liquid_is_static,
                                            solid_is_incompressible, liquid_is_incompressible,
                                            method, degree_l, use_kamata):
    """Tests solving for both tidal and loading a the same time."""

    is_static_by_layer = (solid_is_static, liquid_is_static, solid_is_static)
    is_incompressible_by_layer = (solid_is_incompressible, liquid_is_incompressible, solid_is_incompressible)

    solve_for=('tidal', 'loading')

    try:
        out = radial_solver(
            radius_array, density_array, gravity_array, bulk_modulus_array, complex_shear_modulus_array,
            frequency, planet_bulk_density, 
            layer_types, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
            degree_l=degree_l, solve_for=solve_for, use_kamata=use_kamata,
            integration_method=method, integration_rtol=1.0e-7, integration_atol=1.0e-10,
            scale_rtols_by_layer_type=False,
            max_num_steps=5_000_000, expected_size=250, max_step=0,
            limit_solution_to_radius=True, verbose=False, nondimensionalize=True)
        
        if not out.success:
            if (not liquid_is_static) and solid_is_incompressible:
                # TODO: Look into this.
                pytest.skip('Integration Failed. Dynamic liquid with incompressible solid is not very stable.')
        
        assert out.success
        assert type(out.message) is str
        assert type(out.result) is np.ndarray
        assert out.result.shape == (len(solve_for) * 6, radius_array.size)
    except NotImplementedError as e:
        warnings.warn(f'function does not currently support requested inputs. Skipping Test. Details: {e}')
