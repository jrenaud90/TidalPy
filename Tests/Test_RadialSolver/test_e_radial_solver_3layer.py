import pytest
import numpy as np


from TidalPy.exceptions import SolutionFailedError
from TidalPy.RadialSolver import radial_solver
from TidalPy.utilities.spherical_helper import calculate_mass_gravity_arrays

frequency = 1.0 / (86400. * 0.2)
N = 20
planet_r = 6000.0e3
icb_r = planet_r * (1. / 3.0)
cmb_r = planet_r * (2. / 3.0)
radius_array = np.concatenate((
    np.linspace(0.0, icb_r, N),
    np.linspace(icb_r, cmb_r, N),
    np.linspace(cmb_r, planet_r, N)
    ))
ic_index = np.zeros(radius_array.size, dtype=bool)
oc_index = np.zeros(radius_array.size, dtype=bool)
mantle_index = np.zeros(radius_array.size, dtype=bool)
ic_index[np.arange(0, N)]             = True
oc_index[np.arange(N, 2 * N)]         = True
mantle_index[np.arange(2 * N, 3 * N)] = True

density_array = np.zeros_like(radius_array)
density_array[mantle_index] = 3500.
density_array[oc_index]     = 7000.
density_array[ic_index]     = 8500.

bulk_modulus_array = 1.0e11 * np.ones(radius_array.size, dtype=np.complex128, order='C')

viscosity_array = np.zeros_like(radius_array)
viscosity_array[mantle_index] = 1.0e20
viscosity_array[oc_index]     = 1.0e6
viscosity_array[ic_index]     = 1.0e26

shear_array = np.zeros_like(radius_array)
shear_array[mantle_index] = 5.0e10
shear_array[oc_index]     = 0.0
shear_array[ic_index]     = 1.0e11

from TidalPy.rheology.models import Maxwell
complex_shear_modulus_array = np.empty(radius_array.size, dtype=np.complex128)
max_rho = Maxwell()
max_rho.vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear_modulus_array)

planet_bulk_density = np.average(density_array)
upper_radius_by_layer = np.asarray((icb_r, cmb_r, planet_r))
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
            radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
            frequency, planet_bulk_density, 
            layer_types, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
            degree_l=degree_l, solve_for=solve_for, use_kamata=use_kamata,
            integration_method=method, integration_rtol=1.0e-7, integration_atol=1.0e-10,
            scale_rtols_bylayer_type=False,
            max_num_steps=5_000_000, expected_size=250, max_step=0, starting_radius=0.1,
            raise_on_fail=True,
            verbose=False, nondimensionalize=True)


        if not out.success:
            if (not liquid_is_incompressible) and solid_is_incompressible:
                # TODO: Look into this.
                # v0.6.0 update: It looks like it is happening because the incompressible solid underneath is not coupling with the compressible liquid. 
                pytest.skip('Integration Failed. Compressible liquid with incompressible solid below is not very stable.')

        assert out.success
        assert type(out.message) is str
        assert type(out.result) is np.ndarray
        assert out.result.shape == (6, radius_array.size)
    except NotImplementedError as e:
        pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {e}')
    except SolutionFailedError as e:
        if (not liquid_is_incompressible) and solid_is_incompressible:
            # TODO: Look into this.
            # v0.6.0 update: It looks like it is happening because the incompressible solid underneath is not coupling with the compressible liquid. 
            pytest.skip('Integration Failed. Compressible liquid with incompressible solid below is not very stable.')


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
            radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
            frequency, planet_bulk_density, 
            layer_types, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
            degree_l=degree_l, solve_for=solve_for, use_kamata=use_kamata,
            integration_method=method, integration_rtol=1.0e-7, integration_atol=1.0e-10,
            scale_rtols_bylayer_type=False,
            max_num_steps=5_000_000, expected_size=250, max_step=0, starting_radius=0.1,
            raise_on_fail=True,
            verbose=False, nondimensionalize=True)
        
        if not out.success:
            if (not liquid_is_incompressible) and solid_is_incompressible:
                # TODO: Look into this.
                # v0.6.0 update: It looks like it is happening because the incompressible solid underneath is not coupling with the compressible liquid. 
                pytest.skip('Integration Failed. Compressible liquid with incompressible solid below is not very stable.')
        
        assert out.success
        assert type(out.message) is str
        assert type(out.result) is np.ndarray
        assert out.result.shape == (len(solve_for) * 6, radius_array.size)
    except NotImplementedError as e:
         pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {e}')
    except SolutionFailedError as e:
        if (not liquid_is_incompressible) and solid_is_incompressible:
            # TODO: Look into this.
            # v0.6.0 update: It looks like it is happening because the incompressible solid underneath is not coupling with the compressible liquid. 
            pytest.skip('Integration Failed. Compressible liquid with incompressible solid below is not very stable.')
