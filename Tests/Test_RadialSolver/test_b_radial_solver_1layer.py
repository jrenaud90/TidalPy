import pytest
import numpy as np


from TidalPy.RadialSolver import radial_solver

frequency = 1.0 / (86400. * 0.2)
N = 100
radius_array = np.linspace(0.0, 6000.e3, N)
density_array = 3500. * np.ones_like(radius_array)
bulk_modulus_array = 1.0e11 * np.ones(radius_array.size, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)

from TidalPy.rheology.models import Maxwell
complex_shear_modulus_array = np.empty(N, dtype=np.complex128)
max_rho = Maxwell()
max_rho.vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear_modulus_array)


from TidalPy.utilities.spherical_helper import calculate_mass_gravity_arrays
volume_array, mass_array, gravity_array = calculate_mass_gravity_arrays(radius_array, density_array)

planet_bulk_density = np.sum(mass_array) / np.sum(volume_array)
upper_radius_by_layer = np.asarray((radius_array[-1],))

@pytest.mark.parametrize('layer_type', ("solid", "liquid"))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('solve_for', (('free',), ('tidal',), ('loading',)))
def test_radial_solver_1layer(layer_type, is_static, is_incompressible, method, degree_l,
                              use_kamata, solve_for):

    layer_type_by_layer = (layer_type,)
    is_static_by_layer = (is_static,)
    is_incompressible_by_layer = (is_incompressible,)

    # TODO: Currently very unstable for 1-layer planets that are all liquid. For now, skip.
    if layer_type != 'solid':
        pytest.skip(f'Planets with 1-layer liquid are not currently very stable. Skipping tests.')
    else:
        try:
            out = radial_solver(
                radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
                frequency, planet_bulk_density, 
                layer_type_by_layer, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
                degree_l=degree_l, solve_for=solve_for, use_kamata=use_kamata,
                integration_method=method, integration_rtol=1.0e-7, integration_atol=1.0e-10,
                scale_rtols_bylayer_type=False,
                max_num_steps=5_000_000, expected_size=250, max_step=0,
                verbose=False, nondimensionalize=True, starting_radius=0.1,
                raise_on_fail=True,
                log_info=True  # For this test lets also check that logging info kwarg works.
            )

            assert out.success
            assert type(out.message) is str
            assert type(out.result) is np.ndarray
            assert out.result.shape == (6, N)
        except NotImplementedError as e:
            pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {e}')

@pytest.mark.parametrize('layer_type', ("solid", "liquid"))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('method', ("rk23", "rk45", "dop853"))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_kamata', (True, False))
def test_radial_solver_1layer_solve_for_both(layer_type, is_static, is_incompressible,
                                             method, degree_l, use_kamata):
    """Tests solving for both tidal and loading a the same time."""

    layer_type_by_layer = (layer_type,)
    is_static_by_layer = (is_static,)
    is_incompressible_by_layer = (is_incompressible,)

    solve_for=('tidal', 'loading')

    # TODO: Currently very unstable for 1-layer planets that are all liquid. For now, skip.
    if layer_type != 'solid':
        pytest.skip(f'Planets with 1-layer liquid are not currently very stable. Skipping tests.')
    else:
        try:
            out = radial_solver(
                radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
                frequency, planet_bulk_density, 
                layer_type_by_layer, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
                degree_l=degree_l, solve_for=solve_for, use_kamata=use_kamata,
                integration_method=method, integration_rtol=1.0e-7, integration_atol=1.0e-10,
                scale_rtols_bylayer_type=False,
                max_num_steps=5_000_000, expected_size=250, max_step=0, starting_radius=0.1,
                raise_on_fail=True,
                verbose=False, nondimensionalize=True)

            assert out.success
            assert type(out.message) is str
            assert type(out.result) is np.ndarray
            assert out.result.shape == (len(solve_for) * 6, N)
        except NotImplementedError as e:
            pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {e}')
