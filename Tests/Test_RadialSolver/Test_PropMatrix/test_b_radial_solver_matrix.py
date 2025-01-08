import pytest
import numpy as np

from TidalPy.RadialSolver import radial_solver

frequency = 1.0 / (86400. * 1.5)
N = 25
radius_array = np.linspace(0.0, 6000.e3, N)
density_array = 5500. * np.ones_like(radius_array)
bulk_modulus_array = 1.0e14 * np.ones(radius_array.size, dtype=np.complex128, order='C')
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

layer_type_by_layer = ('solid',)
is_static_by_layer = (True,)
is_incompressible_by_layer = (True,)
is_incompressible2_by_layer = (False,)

@pytest.mark.parametrize('core_model', (0,1,2,3))
@pytest.mark.parametrize('nondimensionalize', (True, False,))
@pytest.mark.parametrize('degree_l', (2,3))
@pytest.mark.parametrize('solve_for', (('free',), ('tidal',), ('loading',)))
def test_radial_solver_matrix_1layer(core_model, nondimensionalize, degree_l, solve_for):

    try:
        out = radial_solver(
            radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
            frequency, planet_bulk_density, 
            layer_type_by_layer, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
            degree_l=degree_l, solve_for=solve_for, core_model=core_model,
            use_prop_matrix=True,
            verbose=False, nondimensionalize=nondimensionalize, raise_on_fail=True)

        assert out.success
        assert type(out.message) is str
        assert type(out.result) is np.ndarray
        assert out.result.shape == (6, N)
    except NotImplementedError as e:
        pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {e}')


@pytest.mark.parametrize('core_model', (0,1,2,3))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_radial_solver_matrix_1layer_solve_for_both(core_model, degree_l):
    """Tests solving for both tidal and loading a the same time."""

    solve_for=('tidal', 'loading')

    try:
        out = radial_solver(
            radius_array, density_array, bulk_modulus_array, complex_shear_modulus_array,
            frequency, planet_bulk_density, 
            layer_type_by_layer, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
            degree_l=degree_l, solve_for=solve_for, core_model=core_model, starting_radius=0.1,
            use_prop_matrix=True,
            verbose=False, nondimensionalize=False, raise_on_fail=True)

        assert out.success
        assert type(out.message) is str
        assert type(out.result) is np.ndarray
        assert out.result.shape == (len(solve_for) * 6, N)
    except NotImplementedError as e:
        pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {e}')
