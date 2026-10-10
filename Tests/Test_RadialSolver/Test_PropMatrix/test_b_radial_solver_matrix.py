"""The propagation matrix method runs on a 1-layer static incompressible solid for each core model and solve type."""
import numpy as np
import pytest

from TidalPy.RadialSolver.solver import radial_solver
from TidalPy.Rheology import Maxwell

frequency = 1.0 / (86400. * 1.5)
N = 25
radius_array = np.linspace(0.0, 6000.e3, N)
density_array = 5500. * np.ones_like(radius_array)
bulk_modulus_array = 1.0e14 * np.ones(radius_array.size, dtype=np.complex128, order='C')
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)
complex_shear_modulus_array = Maxwell().calc_complex_modulus_vectorize_modulus(shear_array, viscosity_array, frequency)
planet_bulk_density = float(density_array[0])
upper_radius_by_layer = np.asarray((radius_array[-1],))


def _check_matrix_solve(solve_for, **kwargs):
    """Run the matrix method and check the output types and shape; skip unsupported inputs."""
    try:
        out = radial_solver(
            radius_array,
            density_array,
            bulk_modulus_array,
            complex_shear_modulus_array,
            frequency,
            planet_bulk_density,
            ('solid',),
            (True,),
            (True,),
            upper_radius_by_layer,
            solve_for=solve_for,
            love_method='propagation_matrix',
            verbose=False,
            raise_on_fail=True,
            **kwargs)

        assert out.success
        assert type(out.message) is str
        assert type(out.result) is np.ndarray
        assert out.result.shape == (len(solve_for) * 6, N)
    except NotImplementedError as error:
        pytest.skip(f'function does not currently support requested inputs. Skipping Test. Details: {error}')


@pytest.mark.parametrize('core_model', (0, 1, 2, 3))
@pytest.mark.parametrize('nondimensionalize', (True, False))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('solve_for', (('free',), ('tidal',), ('loading',)))
def test_radial_solver_matrix_1layer(core_model, nondimensionalize, degree_l, solve_for):
    """A single solve type succeeds with a (6, N) result."""
    _check_matrix_solve(solve_for, degree_l=degree_l, core_model=core_model, nondimensionalize=nondimensionalize)


@pytest.mark.parametrize('core_model', (0, 1, 2, 3))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_radial_solver_matrix_1layer_solve_for_both(core_model, degree_l):
    """Solving for tidal and loading together succeeds with a (12, N) result."""
    # Core models other than 0 require the automatic starting radius (a manual one fails with error -22).
    _check_matrix_solve(
        ('tidal', 'loading'),
        degree_l=degree_l,
        core_model=core_model,
        starting_radius=0.2 * radius_array[-1] if core_model == 0 else 0.0,
        nondimensionalize=False)
