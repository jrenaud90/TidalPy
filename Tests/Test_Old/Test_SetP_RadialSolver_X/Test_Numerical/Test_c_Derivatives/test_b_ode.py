""" Test the `TidalPy.radial_solver.ode` (cython extension) functionality. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.radial_solver.numerical.derivatives import get_radial_solver_class

N = 10
radius_array = np.linspace(0.1, 10., N)
shear_modulus_array = 1.0e10 * np.ones(radius_array.shape, dtype=np.complex128)
density_array = 3000. * np.ones(radius_array.shape, dtype=np.float64)
gravity_array = np.linspace(0.1, 2., N, dtype=np.float64)
bulk_modulus_array = 1.0e11 * np.ones(radius_array.shape, dtype=np.float64)

frequency = 1.
degree_l = 2
G_to_use = 1.0e-12

t_span = (0.1, 10.)

y0 = np.asarray((
    1.0, -0.2,
    2.0, -0.3,
    1.0, -0.2,
    2.0, -0.3,
    1.0, -0.2,
    2.0, -0.3), dtype=np.float64)

@pytest.mark.parametrize('rk_method', (1, 2))
@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('limit_solution_to_radius', (False, True))
@pytest.mark.parametrize('is_incompressible', (False, True))
@pytest.mark.parametrize('is_static', (False, True))
@pytest.mark.parametrize('is_solid', (False, True))
def test_solid_static_compressible_solver(
        is_solid,
        is_static,
        is_incompressible,
        limit_solution_to_radius,
        degree_l,
        rk_method):

    num_y, SolverClass = get_radial_solver_class(is_solid, is_static, is_incompressible)
    y0_reduced = y0[:num_y]

    solver = SolverClass(
        radius_array, density_array, gravity_array, shear_modulus_array, bulk_modulus_array, frequency, degree_l,
        G_to_use, t_span, y0_reduced, rk_method=rk_method, auto_solve=True, limit_solution_to_radius=limit_solution_to_radius
        )
    solver.solve()
    assert solver.success

    if limit_solution_to_radius:
        assert solver.t.shape == (N,)
        assert solver.y.shape == (num_y, N)
    else:
        assert solver.y.shape[0] == num_y
        assert solver.t.size > N

    assert type(solver.t[0]) is np.float64
    assert type(solver.y[0, 0]) is np.float64
    assert solver.t[0] == radius_array[0]
    assert solver.t[-1] == radius_array[-1]
    assert np.all(np.isfinite(solver.t))
    assert np.all(np.isfinite(solver.y))
