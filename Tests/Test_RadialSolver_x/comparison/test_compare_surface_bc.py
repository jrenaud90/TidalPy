"""Compare RadialSolver_x boundary conditions against original RadialSolver."""
import pytest
import numpy as np

from TidalPy.RadialSolver.boundaries.surface_bc import get_surface_bc as get_surface_bc_old
from TidalPy.RadialSolver_x.boundaries.surface_bc import get_surface_bc as get_surface_bc_new


@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('bc_models', ((0,), (1,), (2,), (0, 1), (1, 2), (0, 1, 2)))
def test_compare_surface_bc(degree_l, bc_models):
    """Surface boundary conditions should match exactly between old and new implementations."""
    radius = 1000.
    density = 2000.
    bc_arr = np.asarray(bc_models, dtype=np.intc)

    old_result = get_surface_bc_old(bc_arr, radius, density, degree_l)
    new_result = get_surface_bc_new(bc_arr, radius, density, degree_l)

    # These should be exact since it's just simple arithmetic.
    np.testing.assert_allclose(new_result, old_result, rtol=1e-15)
