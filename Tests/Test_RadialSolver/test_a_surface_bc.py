"""Surface boundary conditions for one, two, and three simultaneous boundary models."""
import numpy as np
import pytest

from TidalPy.RadialSolver.boundaries.surface_bc import get_surface_bc

RADIUS = 1000.
DENSITY = 2000.


def _expected_bc(model_type, degree_l):
    """The three boundary values for one model: 0 free, 1 tidal, 2 loading."""
    if model_type == 0:
        return (0., 0., 0.)
    if model_type == 1:
        return (0., 0., (2. * degree_l + 1.) / RADIUS)
    return ((-1. / 3.) * (2. * degree_l + 1.) * DENSITY, 0., (2. * degree_l + 1.) / RADIUS)


@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('bc_models', (
    (0,), (1,), (2,),
    (0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2),
    (0, 0, 0), (0, 1, 2), (1, 1, 1), (2, 2, 2), (2, 1, 0),
))
def test_get_surface_bc(degree_l, bc_models):
    """Each model's three boundary values are exact."""
    boundary_condition_array = get_surface_bc(np.asarray(bc_models, dtype=np.intc), RADIUS, DENSITY, degree_l)
    for index, model_type in enumerate(bc_models):
        values = boundary_condition_array[3 * index:3 * index + 3]
        for value, expected in zip(values, _expected_bc(model_type, degree_l)):
            assert value == expected
