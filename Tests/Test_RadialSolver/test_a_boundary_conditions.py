import pytest
from math import isnan

import numpy as np

from TidalPy.RadialSolver.boundaries.surface_bc import get_surface_bc


@pytest.mark.parametrize('degree_l', (2, 3))
def test_get_surface_bc(degree_l):
    """ Test boundary condition finder function. """

    radius = 1000.
    density = 2000.

    def check_val(array, model_type):
        if model_type == 0:
            for i in range(3):
                assert array[i] == 0.
        elif model_type == 1:
            for i in range(2):
                assert array[i] == 0.
            assert array[2] == (2. * degree_l + 1.) / radius
        elif model_type == 2:
            assert array[0] == (-1. / 3.) * (2. * degree_l + 1.) * density
            assert array[1] == 0.
            assert array[2] == (2. * degree_l + 1.) / radius
        else:
            raise NotImplementedError

        return True

    # Test with 1 model
    for model_type in (0, 1, 2):
        bc_models = np.asarray((model_type,), dtype=np.intc)
        boundary_condition_array = get_surface_bc(
            bc_models,
            radius,
            density,
            degree_l,
            )
        assert check_val(boundary_condition_array[:3], model_type=model_type)
        for value in boundary_condition_array[3:]:
            assert isnan(value)
    
    # Test with 2 models
    for model_pair in ((0, 0), (0, 1), (1, 0), (1, 1,), (1, 2), (2, 1), (2, 2)):
        bc_models = np.asarray(model_pair, dtype=np.intc)
        boundary_condition_array = get_surface_bc(
            bc_models,
            radius,
            density,
            degree_l,
            )
        for i, model_type in enumerate(model_pair):
            assert check_val(boundary_condition_array[3*i:(3*i)+3], model_type=model_type)
    
    # Test with 3 models
    for model_pair in ((0, 0, 0), (0, 1, 2), (1, 1, 1), (2, 2, 2), (2, 1, 0)):
        bc_models = np.asarray(model_pair, dtype=np.intc)
        boundary_condition_array = get_surface_bc(
            bc_models,
            radius,
            density,
            degree_l,
            )
        for i, model_type in enumerate(model_pair):
            assert check_val(boundary_condition_array[3*i:(3*i)+3], model_type=model_type)
