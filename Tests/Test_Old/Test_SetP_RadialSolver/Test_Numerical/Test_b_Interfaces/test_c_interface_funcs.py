""" Test the `TidalPy.radial_solver.numerical.interfaces` functionality. All interface types will be tested. """

import numpy as np
import pytest

import TidalPy


from TidalPy.constants import G
from TidalPy.utilities.performance import nbList
from TidalPy.radial_solver.numerical.interfaces import find_interface_func

static_liquid_density, interface_gravity, G_to_use = 4000., 7.5, G

lower_layer_ys_dict = {
    # Stored by is_lower_solid, is_lower_static
    (True, True):   np.asarray(((10. - 1.j) * np.ones(6, dtype=np.complex128),
                                (15. - 2.j) * np.ones(6, dtype=np.complex128),
                                (12. - 3.j) * np.ones(6, dtype=np.complex128)
                                )),
    (True, False):  np.asarray(((13. - 4.j) * np.ones(6, dtype=np.complex128),
                                (18. - 5.j) * np.ones(6, dtype=np.complex128),
                                (15. - 6.j) * np.ones(6, dtype=np.complex128)
                                )),
    # Dynamic liquid layer only has 2 solutions with 4 ys defined.
    (False, False): np.asarray(((112. - 7.j) * np.ones(4, dtype=np.complex128),
                                (17. - 8.j) * np.ones(4, dtype=np.complex128)
                                )),
    # Static liquid layer only has 1 solution with 2 ys defined.
    (False, True):  np.asarray(((183. - 9.j) * np.ones(2, dtype=np.complex128),
                                ))
    }

@pytest.mark.parametrize('is_lower_solid', (True, False))
@pytest.mark.parametrize('is_upper_solid', (True, False))
@pytest.mark.parametrize('is_lower_static', (True, False))
@pytest.mark.parametrize('is_upper_static', (True, False))
def test_interface_functions(is_lower_solid, is_upper_solid, is_lower_static, is_upper_static):
    """ Test all interface functions from the `TidalPy.radial_solver.numerical.interfaces` module. """

    # Get radial solution inputs
    lower_layer_ys_by_solution = lower_layer_ys_dict[(is_lower_solid, is_lower_static)]
    num_lower_solutions = lower_layer_ys_by_solution.shape[0]

    # Get expected shapes and sizes for outputs
    upper_layer_ys_check_by_solution = lower_layer_ys_dict[(is_upper_solid, is_upper_static)]
    num_upper_solutions = upper_layer_ys_check_by_solution.shape[0]

    # Find interface function
    interface_func, extra_inputs = \
        find_interface_func(
            is_lower_solid, is_lower_static, is_upper_solid, is_upper_static,
            static_liquid_density, interface_gravity, G_to_use=G)

    # Calculate upper layer solution based on the interface function
    upper_layer_ys_by_solution = interface_func(lower_layer_ys_by_solution, *extra_inputs)
    num_upper_solutions_calculated = upper_layer_ys_by_solution.shape[0]

    # Check number of solutions
    assert num_upper_solutions_calculated == num_upper_solutions
    assert upper_layer_ys_by_solution.shape == upper_layer_ys_check_by_solution.shape
    assert upper_layer_ys_by_solution.dtype == np.complex128

    # Perform more specialized checks
    if is_lower_static and is_upper_static:
        if is_lower_solid and is_upper_solid:
            # The lower and upper ys should be the same
            assert np.all(upper_layer_ys_by_solution == lower_layer_ys_by_solution)

        if not is_lower_solid and not is_upper_solid:
            # The lower and upper ys should be the same
            assert np.all(upper_layer_ys_by_solution == lower_layer_ys_by_solution)
