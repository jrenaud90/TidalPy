""" Test the `TidalPy.radial_solver.interfaces` `solution_num` func (cython extension) functionality. """

import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.radial_solver.numerical.interfaces.interfaces_x import find_solution_num


@pytest.mark.parametrize('is_solid', (True, False))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (False, True))
def test_solution_num(is_solid, is_static, is_incompressible):
    """ Test solution number function """

    sol_num_check = find_solution_num(is_solid, is_static, is_incompressible)

    if is_solid:
        if is_static:
            if is_incompressible:
                num_sols = 3
            else:
                num_sols = 3
        else:
            # Dynamic
            if is_incompressible:
                num_sols = 3
            else:
                num_sols = 3
    else:
        # Liquid
        if is_static:
            if is_incompressible:
                num_sols = 1
            else:
                num_sols = 1
        else:
            # Dynamic
            if is_incompressible:
                num_sols = 2
            else:
                num_sols = 2
    assert sol_num_check == num_sols
