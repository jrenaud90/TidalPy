""" Test the `TidalPy.radial_solver.interfaces` `solution_num` func (cython extension) functionality. """

import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.radial_solver.numerical.interfaces.interfaces_x import solution_num


@pytest.mark.parametrize('is_solid', (True, False))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_compressible', (True, False))
def test_solution_num(is_solid, is_static, is_compressible):
    """ Test solution number function """

    sol_num_check = solution_num(is_solid, is_static, is_compressible)

    if is_solid:
        if is_static:
            if is_compressible:
                num_sols = 3
            else:
                # TODO: Confirm
                num_sols = 3
        else:
            # Dynamic
            if is_compressible:
                num_sols = 3
            else:
                # TODO: Confirm
                num_sols = 3
    else:
        # Liquid
        if is_static:
            if is_compressible:
                num_sols = 1
            else:
                # TODO: Confirm
                num_sols = 1
        else:
            # Dynamic
            if is_compressible:
                num_sols = 2
            else:
                # TODO: Confirm
                num_sols = 2
    assert sol_num_check == num_sols
