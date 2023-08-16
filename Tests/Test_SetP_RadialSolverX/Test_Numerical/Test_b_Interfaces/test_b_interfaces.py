""" Test the `TidalPy.radial_solver.interfaces` (cython extension) functionality. """

import pytest
import numpy as np

import TidalPy

TidalPy.test_mode()

from TidalPy.radial_solver.numerical.interfaces.interfaces_x import find_solution_num, interface_x


@pytest.mark.parametrize('lower_is_solid', (True, False))
@pytest.mark.parametrize('lower_is_static', (True, False))
@pytest.mark.parametrize('lower_is_compressible', (True, False))
@pytest.mark.parametrize('upper_is_solid', (True, False))
@pytest.mark.parametrize('upper_is_static', (True, False))
@pytest.mark.parametrize('upper_is_compressible', (True, False))
def test_interfaces(lower_is_solid, lower_is_static, lower_is_compressible,
                    upper_is_solid, upper_is_static, upper_is_compressible):
    """ Test interface function """

    upper_equals_lower = False
    if lower_is_solid:
        if upper_is_solid:
            upper_equals_lower = True
        if lower_is_static:
            if lower_is_compressible:
                y_lower = np.asarray((
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        )
                    ), dtype=np.complex128)
            else:
                # TODO: Confirm
                y_lower = np.asarray((
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        )
                    ), dtype=np.complex128)
        else:
            # Dynamic
            if lower_is_compressible:
                y_lower = np.asarray((
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        )
                    ), dtype=np.complex128)
            else:
                # TODO: Confirm
                y_lower = np.asarray((
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        )
                    ), dtype=np.complex128)
    else:
        # Liquid
        if not upper_is_solid:
            upper_equals_lower = True
        if lower_is_static:
            if lower_is_compressible:
                y_lower = np.asarray((
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        )
                    ), dtype=np.complex128)
            else:
                # TODO: Confirm
                y_lower = np.asarray((
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        )
                    ), dtype=np.complex128)
        else:
            # Dynamic
            if lower_is_compressible:
                y_lower = np.asarray((
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        )
                    ), dtype=np.complex128)
            else:
                # TODO: Confirm
                y_lower = np.asarray((
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        ),
                    (
                        1.0 + 2.0j,
                        3.0 + 2.0j,
                        1.0 + 2.0j,
                        3.0 + 2.0j
                        )
                    ), dtype=np.complex128)

    interface_gravity = 4.0
    liquid_density = 3000.

    upper_sols = find_solution_num(upper_is_solid, upper_is_static, upper_is_compressible)
    y_upper = np.empty((upper_sols, upper_sols*2), dtype=np.complex128)

    interface_x(y_lower, y_upper,
                lower_is_solid, lower_is_static, lower_is_compressible,
                upper_is_solid, upper_is_static, upper_is_compressible,
                interface_gravity, liquid_density)

    if upper_equals_lower:
        assert np.allclose(y_lower, y_upper)
    elif lower_is_solid and not upper_is_solid:
        # Lower solid / upper liquid
        if upper_is_static:
            # TODO LEFT OFF


    assert y_upper
    assert sol_num_check == num_sols
