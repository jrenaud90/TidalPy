""" Test the `TidalPy.radial_solver.numerical.initial` initial value calculations. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.toolbox.conversions import days2rads
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.radial_solver.numerical.initial import find_initial_guess

radius    = 100.
shear     = 50.e9
viscosity = 1.0e19
bulk      = 150.e9
density   = 5500.
frequency = days2rads(6.0)
complex_shear = maxwell(frequency, shear**(-1), viscosity)**(-1)


@pytest.mark.parametrize('is_solid', (True, False))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('is_kamata', (True, False))
@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_known_initial_funcs(is_solid, is_static, is_incompressible, is_kamata, order_l):
    """ Test the `radial_solver.initial` initial solution functions. """

    # Skip the following tests since these functions have not been implemented yet.
    skip = False
    if is_incompressible:
        if is_solid:
            if is_static:
                skip = True
            else:
                if not is_kamata:
                    skip = True
        else:
            if not is_static:
                if not is_kamata:
                    skip = True
    if skip:
        pytest.skip("Functionality has not yet been implemented.")

    # Perform initial value calculations
    initial_solutions = find_initial_guess(
        is_solid, is_static, is_incompressible, is_kamata,
        radius, complex_shear, bulk, density, frequency,
        order_l=order_l, G_to_use=G)

    # Check if there are the expected number of solutions.
    if is_solid:
        # Solid (dynamic or static) will have three solutions.
        assert len(initial_solutions) == 3
    elif is_static:
        # Static liquid has one solution
        assert len(initial_solutions) == 1
    else:
        # Dynamic liquid has two solution
        assert len(initial_solutions) == 2

    # For each solution...
    for solution in initial_solutions:
        # Check if there is the right number of y-variables.
        if is_solid:
            # Solid (dynamic or static) has 6 y-values
            assert solution.shape == (6,)
        elif is_static:
            # Static liquid has 2 y-values
            assert solution.shape == (2,)
        else:
            # Dynamic liquid has 4 y-values
            assert solution.shape == (4,)

        # Check type
        assert solution.dtype == np.complex128
