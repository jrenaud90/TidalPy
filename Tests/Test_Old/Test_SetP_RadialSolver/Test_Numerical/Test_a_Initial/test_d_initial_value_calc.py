""" Test the `TidalPy.radial_solver.numerical.initial` initial value calculations. """

import numpy as np
import pytest

import TidalPy


from TidalPy.constants import G
from TidalPy.utilities.conversions import days2rads
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

    assert type(initial_solutions) == np.ndarray

    # Check if there are the expected number of solutions.
    if is_solid:
        # Solid (dynamic or static) will have three solutions.
        assert initial_solutions.shape[0] == 3
        # Solid (dynamic or static) has 6 y-values
        assert initial_solutions.shape[1] == 6
    elif is_static:
        # Static liquid has one solution
        assert initial_solutions.shape[0] == 1
        # Static liquid has 2 y-values
        assert initial_solutions.shape[1] == 2
    else:
        # Dynamic liquid has two solution
        assert initial_solutions.shape[0] == 2
        # Dynamic liquid has 4 y-values
        assert initial_solutions.shape[1] == 4

    # Check type
    assert initial_solutions.dtype == np.complex128
