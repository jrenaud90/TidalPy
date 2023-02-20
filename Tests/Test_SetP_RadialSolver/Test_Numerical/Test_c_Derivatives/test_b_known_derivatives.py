""" Test the known functions in the `TidalPy.radial_solver.numerical.derivatives` module. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.toolbox.conversions import days2rads
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.radial_solver.numerical.derivatives import known_multilayer_odes, find_ode
from TidalPy.radial_solver.numerical.derivatives.odes import (
    dynamic_liquid_ode, dynamic_solid_ode, static_liquid_ode, static_solid_ode)

# Known initial guess functions stored by: is_solid, is_static
# Incompressible flag is built into the outer odes.
correct_functions = {
    (True, True)  : static_solid_ode,
    (True, False) : dynamic_solid_ode,
    (False, True) : static_liquid_ode,
    (False, False): dynamic_liquid_ode
    }

N = 10
radius        = np.linspace(100., 1.0e6, N)
shear         = np.linspace(10.e9, 50.e9, N)
viscosity     = np.logspace(24., 18., N)
bulk          = np.linspace(1000.0e9, 100.0e9, N)
density       = np.linspace(9000., 2500., N)
gravity       = np.linspace(0.1, 10., N)
frequency     = days2rads(6.0)
complex_shear = maxwell(frequency, shear**(-1), viscosity)**(-1)

@pytest.mark.parametrize('is_solid', (True, False))
@pytest.mark.parametrize('is_static', (True, False))
def test_known_derivatives(is_solid, is_static):
    """ Test the known derivative functions. """

    # Get correct function and the one from TidalPy
    ode         = correct_functions[(is_solid, is_static)]
    ode_tidalpy = known_multilayer_odes[(is_solid, is_static)]

    # Check they are the same
    assert ode_tidalpy is ode

@pytest.mark.parametrize('is_solid', (True, False))
@pytest.mark.parametrize('is_static', (True, False))
def test_find_ode(is_solid, is_static):
    """ Test the `radial_solver.numerical.derivative` `find_ode` function. """

    # Find the ode and additional input from the function
    ode_tidalpy, additional_input = \
        find_ode(is_solid, is_static, False, radius, complex_shear, bulk, density, gravity, frequency,
                 order_l=2, G_to_use=G)

    # Compared to expected value
    ode = correct_functions[(is_solid, is_static)]
    assert ode_tidalpy is ode

    # Check the additional inputs
    if is_solid:
        if is_static:
            assert len(additional_input) == 8
        else:
            assert len(additional_input) == 9
    else:
        if is_static:
            assert len(additional_input) == 6
        else:
            assert len(additional_input) == 8
