""" Test the `TidalPy.radial_solver.numerical.derivatives` functionality. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.toolbox.conversions import days2rads
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.radial_solver.numerical.initial import find_initial_guess
from TidalPy.radial_solver.numerical.derivatives import find_ode

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
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_derivatives(is_solid, is_static, is_incompressible, order_l):
    """ Test all `radial_solver.numerical.derivatives` ODE functions. """

    # Get correct ode function
    ode_func, additional_input = \
        find_ode(
            is_solid, is_static, False, radius, complex_shear, bulk, density, gravity, frequency,
            order_l=2, G_to_use=G)

    # Get initial y value
    if is_solid:
        y_init_solutions = find_initial_guess(True, True, False, False,
                                            radius[0], complex_shear[0], bulk[0], density[0], frequency,
                                            order_l=order_l, G_to_use=G)
    elif is_static:
        y_init_solutions = find_initial_guess(
            False, True, False, False,
            radius[0], complex_shear[0], bulk[0], density[0], frequency,
            order_l=order_l, G_to_use=G)
    else:
        y_init_solutions = find_initial_guess(
                False, False, False, False,
                radius[0], complex_shear[0], bulk[0], density[0], frequency,
                order_l=order_l, G_to_use=G)

    # Convert M complex-valued initial conditions to 2M floats
    y_init = y_init_solutions[0]
    y_init_float = np.empty(2 * len(y_init), dtype=np.float64)

    for i in range(len(y_init)):
        y_init_float[2 * i]     = np.real(y_init[i])
        y_init_float[2 * i + 1] = np.imag(y_init[i])

    # Perform calculation at multiple r values
    import pdb;
    for r_i in (2, 5, 8):
        assert r_i < N
        ode_result = ode_func(radius[r_i], y_init_float, *additional_input)

        # Check shape
        if is_solid:
            # Solid radial solution has 6 complex ys x 2 for conversion to floats
            assert ode_result.shape == (2 * 6,)
        elif is_static:
            # Static liquid radial solution has 2 complex ys x 2 for conversion to floats
            assert ode_result.shape == (2 * 2,)
        else:
            # Dynamic liquid radial solution has 4 complex ys x 2 for conversion to floats
            assert ode_result.shape == (2 * 4,)

        # Check type
        assert ode_result.dtype == np.float64

        # Check specific values
        if order_l == 2 and r_i == 2:
            if is_solid:
                if is_static:
                    if is_incompressible:
                        expected_value = np.asarray([
                            3064.809673825056, -1.834800256574697e-06, 108011850951755.95, 47096.47434044733,
                            6559.012727857043, -0.00021463894247496384, -2786797157.2015085, 15.165539584674478,
                            -14295381529.773653, -6.237529599143017, -64324.07099404163, -2.8066637864772505e-05],
                                dtype=np.float64)
                    else:
                        expected_value = np.asarray([
                            3064.809673825056, -1.834800256574697e-06, 108011850951755.95, 47096.47434044733,
                            6559.012727857043, -0.00021463894247496384, -2786797157.2015085, 15.165539584674478,
                            -14295381529.773653, -6.237529599143017, -64324.07099404163, -2.8066637864772505e-05],
                                dtype=np.float64)
                else:
                    if is_incompressible:
                        expected_value = np.asarray([
                            3064.809673825056, -1.834800256574697e-06, 108011850951755.95, 47096.47434044733,
                            6559.012727857043, -0.00021463894247496384, -2786797157.2015085, 15.165539584674478,
                            -14295381529.773653, -6.237529599143017, -64324.07099404163, -2.8066637864772505e-05],
                                dtype=np.float64)
                    else:
                        expected_value = np.asarray([
                            3064.809673825056, -1.834800256574697e-06, 108011850951755.95, 47096.47434044733,
                            6559.012727857043, -0.00021463894247496384, -2786797157.2015085, 15.165539584674478,
                            -14295381529.773653, -6.237529599143017, -64324.07099404163, -2.8066637864772505e-05],
                                dtype=np.float64)
            else:
                if is_static:
                    if is_incompressible:
                        expected_value = np.asarray([
                            199.89259929850508, 0.000e+00, 0.00034889169130223523, 0.000e+00], dtype=np.float64)
                    else:
                        expected_value = np.asarray([
                            199.89259929850508, 0.000e+00, 0.00034889169130223523, 0.000e+00], dtype=np.float64)
                else:
                    if is_incompressible:
                        expected_value = np.asarray([
                            -9435468.451139623, 0.000e+00, -139467332446.44806, 0.000e+00, -3242673.621620953,
                            0.000e+00, 45.20146789341666, 0.000e+00], dtype=np.float64)
                    else:
                        expected_value = np.asarray([
                            -9435468.451139623, 0.000e+00, -139467332446.44806, 0.000e+00, -3242673.621620953,
                            0.000e+00, 45.20146789341666, 0.000e+00], dtype=np.float64)

            assert np.allclose(ode_result, expected_value)
