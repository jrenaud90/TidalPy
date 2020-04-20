import numpy as np
import TidalPy
from numba import njit


TidalPy.verbose_level = 0
TidalPy.logging_level = 0
TidalPy.use_disk = False


def test_integration():

    initial_conditions = (0., )
    @njit
    def simple_diffeq(time, variables):
        # dx_dt = 2x ->> x(t) = x^2 + 0

        return [2. * time]

    from TidalPy.integration_old.integration import ivp_integration

    _, integration_results_lsoda = ivp_integration(simple_diffeq, (0, 10.), initial_conditions,
                                                dependent_variable_names=('x',),
                                                integration_method='LSODA', integration_rtol=1e-4)
    _, integration_results_rk23 = ivp_integration(simple_diffeq, (0, 10.), initial_conditions,
                                               dependent_variable_names=('x',),
                                               integration_method='RK23', integration_rtol=1e-4)
    _, integration_results_rk45 = ivp_integration(simple_diffeq, (0, 10.), initial_conditions,
                                               dependent_variable_names=('x',),
                                               integration_method='RK45', integration_rtol=1e-4)

    np.testing.assert_array_almost_equal(integration_results_lsoda['x'], integration_results_lsoda['time_domain']**2)
    np.testing.assert_array_almost_equal(integration_results_rk23['x'], integration_results_rk23['time_domain']**2)
    np.testing.assert_array_almost_equal(integration_results_rk45['x'], integration_results_rk45['time_domain']**2)
