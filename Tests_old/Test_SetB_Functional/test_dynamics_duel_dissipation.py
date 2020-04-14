import numpy as np
from scipy.constants import G

import TidalPy


TidalPy.verbose_level = 0
TidalPy.logging_level = 0
TidalPy.use_disk = False

from TidalPy.dynamics import semi_major_axis_derivative_duel


MASS_1 = 5.972e24
MASS_2 = 7.3459e22
SEMI_A = 3.844e8
TIDAL_HEAT_1 = 1.e10
TIDAL_HEAT_2 = 1.e14
TIDAL_TORQUE_1 = 1.e15
TIDAL_TORQUE_2 = -1.e20
SPIN_1 = 2. * np.pi / (1.0 * 24 * 60 * 60)
SPIN_2 = 2. * np.pi / (27.322 * 24 * 60 * 60)


def test_semi_major_axis_derivative():
    # Test for correct logic
    np.testing.assert_approx_equal(semi_major_axis_derivative_duel(G**(1 / 2), 1., 1., 1., 1., 1., 1., 1., 1.), -8.)
    np.testing.assert_approx_equal(semi_major_axis_derivative_duel(-G**(1 / 2), 1., 1., 1., 1., 1., 1., 1., 1.), -8.)
    np.testing.assert_approx_equal(semi_major_axis_derivative_duel(G**(1 / 2), 1., 1., 1., 0., 0., 1., 0., 0.), 0.)
    np.testing.assert_approx_equal(semi_major_axis_derivative_duel(G**(1 / 2), 1., 1., 1., -1., -1., 1., -1., -1.), 8.)
    np.testing.assert_approx_equal(semi_major_axis_derivative_duel(G**(1 / 2), 1., 1., 1., -1., -1., 1., 1., 1.), 0.)

    # Test for correct calculations
    np.testing.assert_approx_equal(
            semi_major_axis_derivative_duel(SEMI_A, MASS_1, MASS_2, SPIN_1, TIDAL_TORQUE_1, TIDAL_HEAT_1,
                                            SPIN_2, TIDAL_TORQUE_2, TIDAL_HEAT_2),
            1.6763632928138172e-06)
    np.testing.assert_approx_equal(
            semi_major_axis_derivative_duel(SEMI_A, MASS_1, MASS_2, SPIN_1, TIDAL_TORQUE_1, TIDAL_HEAT_1,
                                            -SPIN_2, TIDAL_TORQUE_2, TIDAL_HEAT_2),
            -3.696727374620553e-06)
    np.testing.assert_approx_equal(
            semi_major_axis_derivative_duel(SEMI_A, MASS_1, MASS_2, -SPIN_1, TIDAL_TORQUE_1, TIDAL_HEAT_1,
                                            -SPIN_2, TIDAL_TORQUE_2, TIDAL_HEAT_2),
            -3.6952593387883963e-06)
