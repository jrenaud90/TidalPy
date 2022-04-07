""" Tests for propagating the multilayer dissipation model using the fundamental matrix and its inverse
"""

import numpy as np

import TidalPy
from TidalPy.constants import G
from TidalPy.tides.multilayer.matrix.fundamental_solid import fundamental_matrix_generic, fundamental_matrix_orderl2
from TidalPy.tides.multilayer.matrix.propagate import propagate

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex128)


def test_calc_fundamental_order2():
    # Calculate the fundamental matrix and its inverse
    F, F_inv, deriv_mtx = fundamental_matrix_orderl2(
        radius_array[1:], shear_array,
        density_array, gravity_array
        )

    # Central boundary condition
    ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
    core_condition = np.zeros((6, 3), dtype=np.complex128)
    # Roberts & Nimmo (2000): Liquid innermost zone.
    core_condition[2, 0] = 1.0
    core_condition[3, 1] = 1.0
    core_condition[5, 2] = 1.0

    # Find tidal solution
    tidal_y, tidal_y_derivative = \
        propagate(F, F_inv, deriv_mtx, core_condition, world_radius=radius_array[-1], order_l=2)

    # See if shape matches expectations
    assert tidal_y.shape[0] == 6
    assert tidal_y.shape[1] == 10
    assert tidal_y_derivative.shape[0] == 6
    assert tidal_y_derivative.shape[1] == 10

    # See if the types make sense
    for i in range(6):
        assert type(tidal_y[i, 0]) in [np.complex128, complex]
        assert type(tidal_y_derivative[i, 0]) in [np.complex128, complex]


def test_calc_fundamental_order3():
    # Calculate the fundamental matrix and its inverse
    F, F_inv, deriv_mtx = fundamental_matrix_generic(
        radius_array[1:], shear_array,
        density_array, gravity_array, order_l=3
        )

    # Central boundary condition
    ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
    core_condition = np.zeros((6, 3), dtype=np.complex128)
    # Roberts & Nimmo (2000): Liquid innermost zone.
    core_condition[2, 0] = 1.0
    core_condition[3, 1] = 1.0
    core_condition[5, 2] = 1.0

    # Find tidal solution
    tidal_y, tidal_y_derivative = \
        propagate(F, F_inv, deriv_mtx, core_condition, world_radius=radius_array[-1], order_l=3)

    # See if shape matches expectations
    assert tidal_y.shape[0] == 6
    assert tidal_y.shape[1] == 10
    assert tidal_y_derivative.shape[0] == 6
    assert tidal_y_derivative.shape[1] == 10

    # See if the types make sense
    for i in range(6):
        assert type(tidal_y[i, 0]) in [np.complex128, complex]
        assert type(tidal_y_derivative[i, 0]) in [np.complex128, complex]
