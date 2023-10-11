""" Test the `TidalPy.radial_solver.matrix.propagate` functionality. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.radial_solver.matrix import fundamental_matrix_generic
from TidalPy.radial_solver.matrix import matrix_propagate

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
radius_array_reduced = radius_array[1:]
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
complex_shear_array = 5.e10 * np.ones(10, dtype=np.complex128)

@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_matrix_propagate(order_l):
    """ Test `matrix_propagate` function for multiple ls. """

    # Calculate the fundamental matrix and its inverse
    fundamental_mtx, inverse_fundamental_mtx, derivative_mtx = \
        fundamental_matrix_generic(
            radius_array_reduced, complex_shear_array, density_array, gravity_array,
            order_l=order_l)

    # Central boundary condition
    ## From IcyDwarf: "They are inconsequential on the rest of the solution, so false assumptions are OK."
    core_condition = np.zeros((6, 3), dtype=np.complex128)
    # Roberts & Nimmo (2000): Liquid innermost zone.
    core_condition[2, 0] = 1.0
    core_condition[3, 1] = 1.0
    core_condition[5, 2] = 1.0

    # Find tidal solution
    radial_solutions = \
        matrix_propagate(fundamental_mtx, inverse_fundamental_mtx, derivative_mtx,
                         core_condition, world_radius=radius_array[-1], order_l=order_l)

    # See if shape matches expectations
    assert radial_solutions.shape[0] == 6
    assert radial_solutions.shape[1] == 10

    # See if the types make sense
    for i in range(6):
        assert type(radial_solutions[i, 0]) in [np.complex128, complex]
