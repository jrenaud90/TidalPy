""" Test the `TidalPy.radial_solver.matrix.fundamental_solid` functionality. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.radial_solver.matrix.fundamental_solid import fundamental_matrix_generic, fundamental_matrix_orderl2

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

def test_fundamental_matrix_orderl2():
    """ Test `fundamental_matrix_orderl2` function. """

    fundamental_mtx, inverse_fundamental_mtx, derivative_mtx = \
        fundamental_matrix_orderl2(radius_array_reduced, complex_shear_array, density_array, gravity_array)

    # Check that the shapes are correct
    assert fundamental_mtx.shape[0] == 6
    assert fundamental_mtx.shape[1] == 6
    assert fundamental_mtx.shape[2] == 10

    assert inverse_fundamental_mtx.shape[0] == 6
    assert inverse_fundamental_mtx.shape[1] == 6
    assert inverse_fundamental_mtx.shape[2] == 10

    assert derivative_mtx.shape[0] == 6
    assert derivative_mtx.shape[1] == 6
    assert derivative_mtx.shape[2] == 10

    # See if the types make sense
    assert type(fundamental_mtx[0, 0, 0]) in [np.complex128, complex]
    assert type(inverse_fundamental_mtx[0, 0, 0]) in [np.complex128, complex]
    assert type(derivative_mtx[0, 0, 0]) in [np.complex128, complex]

    # Check that the inverse is correct
    for i in range(10):
        identity = fundamental_mtx[:, :, i] @ inverse_fundamental_mtx[:, :, i]
        # TODO: The generic version fails if atol is < 1e-6; is this an acceptable error level?
        assert np.allclose(identity, np.identity(6), atol=1.e-6)

@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_fundamental_matrix_generic(order_l):
    """ Test `fundamental_matrix_generic` function using multiple ls. """

    fundamental_mtx, inverse_fundamental_mtx, derivative_mtx = \
        fundamental_matrix_generic(radius_array_reduced, complex_shear_array, density_array, gravity_array,
                                   order_l=order_l)

    # Check that the shapes are correct
    assert fundamental_mtx.shape[0] == 6
    assert fundamental_mtx.shape[1] == 6
    assert fundamental_mtx.shape[2] == 10

    assert inverse_fundamental_mtx.shape[0] == 6
    assert inverse_fundamental_mtx.shape[1] == 6
    assert inverse_fundamental_mtx.shape[2] == 10

    assert derivative_mtx.shape[0] == 6
    assert derivative_mtx.shape[1] == 6
    assert derivative_mtx.shape[2] == 10

    # See if the types make sense
    assert type(fundamental_mtx[0, 0, 0]) in [np.complex128, complex]
    assert type(inverse_fundamental_mtx[0, 0, 0]) in [np.complex128, complex]
    assert type(derivative_mtx[0, 0, 0]) in [np.complex128, complex]

    # Check that the inverse is correct
    for i in range(10):
        identity = fundamental_mtx[:, :, i] @ inverse_fundamental_mtx[:, :, i]
        # TODO: The generic version fails if atol is < 1e-6; is this an acceptable error level?
        assert np.allclose(identity, np.identity(6), atol=1.e-6)

def test_compare_generic_to_fixed():
    """ Compare the generic to the l=2 version of the matrix. """

    # Calculate generic
    fundamental_mtx, inverse_fundamental_mtx, derivative_mtx = \
        fundamental_matrix_generic(radius_array_reduced, complex_shear_array, density_array, gravity_array,
                                   order_l=2)

    # Calculate l=2
    fundamental_mtx_l2, inverse_fundamental_mtx_l2, derivative_mtx_l2 = \
        fundamental_matrix_orderl2(radius_array_reduced, complex_shear_array, density_array, gravity_array)

    # Check that the values match the hardcoded l=2 version
    for i in range(10):
        assert np.allclose(fundamental_mtx[:, :, i],         fundamental_mtx_l2[:, :, i])
        assert np.allclose(inverse_fundamental_mtx[:, :, i], inverse_fundamental_mtx_l2[:, :, i])
        assert np.allclose(derivative_mtx[:, :, i],          derivative_mtx_l2[:, :, i])
