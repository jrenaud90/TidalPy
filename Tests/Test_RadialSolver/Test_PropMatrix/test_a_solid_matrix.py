import os
import pathlib

import pytest
import numpy as np

from TidalPy.constants import G
from TidalPy.RadialSolver.matrix_types.solid_matrix import fundamental_matrix

file_path = pathlib.Path(__file__).parent.resolve()

# Model planet - 2layers
N = 100
# The pre-calculated arrays are based on a radius that starts at r>0; New method starts at r=0, so we want to use one extra N and then compare the two for 1:N+1
N += 1

density_array = 5000. * np.ones(N)
radius_array = np.linspace(0., 1.e6, N)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array[1:]
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(N-1)])
gravity_array = G * mass_below / (radius_array[1:]**2)
gravity_array = np.insert(gravity_array, 0, 0.0, axis=0)
complex_shear_array = 5.e10 * np.ones(N, dtype=np.complex128)


expected_fundamental_mtx = np.load(os.path.join(file_path, 'solid_fundamental_matrix.npy'))
expected_inverse_fundamental_mtx = np.load(os.path.join(file_path, 'solid_inverse_fundamental_matrix.npy'))
expected_derivative_mtx = np.load(os.path.join(file_path, 'solid_derivative_matrix.npy'))

expected_fundamental_mtx_l3 = np.load(os.path.join(file_path, 'solid_fundamental_matrix_l3.npy'))
expected_inverse_fundamental_mtx_l3 = np.load(os.path.join(file_path, 'solid_inverse_fundamental_matrix_l3.npy'))
expected_derivative_mtx_l3 = np.load(os.path.join(file_path, 'solid_derivative_matrix_l3.npy'))

@pytest.mark.parametrize('degree_l', (2, 3))
def test_solid_fundamental_matrix(degree_l):
    """ Test the functions to build the fundamental matrices for a solid spherical shell. """

    # For loop checks if there are any rare memory-related crashes
    for i in range(10):
        fundamental_mtx, inverse_fundamental_mtx, derivative_mtx = fundamental_matrix(
            radius_array,
            density_array,
            gravity_array,
            complex_shear_array,
            degree_l=degree_l)
    
    # When comparing to the saved matrices, start at index 1 for the new solution since the saved solution was not defined at r=0
    if degree_l == 2:
        np.testing.assert_allclose(fundamental_mtx[1:, :, :], expected_fundamental_mtx)
        np.testing.assert_allclose(inverse_fundamental_mtx[1:, :, :], expected_inverse_fundamental_mtx)
        np.testing.assert_allclose(derivative_mtx[1:, :, :], expected_derivative_mtx)
    elif degree_l == 3:
        np.testing.assert_allclose(fundamental_mtx[1:, :, :], expected_fundamental_mtx_l3)
        np.testing.assert_allclose(inverse_fundamental_mtx[1:, :, :], expected_inverse_fundamental_mtx_l3)
        np.testing.assert_allclose(derivative_mtx[1:, :, :], expected_derivative_mtx_l3)
    else:
        raise NotImplementedError

def test_solid_fundamental_matrix_errors():
    """ Test that the correct errors are raised when bad inputs are provided. """

    with pytest.raises(ValueError) as excinfo:
        fundamental_mtx, inverse_fundamental_mtx, derivative_mtx = fundamental_matrix(
            radius_array,
            density_array[2:],
            gravity_array,
            complex_shear_array,
            degree_l=2)

    assert str(excinfo.value) == "Unexpected size encountered for density array."
