import os
import pathlib

import pytest
import numpy as np

import TidalPy
TidalPy.test_mode = True


from TidalPy.constants import G
from TidalPy.RadialSolver.matrix_types.solid_matrix import fundamental_matrix

file_path = pathlib.Path(__file__).parent.resolve()

# Model planet - 2layers
N = 100
density_array = 5000. * np.ones(N)
radius_array = np.linspace(0., 1.e6, N+1)
radius_array_reduced = radius_array[1:]
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(N)])
gravity_array = G * mass_below / (radius_array[1:]**2)
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
    # Create empty arrays
    fundamental_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')
    inverse_fundamental_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')
    derivative_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')

    # For loop checks if there are any rare memory-related crashes
    for i in range(10):
        # Create empty arrays
        fundamental_matrix(
            radius_array_reduced,
            density_array,
            gravity_array,
            complex_shear_array,
            fundamental_mtx,
            inverse_fundamental_mtx,
            derivative_mtx,
            degree_l=degree_l)
    
    if degree_l == 2:
        assert np.allclose(fundamental_mtx, expected_fundamental_mtx)
        assert np.allclose(inverse_fundamental_mtx, expected_inverse_fundamental_mtx)
        assert np.allclose(derivative_mtx, expected_derivative_mtx)
    elif degree_l == 3:
        assert np.allclose(fundamental_mtx, expected_fundamental_mtx_l3)
        assert np.allclose(inverse_fundamental_mtx, expected_inverse_fundamental_mtx_l3)
        assert np.allclose(derivative_mtx, expected_derivative_mtx_l3)
    else:
        raise NotImplementedError

def test_solid_fundamental_matrix_errors():
    """ Test that the correct errors are raised when bad inputs are provided. """

    with pytest.raises(ValueError) as excinfo:
        fundamental_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')
        inverse_fundamental_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')
        derivative_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')

        fundamental_matrix(
            radius_array_reduced,
            density_array[2:],
            gravity_array,
            complex_shear_array,
            fundamental_mtx,
            inverse_fundamental_mtx,
            derivative_mtx,
            degree_l=2) 
        
    assert str(excinfo.value) == "Unexpected size encountered for density array."

    for j in range(2):
        for i in range(3):
            with pytest.raises(ValueError) as excinfo:
                fundamental_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')
                inverse_fundamental_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')
                derivative_mtx = np.empty((N, 6, 6), dtype=np.complex128, order='C')
                if i == 0:
                    if j == 0:
                        fundamental_mtx = np.empty((N-10, 6, 6), dtype=np.complex128, order='C')
                    elif j == 1:
                        fundamental_mtx = np.empty((N, 3, 6), dtype=np.complex128, order='C')
                elif i == 1:
                    if j == 0:
                        inverse_fundamental_mtx = np.empty((N-10, 6, 6), dtype=np.complex128, order='C')
                    elif j == 1:
                        inverse_fundamental_mtx = np.empty((N, 3, 6), dtype=np.complex128, order='C')
                elif i == 2:
                    if j == 0:
                        derivative_mtx = np.empty((N-10, 6, 6), dtype=np.complex128, order='C')
                    elif j == 1:
                        derivative_mtx = np.empty((N, 3, 6), dtype=np.complex128, order='C')

                fundamental_matrix(
                    radius_array_reduced,
                    density_array,
                    gravity_array,
                    complex_shear_array,
                    fundamental_mtx,
                    inverse_fundamental_mtx,
                    derivative_mtx,
                    degree_l=2) 
            
            if i == 0:
                assert str(excinfo.value) == "Unexpected shape encountered for Fundamental Matrix."
            elif i == 1:
                assert str(excinfo.value) == "Unexpected shape encountered for Fundamental Matrix (inv)."
            elif i == 2:
                assert str(excinfo.value) == "Unexpected shape encountered for Derivative Matrix."
