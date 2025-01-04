import pytest
import numpy as np

from TidalPy.exceptions import UnknownModelError, ArgumentException
from TidalPy.RadialSolver import radial_solver


def test_invalid_density_array_size():
    radius_array = np.linspace(0, 1000, 10, dtype=np.float64)
    density_array = np.linspace(1000, 2000, 9, dtype=np.float64)  # Incorrect size
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("solid",)
    is_static_bylayer = (True,)
    is_incompressible_bylayer = (True,)
    upper_radius_bylayer_array = np.array([1000.0], dtype=np.float64)
    
    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )

def test_invalid_upper_radius_bylayer_order():
    radius_array = np.concatenate((
        np.linspace(0, 500, 5, dtype=np.float64),
        np.linspace(500, 1000, 5, dtype=np.float64)
    ))
    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("solid", "liquid")
    is_static_bylayer = (True, False)
    is_incompressible_bylayer = (True, False)
    upper_radius_bylayer_array = np.array([1000.0, 500.0], dtype=np.float64)  # Incorrect order

    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )

def test_invalid_frequency_range():
    radius_array = np.linspace(0, 1000, 10, dtype=np.float64)
    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("solid",)
    is_static_bylayer = (False,)
    is_incompressible_bylayer = (False,)
    upper_radius_bylayer_array = np.array([1000.0], dtype=np.float64)

    with pytest.raises(ValueError):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1e-20, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )
    
    with pytest.raises(ValueError):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1e10, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )

def test_invalid_radius_array_start():
    radius_array = np.linspace(10, 1000, 10, dtype=np.float64)  # Starts at 10, not 0
    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("solid",)
    is_static_bylayer = (True,)
    is_incompressible_bylayer = (True,)
    upper_radius_bylayer_array = np.array([1000.0], dtype=np.float64)

    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )

def test_invalid_radius_array():
    radius_array = np.linspace(0, 1000, 10, dtype=np.float64)
    
    tmp = radius_array[4]
    radius_array[4] = radius_array[5]
    radius_array[5] = tmp

    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("solid",)
    is_static_bylayer = (True,)
    is_incompressible_bylayer = (True,)
    upper_radius_bylayer_array = np.array([1000.0], dtype=np.float64)

    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )

    radius_array = np.linspace(0, 1000, 10, dtype=np.float64)
    radius_array[4] = -radius_array[4]
    
    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )

def test_invalid_layer_type():
    radius_array = np.concatenate((
        np.linspace(0, 500, 5, dtype=np.float64),
        np.linspace(500, 1000, 5, dtype=np.float64)
    ))
    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("solid", "unknown")  # Invalid layer type
    is_static_bylayer = (False, False)
    is_incompressible_bylayer = (False, False)
    upper_radius_bylayer_array = np.array([500.0, 1000.0], dtype=np.float64)

    with pytest.raises(UnknownModelError):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )

def test_layer_missing_interface_value():
    radius_array = np.concatenate((
        np.linspace(0, 500, 5, dtype=np.float64),
        np.linspace(510, 1000, 5, dtype=np.float64) # Should start at 500
    ))
    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("liquid", "solid")
    is_static_bylayer = (False, False)
    is_incompressible_bylayer = (False, False)
    upper_radius_bylayer_array = np.array([500.0, 1000.0], dtype=np.float64)

    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )
    
    radius_array = np.concatenate((
        np.linspace(0, 490, 5, dtype=np.float64), # Should end at 500
        np.linspace(500, 1000, 5, dtype=np.float64)
    ))
    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )
    
    radius_array = np.concatenate((
        np.linspace(0, 490, 5, dtype=np.float64), # Should end at 500
        np.linspace(510, 1000, 5, dtype=np.float64) # Should start at 500
    ))
    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )

def test_layer_too_few_slices():
    radius_array = np.concatenate((
        np.linspace(0, 500, 2, dtype=np.float64),  # Needs to be at least 5
        np.linspace(510, 1000, 3, dtype=np.float64) # Needs to be at least 5
    ))
    density_array = np.linspace(1000, 2000, 5, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(5, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(5, dtype=np.complex128)
    layer_types = ("liquid", "solid")
    is_static_bylayer = (False, False)
    is_incompressible_bylayer = (False, False)
    upper_radius_bylayer_array = np.array([500.0, 1000.0], dtype=np.float64)

    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, raise_on_fail=True
        )
    
def test_prop_matrix_limitations_too_many_layers():
    radius_array = np.concatenate((
        np.linspace(0, 500, 5, dtype=np.float64),
        np.linspace(510, 1000, 5, dtype=np.float64)
    ))
    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("solid", "solid")
    is_static_bylayer = (True, True)
    is_incompressible_bylayer = (True, True)
    upper_radius_bylayer_array = np.array([500.0, 1000.0], dtype=np.float64)

    with pytest.raises(NotImplementedError):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, use_prop_matrix=True, raise_on_fail=True
        )

def test_prop_matrix_limitations_layer_assumptions():
    radius_array = np.linspace(0, 1000, 10, dtype=np.float64)
    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("liquid",)  # Must be solid for prop matrix
    is_static_bylayer = (True,)
    is_incompressible_bylayer = (True,)
    upper_radius_bylayer_array = np.array([1000.0], dtype=np.float64)

    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, use_prop_matrix=True, raise_on_fail=True
        )
    
    layer_types = ("solid",)
    is_static_bylayer = (False,) # Must be static for prop matrix
    is_incompressible_bylayer = (True,)
    upper_radius_bylayer_array = np.array([1000.0], dtype=np.float64)

    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, use_prop_matrix=True, raise_on_fail=True
        )
    
    layer_types = ("solid",)
    is_static_bylayer = (True,) 
    is_incompressible_bylayer = (False,)  # Must be incompressible for prop matrix
    upper_radius_bylayer_array = np.array([1000.0], dtype=np.float64)

    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array, use_prop_matrix=True, raise_on_fail=True
        )

def test_bad_starting_radius():
    radius_array = np.linspace(0, 1000, 10, dtype=np.float64)
    density_array = np.linspace(1000, 2000, 10, dtype=np.float64)
    complex_bulk_modulus_array = np.zeros(10, dtype=np.complex128)
    complex_shear_modulus_array = np.zeros(10, dtype=np.complex128)
    layer_types = ("solid",)
    is_static_bylayer = (True,)
    is_incompressible_bylayer = (True,)
    upper_radius_bylayer_array = np.array([1000.0], dtype=np.float64)
    
    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array,
            starting_radius = 0.91 * 1000,  # Must be less than 90% total radius.
            raise_on_fail=True
        )
    
    with pytest.raises(ArgumentException):
        radial_solver(
            radius_array, density_array, complex_bulk_modulus_array,
            complex_shear_modulus_array, 1.0, 3000.0, layer_types,
            is_static_bylayer, is_incompressible_bylayer,
            upper_radius_bylayer_array,
            starting_radius = 0.91 * 1000,  # Must be less than 90% total radius.
            raise_on_fail=True
        )
