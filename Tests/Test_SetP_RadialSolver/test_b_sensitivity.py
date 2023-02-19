""" Test the `TidalPy.radial_solver.sensitivity` functionality. """
import pdb

import numpy as np
import pytest

import TidalPy

TidalPy.test_mode()

N = 5
radius_array = np.linspace(0., 1.e6, N)
radial_solutions = np.asarray(
    (
        (10. -5.j) * radius_array,  # y1
        (10. -5.j) * radius_array,  # y2
        (10. -5.j) * radius_array,  # y3
        (10. -5.j) * radius_array,  # y4
        (10. -5.j) * radius_array,  # y5
        (10. -5.j) * radius_array,  # y6
    ), dtype=np.complex128
)


def test_file_shear():
    """ Test the importing of the function from the root of the module. """

    from TidalPy.radial_solver import sensitivity_to_shear
    from TidalPy.radial_solver import sensitivity_to_bulk
    assert True


@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_sensitivity_to_shear_all_real(order_l):
    """ Test the sensitivity to shear modulus function. Test with both shear and bulk being real floats. """
    from TidalPy.radial_solver import sensitivity_to_shear

    # Make real float arrays for the viscoelastic properties
    shear_array = 1.0e10  * np.ones(N, dtype=np.float64)
    bulk_array  = 10.0e10 * np.ones(N, dtype=np.float64)

    shear_sensitivity = sensitivity_to_shear(radial_solutions, radius_array, shear_array, bulk_array, order_l=order_l)

    # Check shape
    assert shear_sensitivity.shape == (N,)

    # Check types
    # Note: Even though the shear and bulk are floats and real, the radial solutions are not. However, the sensitivity
    #  should always be real-valued.
    assert shear_sensitivity.dtype == np.float64

    # The shear sensitivity is not defined at r = 0. Check that it is nan.
    assert np.isnan(shear_sensitivity[0])

    # Check the last most value (which is just from a pre-calculated list)
    expected_last_val = [6141878208205017.0, 3.21367024371583e+16, 9.81228809938711e+16][order_l - 2]
    assert np.isclose(shear_sensitivity[-1], expected_last_val)

@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_sensitivity_to_shear_all_complex(order_l):
    """ Test the sensitivity to shear modulus function. Test with both shear and bulk being complex values. """
    from TidalPy.radial_solver import sensitivity_to_shear

    # Make real float arrays for the viscoelastic properties
    shear_array = (1.0e10 + 100j)   * np.ones(N, dtype=np.complex128)
    bulk_array  = (1.0e100 + 1000j) * np.ones(N, dtype=np.complex128)

    shear_sensitivity = sensitivity_to_shear(radial_solutions, radius_array, shear_array, bulk_array, order_l=order_l)

    # Check shape
    assert shear_sensitivity.shape == (N,)

    # Check types
    # Note: The sensitivity should always be real-valued even if bulk and shear are complex.
    assert shear_sensitivity.dtype == np.float64

    # The shear sensitivity is not defined at r = 0. Check that it is nan.
    assert np.isnan(shear_sensitivity[0])

    # Check the last most value (which is just from a pre-calculated list)
    expected_last_val = [7000000007500000.0, 3.7500000015e+16, 1.1550000002500002e+17][order_l - 2]
    assert np.isclose(shear_sensitivity[-1], expected_last_val)

@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_sensitivity_to_bulk_all_real(order_l):
    """ Test the sensitivity to bulk modulus function. Test with both shear and bulk being real floats. """
    from TidalPy.radial_solver import sensitivity_to_bulk

    # Make real float arrays for the viscoelastic properties
    shear_array = 1.0e10  * np.ones(N, dtype=np.float64)
    bulk_array  = 10.0e10 * np.ones(N, dtype=np.float64)

    bulk_sensitivity = sensitivity_to_bulk(radial_solutions, radius_array, shear_array, bulk_array, order_l=order_l)

    # Check shape
    assert bulk_sensitivity.shape == (N,)

    # Check types
    # Note: Even though the shear and bulk are floats and real, the radial solutions are not. However, the sensitivity
    #  should always be real-valued.
    assert bulk_sensitivity.dtype == np.float64

    # The shear sensitivity is not defined at r = 0. Check that it is nan.
    assert np.isnan(bulk_sensitivity[0])

    # Check the last most value (which is just from a pre-calculated list)
    expected_last_val = [2356408650528763.0, 1.8477526816618724e+16, 6.346716072665333e+16][order_l - 2]
    assert np.isclose(bulk_sensitivity[-1], expected_last_val)

@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_sensitivity_to_bulk_all_complex(order_l):
    """ Test the sensitivity to bulk modulus function. Test with both shear and bulk being complex values. """
    from TidalPy.radial_solver import sensitivity_to_bulk

    # Make real float arrays for the viscoelastic properties
    shear_array = (1.0e10 + 100j)   * np.ones(N, dtype=np.complex128)
    bulk_array  = (1.0e100 + 1000j) * np.ones(N, dtype=np.complex128)

    bulk_sensitivity = sensitivity_to_bulk(radial_solutions, radius_array, shear_array, bulk_array, order_l=order_l)

    # Check shape
    assert bulk_sensitivity.shape == (N,)

    # Check types
    # Note: The sensitivity should always be real-valued even if bulk and shear are complex.
    assert bulk_sensitivity.dtype == np.float64

    # The shear sensitivity is not defined at r = 0. Check that it is nan.
    assert np.isnan(bulk_sensitivity[0])

    # Check the last most value (which is just from a pre-calculated list)
    expected_last_val = [3000000000000000.5, 2.25e+16, 7.650000000000002e+16][order_l - 2]
    assert np.isclose(bulk_sensitivity[-1], expected_last_val)
