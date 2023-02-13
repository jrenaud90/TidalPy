""" Tests for calculating the derivatives of the radial functions throughout a liquid or solid layer.
"""

import numpy as np


import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.tides.multilayer.numerical_int.derivatives import (radial_derivatives_liquid_dynamic,
                                                                radial_derivatives_liquid_static,
                                                                radial_derivatives_solid_dynamic,
                                                                radial_derivatives_solid_static)

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex128)
bulk_array = 10.e10 * np.ones(10, dtype=np.float64)
radius_array_to_use = radius_array[1:]
frequency = 2. * np.pi / (86400. * 1.)


def test_derivatives_solid_static():
    """ Tests the calculation of radial derivatives for a solid layer using the static assumption """

    # Initialize 6 fake radial function solutions
    radial_funcs_complex = np.ones(6, dtype=np.complex128)
    radial_funcs = np.ones(12, dtype=np.float64)

    # Test for l=2
    solid_derivative = radial_derivatives_solid_static(
        radius_array_to_use[0], radial_funcs, shear_array[0], bulk_array[0],
        density_array[0], gravity_array[0], order_l=2
        )
    assert type(solid_derivative) == np.ndarray
    assert solid_derivative.size == 12
    assert solid_derivative.dtype == np.float64

    # Test for l=3
    solid_derivative = radial_derivatives_solid_static(
        radius_array_to_use[0], radial_funcs, shear_array[0], bulk_array[0],
        density_array[0], gravity_array[0], order_l=3
        )
    assert type(solid_derivative) == np.ndarray
    assert solid_derivative.size == 12
    assert solid_derivative.dtype == np.float64


def test_derivatives_solid_dynamic():
    """ Tests the calculation of radial derivatives for a solid layer using the dynamic assumption """

    # Initialize 6 fake radial function solutions
    radial_funcs_complex = np.ones(6, dtype=np.complex128)
    radial_funcs = np.ones(12, dtype=np.float64)

    # Test for l=2
    solid_derivative = radial_derivatives_solid_dynamic(
        radius_array_to_use[0], radial_funcs, shear_array[0], bulk_array[0],
        density_array[0], gravity_array[0], frequency, order_l=2
        )
    assert type(solid_derivative) == np.ndarray
    assert solid_derivative.size == 12
    assert solid_derivative.dtype == np.float64

    # Test for l=3
    solid_derivative = radial_derivatives_solid_dynamic(
        radius_array_to_use[0], radial_funcs, shear_array[0], bulk_array[0],
        density_array[0], gravity_array[0], frequency, order_l=3
        )
    assert type(solid_derivative) == np.ndarray
    assert solid_derivative.size == 12
    assert solid_derivative.dtype == np.float64


def test_derivatives_liquid_static():
    """ Tests the calculation of radial derivatives for a liquid layer using the static assumption """

    # Initialize 2 fake radial function solutions
    radial_funcs_complex = np.ones(2, dtype=np.complex128)
    radial_funcs = np.ones(4, dtype=np.float64)

    # Test for l=2
    liquid_derivative = radial_derivatives_liquid_static(
        radius_array_to_use[0], radial_funcs, density_array[0], gravity_array[0],
        order_l=2
        )
    assert type(liquid_derivative) == np.ndarray
    assert liquid_derivative.size == 4
    assert liquid_derivative.dtype == np.float64

    # Test for l=3
    liquid_derivative = radial_derivatives_liquid_static(
        radius_array_to_use[0], radial_funcs, density_array[0], gravity_array[0],
        order_l=3
        )
    assert type(liquid_derivative) == np.ndarray
    assert liquid_derivative.size == 4
    assert liquid_derivative.dtype == np.float64


def test_derivatives_liquid_dynamic():
    """ Tests the calculation of radial derivatives for a liquid layer using the dynamic assumption """

    # Initialize 4 fake radial function solutions
    radial_funcs_complex = np.ones(4, dtype=np.complex128)
    radial_funcs = np.ones(8, dtype=np.complex128)

    # Test for l=2
    liquid_derivative = radial_derivatives_liquid_dynamic(
        radius_array_to_use[0], radial_funcs, bulk_array[0], density_array[0],
        gravity_array[0], frequency, order_l=2
        )
    assert type(liquid_derivative) == np.ndarray
    assert liquid_derivative.size == 8
    assert liquid_derivative.dtype == np.float64

    # Test for l=3
    liquid_derivative = radial_derivatives_liquid_dynamic(
        radius_array_to_use[0], radial_funcs, bulk_array[0],
        density_array[0], gravity_array[0], frequency, order_l=3
        )
    assert type(liquid_derivative) == np.ndarray
    assert liquid_derivative.size == 8
    assert liquid_derivative.dtype == np.float64
