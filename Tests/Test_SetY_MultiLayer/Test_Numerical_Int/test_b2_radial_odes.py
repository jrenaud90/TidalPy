""" Tests for calculating the derivatives of the radial functions throughout a liquid or solid layer.
"""

import numpy as np

import TidalPy

TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.radial_solver.numerical.derivatives import (
    dynamic_liquid_ode, dynamic_solid_ode, static_liquid_ode, static_solid_ode)

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
    """ Tests the calculation of radial ODE for a solid layer using the static assumption """

    # Initialize 6 (increased to 12 for complex->float) fake radial function solutions
    radial_funcs = (10. + 1.0j) * np.ones(12, dtype=np.float64)

    # Test for l=2
    solid_derivative = static_solid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, shear_array, bulk_array,
            density_array, gravity_array, order_l=2
            )
    assert type(solid_derivative) == np.ndarray
    assert len(solid_derivative) == 12
    for y_i in range(12):
        assert type(solid_derivative[y_i]) in [float, np.float64]

    # Test for l=3
    solid_derivative = static_solid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, shear_array, bulk_array,
            density_array, gravity_array, order_l=3
            )
    assert type(solid_derivative) == np.ndarray
    assert len(solid_derivative) == 12
    for y_i in range(12):
        assert type(solid_derivative[y_i]) in [float, np.float64]


def test_derivatives_solid_static_incomp():
    """ Tests the calculation of radial ODE for a solid layer using the static and incompressible assumption """

    # Initialize 6 (increased to 12 for complex->float) fake radial function solutions
    radial_funcs = (10. + 1.0j) * np.ones(12, dtype=np.float64)

    # Test for l=2
    solid_derivative = static_solid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, shear_array, bulk_array,
            density_array, gravity_array, order_l=2, incompressible=True
            )
    assert type(solid_derivative) == np.ndarray
    assert len(solid_derivative) == 12
    for y_i in range(12):
        assert type(solid_derivative[y_i]) in [float, np.float64]

    # Test for l=3
    solid_derivative = static_solid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, shear_array, bulk_array,
            density_array, gravity_array, order_l=3, incompressible=True
            )
    assert type(solid_derivative) == np.ndarray
    assert len(solid_derivative) == 12
    for y_i in range(12):
        assert type(solid_derivative[y_i]) in [float, np.float64]


def test_derivatives_solid_dynamic():
    """ Tests the calculation of radial derivatives for a solid layer using the dynamic assumption """

    # Initialize 6 (increased to 12 for complex->float) fake radial function solutions
    radial_funcs = (10. + 1.0j) * np.ones(12, dtype=np.float64)

    # Test for l=2
    solid_derivative = dynamic_solid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, shear_array, bulk_array,
            density_array, gravity_array, frequency, order_l=2
            )
    assert type(solid_derivative) == np.ndarray
    assert len(solid_derivative) == 12
    for y_i in range(12):
        assert type(solid_derivative[y_i]) in [float, np.float64]

    # Test for l=3
    solid_derivative = dynamic_solid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, shear_array, bulk_array,
            density_array, gravity_array, frequency, order_l=3
            )
    assert type(solid_derivative) == np.ndarray
    assert len(solid_derivative) == 12
    for y_i in range(12):
        assert type(solid_derivative[y_i]) in [float, np.float64]


def test_derivatives_solid_dynamic_incomp():
    """ Tests the calculation of radial derivatives for a solid layer using the dynamic and incompressible assumption """

    # Initialize 6 (increased to 12 for complex->float) fake radial function solutions
    radial_funcs = (10. + 1.0j) * np.ones(12, dtype=np.float64)

    # Test for l=2
    solid_derivative = dynamic_solid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, shear_array, bulk_array,
            density_array, gravity_array, frequency, order_l=2, incompressible=True
            )
    assert type(solid_derivative) == np.ndarray
    assert len(solid_derivative) == 12
    for y_i in range(12):
        assert type(solid_derivative[y_i]) in [float, np.float64]

    # Test for l=3
    solid_derivative = dynamic_solid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, shear_array, bulk_array,
            density_array, gravity_array, frequency, order_l=3, incompressible=True
            )
    assert type(solid_derivative) == np.ndarray
    assert len(solid_derivative) == 12
    for y_i in range(12):
        assert type(solid_derivative[y_i]) in [float, np.float64]


def test_derivatives_liquid_static():
    """ Tests the calculation of radial derivatives for a liquid layer using the static assumption """

    # Initialize 2 (increased to 4 for complex->float) fake radial function solutions
    radial_funcs = (10. + 1.0j) * np.ones(4, dtype=np.float64)

    # Test for l=2
    liquid_derivative = static_liquid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, density_array, gravity_array,
            order_l=2
            )
    assert type(liquid_derivative) == np.ndarray
    assert len(liquid_derivative) == 4
    for y_i in range(4):
        assert type(liquid_derivative[y_i]) in [float, np.float64]

    # Test for l=3
    liquid_derivative = static_liquid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, density_array, gravity_array,
            order_l=3
            )
    assert type(liquid_derivative) == np.ndarray
    assert len(liquid_derivative) == 4
    for y_i in range(4):
        assert type(liquid_derivative[y_i]) in [float, np.float64]


def test_derivatives_liquid_static_incomp():
    """ Tests the calculation of radial derivatives for a liquid layer using the static and incompressible assumption """

    # Initialize 2 (increased to 4 for complex->float) fake radial function solutions
    radial_funcs = (10. + 1.0j) * np.ones(4, dtype=np.float64)

    # Test for l=2
    liquid_derivative = static_liquid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, density_array, gravity_array,
            order_l=2, incompressible=True
            )
    assert type(liquid_derivative) == np.ndarray
    assert len(liquid_derivative) == 4
    for y_i in range(4):
        assert type(liquid_derivative[y_i]) in [float, np.float64]

    # Test for l=3
    liquid_derivative = static_liquid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, density_array, gravity_array,
            order_l=3, incompressible=True
            )
    assert type(liquid_derivative) == np.ndarray
    assert len(liquid_derivative) == 4
    for y_i in range(4):
        assert type(liquid_derivative[y_i]) in [float, np.float64]


def test_derivatives_liquid_dynamic():
    """ Tests the calculation of radial derivatives for a liquid layer using the dynamic assumption """

    # Initialize 4 (increased to 8 for complex->float) fake radial function solutions
    radial_funcs = (10. + 1.0j) * np.ones(8, dtype=np.float64)

    # Test for l=2
    liquid_derivative = dynamic_liquid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, bulk_array, density_array,
            gravity_array, frequency, order_l=2
            )
    assert type(liquid_derivative) == np.ndarray
    assert len(liquid_derivative) == 8
    for y_i in range(8):
        assert type(liquid_derivative[y_i]) in [float, np.float64]

    # Test for l=3
    liquid_derivative = dynamic_liquid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, bulk_array,
            density_array, gravity_array, frequency, order_l=3
            )
    assert type(liquid_derivative) == np.ndarray
    assert len(liquid_derivative) == 8
    for y_i in range(8):
        assert type(liquid_derivative[y_i]) in [float, np.float64]


def test_derivatives_liquid_dynamic_incomp():
    """ Tests the calculation of radial derivatives for a liquid layer using the dynamic and incompressible assumption """

    # Initialize 4 (increased to 8 for complex->float) fake radial function solutions
    radial_funcs = (10. + 1.0j) * np.ones(8, dtype=np.float64)

    # Test for l=2
    liquid_derivative = dynamic_liquid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, bulk_array, density_array,
            gravity_array, frequency, order_l=2, incompressible=True
            )
    assert type(liquid_derivative) == np.ndarray
    assert len(liquid_derivative) == 8
    for y_i in range(8):
        assert type(liquid_derivative[y_i]) in [float, np.float64]

    # Test for l=3
    liquid_derivative = dynamic_liquid_ode(
            radius_array_to_use[2], radial_funcs,
            radius_array_to_use, bulk_array,
            density_array, gravity_array, frequency, order_l=3, incompressible=True
            )
    assert type(liquid_derivative) == np.ndarray
    assert len(liquid_derivative) == 8
    for y_i in range(8):
        assert type(liquid_derivative[y_i]) in [float, np.float64]
