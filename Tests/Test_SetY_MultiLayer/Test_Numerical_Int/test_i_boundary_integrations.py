""" Tests for calculating the radial solution across the interior of a planet
"""

import numpy as np

import TidalPy
from TidalPy.constants import G
from TidalPy.toolbox.multilayer import calculate_sls, calculate_ssls, calculate_homogen_solid

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
R = radius_array[-1]
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex128)
bulk_array = 10.e10 * np.ones(10, dtype=np.float64)
radius_array_to_use = radius_array[1:]
frequency = 2. * np.pi / (86400. * 1.)


def test_calculate_homogen():
    """ Test the solution calculation for homogeneous planet """

    # Test static
    tidal_y, tidal_y_derivative = calculate_homogen_solid(
        radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
        use_static=True,
        use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
        julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]

    # Test dynamic
    tidal_y, tidal_y_derivative = calculate_homogen_solid(
            radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
            use_static=False,
            use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
            julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]

def test_calculate_sls():
    """ Test the solution calculation for solid-liquid-solid planet structure """

    # Test all static
    tidal_y, tidal_y_derivative = calculate_sls(
        radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
        interface_1_radius=(1. / 3.) * R, interface_2_radius=(2. / 3.) * R,
        layer_0_static=True, layer_1_static=True, layer_2_static=True, order_l=2, use_kamata=True,
        use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
        julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]

    # Test all dynamic
    tidal_y, tidal_y_derivative = calculate_sls(
            radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
            interface_1_radius=(1. / 3.) * R, interface_2_radius=(2. / 3.) * R,
            layer_0_static=False, layer_1_static=False, layer_2_static=False, order_l=2, use_kamata=True,
            use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
            julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]

    # Test mix liq=Static
    tidal_y, tidal_y_derivative = calculate_sls(
            radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
            interface_1_radius=(1. / 3.) * R, interface_2_radius=(2. / 3.) * R,
            layer_0_static=False, layer_1_static=True, layer_2_static=False, order_l=2, use_kamata=True,
            use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
            julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]

    # Test mix liq=Dynamic
    tidal_y, tidal_y_derivative = calculate_sls(
            radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
            interface_1_radius=(1. / 3.) * R, interface_2_radius=(2. / 3.) * R,
            layer_0_static=True, layer_1_static=False, layer_2_static=True, order_l=2, use_kamata=True,
            use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
            julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]


def test_calculate_ssls():
    """ Test the solution calculation for solid-solid-liquid-solid planet structure """

    # Test all static
    tidal_y, tidal_y_derivative = calculate_ssls(
        radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
        interface_1_radius=(1. / 4.) * R, interface_2_radius=(2. / 4.) * R, interface_3_radius=(3. / 4.) * R,
        layer_0_static=True, layer_1_static=True, layer_2_static=True, layer_3_static=True, order_l=2, use_kamata=True,
        use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
        julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]

    # Test all dynamic
    tidal_y, tidal_y_derivative = calculate_ssls(
            radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
            interface_1_radius=(1. / 4.) * R, interface_2_radius=(2. / 4.) * R, interface_3_radius=(3. / 4.) * R,
            layer_0_static=False, layer_1_static=False, layer_2_static=False, layer_3_static=False, order_l=2,
            use_kamata=True,
            use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
            julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]

    # Test mix liq=Static
    tidal_y, tidal_y_derivative = calculate_ssls(
            radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
            interface_1_radius=(1. / 4.) * R, interface_2_radius=(2. / 4.) * R, interface_3_radius=(3. / 4.) * R,
            layer_0_static=False, layer_1_static=False, layer_2_static=True, layer_3_static=False, order_l=2,
            use_kamata=True,
            use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
            julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]

    # Test mix liq=Dynamic
    tidal_y, tidal_y_derivative = calculate_ssls(
            radius_array_to_use, shear_array, bulk_array, density_array, gravity_array, frequency,
            interface_1_radius=(1. / 4.) * R, interface_2_radius=(2. / 4.) * R, interface_3_radius=(3. / 4.) * R,
            layer_0_static=True, layer_1_static=True, layer_2_static=False, layer_3_static=True, order_l=2,
            use_kamata=True,
            use_julia=False, verbose=False, int_rtol=1.0e-6, int_atol=1.0e-4, scipy_int_method='RK45',
            julia_int_method='Tsit5')

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex, np.complex128, complex]
    assert tidal_y_derivative.shape == (6, 10)
    assert tidal_y_derivative.dtype in [np.complex, np.complex128, complex]
