""" Tests for calculating the derivatives of the radial functions throughout a liquid or solid layer.
"""

import numpy as np

import TidalPy
from TidalPy.tides.multilayer.numerical_int import radial_derivatives_solid_static, radial_derivatives_solid_dynamic, \
    radial_derivatives_liquid_static, radial_derivatives_liquid_dynamic
from TidalPy.constants import G

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
bulk_array = 10.e10 * np.ones(10, dtype=np.float64)
radius_array_to_use = radius_array[1:]
frequency = 2. * np.pi / (86400. * 1.)


def test_derivatives_solid_static():
    """ Tests the calculation of radial derivatives for a solid layer using the static assumption """

    # Initialize 6 fake radial function solutions
    radial_funcs = tuple([
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype)
    ])

    # Test for l=2
    solid_derivative = radial_derivatives_solid_static(radius_array_to_use, radial_funcs, shear_array, bulk_array,
                                                       density_array, gravity_array, order_l=2)
    assert type(solid_derivative) == tuple
    assert len(solid_derivative) == 6
    for y_i in range(6):
        assert solid_derivative[y_i].dtype == shear_array.dtype
        assert solid_derivative[y_i].shape == (10,)

    # Test for l=3
    solid_derivative = radial_derivatives_solid_static(radius_array_to_use, radial_funcs, shear_array, bulk_array,
                                                       density_array, gravity_array, order_l=3)
    assert type(solid_derivative) == tuple
    assert len(solid_derivative) == 6
    for y_i in range(6):
        assert solid_derivative[y_i].dtype == shear_array.dtype
        assert solid_derivative[y_i].shape == (10,)


def test_derivatives_solid_dynamic():
    """ Tests the calculation of radial derivatives for a solid layer using the dynamic assumption """

    # Initialize 6 fake radial function solutions
    radial_funcs = tuple([
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype),
        np.ones(10, dtype=shear_array.dtype)
    ])

    # Test for l=2
    solid_derivative = radial_derivatives_solid_dynamic(radius_array_to_use, radial_funcs, shear_array, bulk_array,
                                                        density_array, gravity_array, frequency, order_l=2)
    assert type(solid_derivative) == tuple
    assert len(solid_derivative) == 6
    for y_i in range(6):
        assert solid_derivative[y_i].dtype == shear_array.dtype
        assert solid_derivative[y_i].shape == (10,)

    # Test for l=3
    solid_derivative = radial_derivatives_solid_dynamic(radius_array_to_use, radial_funcs, shear_array, bulk_array,
                                                        density_array, gravity_array, frequency, order_l=3)
    assert type(solid_derivative) == tuple
    assert len(solid_derivative) == 6
    for y_i in range(6):
        assert solid_derivative[y_i].dtype == shear_array.dtype
        assert solid_derivative[y_i].shape == (10,)

    # Test for an array in the frequency variable
    # Initialize 6 fake radial function solutions
    radial_funcs_mtx = tuple([
        np.ones((20, 10), dtype=shear_array.dtype),
        np.ones((20, 10), dtype=shear_array.dtype),
        np.ones((20, 10), dtype=shear_array.dtype),
        np.ones((20, 10), dtype=shear_array.dtype),
        np.ones((20, 10), dtype=shear_array.dtype),
        np.ones((20, 10), dtype=shear_array.dtype)
    ])
    freq_domain = np.linspace(-1., 1., 20)
    rad_mtx, freq_mtx = np.meshgrid(radius_array_to_use, freq_domain)
    shear_mtx, _ = np.meshgrid(shear_array, freq_domain)
    bulk_mtx, _ = np.meshgrid(bulk_array, freq_domain)
    den_mtx, _ = np.meshgrid(density_array, freq_domain)
    grav_mtx, _ = np.meshgrid(gravity_array, freq_domain)
    solid_guess = radial_derivatives_solid_dynamic(rad_mtx, radial_funcs_mtx, shear_mtx, bulk_mtx,
                                                   den_mtx, grav_mtx, freq_mtx, order_l=3)
    assert type(solid_guess) == tuple
    assert len(solid_guess) == 6
    for sn in range(6):
        assert type(solid_guess[sn]) == np.ndarray
        assert solid_guess[sn].dtype == shear_mtx.dtype
        # Dynamic solid solution should have 6 y values across the radius domain (10)
        #     and again across the new freq domain (20)
        assert solid_guess[sn].shape == (20, 10)


def test_derivatives_liquid_static():
    """ Tests the calculation of radial derivatives for a liquid layer using the static assumption """

    # Initialize 2 fake radial function solutions
    radial_funcs = tuple([
        np.ones(10, dtype=bulk_array.dtype),
        np.ones(10, dtype=bulk_array.dtype)
    ])

    # Test for l=2
    solid_derivative = radial_derivatives_liquid_static(radius_array_to_use, radial_funcs, density_array, gravity_array,
                                                        order_l=2)
    assert type(solid_derivative) == tuple
    assert len(solid_derivative) == 2
    for y_i in range(2):
        assert solid_derivative[y_i].dtype == bulk_array.dtype
        assert solid_derivative[y_i].shape == (10,)

    # Test for l=3
    solid_derivative = radial_derivatives_liquid_static(radius_array_to_use, radial_funcs, density_array, gravity_array,
                                                        order_l=3)
    assert type(solid_derivative) == tuple
    assert len(solid_derivative) == 2
    for y_i in range(2):
        assert solid_derivative[y_i].dtype == bulk_array.dtype
        assert solid_derivative[y_i].shape == (10,)


def test_derivatives_liquid_dynamic():
    """ Tests the calculation of radial derivatives for a liquid layer using the dynamic assumption """

    # Initialize 4 fake radial function solutions
    radial_funcs = tuple([
        np.ones(10, dtype=bulk_array.dtype),
        np.ones(10, dtype=bulk_array.dtype),
        np.ones(10, dtype=bulk_array.dtype),
        np.ones(10, dtype=bulk_array.dtype)
    ])

    # Test for l=2
    solid_derivative = radial_derivatives_liquid_dynamic(radius_array_to_use, radial_funcs, bulk_array, density_array,
                                                         gravity_array, frequency, order_l=2)
    assert type(solid_derivative) == tuple
    assert len(solid_derivative) == 4
    for y_i in range(4):
        assert solid_derivative[y_i].dtype == bulk_array.dtype
        assert solid_derivative[y_i].shape == (10,)

    # Test for l=3
    solid_derivative = radial_derivatives_liquid_dynamic(radius_array_to_use, radial_funcs, bulk_array,
                                                         density_array, gravity_array, frequency, order_l=3)
    assert type(solid_derivative) == tuple
    assert len(solid_derivative) == 4
    for y_i in range(4):
        assert solid_derivative[y_i].dtype == bulk_array.dtype
        assert solid_derivative[y_i].shape == (10,)

    # Test for an array in the frequency variable
    # Initialize 6 fake radial function solutions
    radial_funcs_mtx = tuple([
        np.ones((20, 10), dtype=bulk_array.dtype),
        np.ones((20, 10), dtype=bulk_array.dtype),
        np.ones((20, 10), dtype=bulk_array.dtype),
        np.ones((20, 10), dtype=bulk_array.dtype)
    ])
    freq_domain = np.linspace(-1., 1., 20)
    rad_mtx, freq_mtx = np.meshgrid(radius_array_to_use, freq_domain)
    bulk_mtx, _ = np.meshgrid(bulk_array, freq_domain)
    den_mtx, _ = np.meshgrid(density_array, freq_domain)
    grav_mtx, _ = np.meshgrid(gravity_array, freq_domain)
    solid_guess = radial_derivatives_liquid_dynamic(rad_mtx, radial_funcs_mtx, bulk_mtx,
                                                    den_mtx, grav_mtx, freq_mtx, order_l=3)
    assert type(solid_guess) == tuple
    assert len(solid_guess) == 4
    for sn in range(4):
        assert type(solid_guess[sn]) == np.ndarray
        assert solid_guess[sn].dtype == bulk_mtx.dtype
        # Dynamic solid solution should have 6 y values across the radius domain (10)
        #     and again across the new freq domain (20)
        assert solid_guess[sn].shape == (20, 10)
