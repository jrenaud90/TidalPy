""" Tests for calculating the initial guess at the bottom of a liquid or solid layer for numerical integration models
"""

import numpy as np

import TidalPy
from TidalPy.constants import G
from TidalPy.tides.multilayer.numerical_int import (liquid_dynamic_guess, liquid_static_guess, solid_dynamic_guess,
                                                    solid_static_guess)

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


def test_static_solid():
    """ Tests the initial guess at the bottom of a static solid layer """

    # Test for order l = 2
    solid_guess = solid_static_guess(radius_array_to_use, shear_array, bulk_array, density_array, order_l=2)
    assert type(solid_guess) == tuple
    assert len(solid_guess) == 3
    for sn in range(3):
        assert type(solid_guess[sn]) == np.ndarray
        assert solid_guess[sn].dtype == shear_array.dtype
        # Static solid solution should have 6 y values across the radius domain (10).
        assert solid_guess[sn].shape == (6, 10)

    # Test for order l = 3
    solid_guess = solid_static_guess(radius_array_to_use, shear_array, bulk_array, density_array, order_l=3)
    assert type(solid_guess) == tuple
    assert len(solid_guess) == 3
    for sn in range(3):
        assert type(solid_guess[sn]) == np.ndarray
        assert solid_guess[sn].dtype == shear_array.dtype
        # Static solid solution should have 6 y values across the radius domain (10).
        assert solid_guess[sn].shape == (6, 10)


def test_dynamic_solid():
    """ Tests the initial guess at the bottom of a dynamic solid layer """

    # Test for order l = 2
    solid_guess = solid_dynamic_guess(radius_array_to_use, shear_array, bulk_array, density_array, frequency, order_l=2)
    assert type(solid_guess) == tuple
    assert len(solid_guess) == 3
    for sn in range(3):
        assert type(solid_guess[sn]) == np.ndarray
        assert solid_guess[sn].dtype == shear_array.dtype
        # Dynamic solid solution should have 6 y values across the radius domain (10).
        assert solid_guess[sn].shape == (6, 10)

    # Test for order l = 3
    solid_guess = solid_dynamic_guess(radius_array_to_use, shear_array, bulk_array, density_array, frequency, order_l=3)
    assert type(solid_guess) == tuple
    assert len(solid_guess) == 3
    for sn in range(3):
        assert type(solid_guess[sn]) == np.ndarray
        assert solid_guess[sn].dtype == shear_array.dtype
        # Dynamic solid solution should have 6 y values across the radius domain (10).
        assert solid_guess[sn].shape == (6, 10)

    # Test for an array in the frequency variable
    freq_domain = np.linspace(-1., 1., 20)
    rad_mtx, freq_mtx = np.meshgrid(radius_array_to_use, freq_domain)
    shear_mtx, _ = np.meshgrid(shear_array, freq_domain)
    bulk_mtx, _ = np.meshgrid(bulk_array, freq_domain)
    den_mtx, _ = np.meshgrid(density_array, freq_domain)
    solid_guess = solid_dynamic_guess(rad_mtx, shear_mtx, bulk_mtx, den_mtx, freq_mtx, order_l=2)
    assert type(solid_guess) == tuple
    assert len(solid_guess) == 3
    for sn in range(3):
        assert type(solid_guess[sn]) == np.ndarray
        assert solid_guess[sn].dtype == shear_mtx.dtype
        # Dynamic solid solution should have 6 y values across the radius domain (10)
        #     and again across the new freq domain (20)
        assert solid_guess[sn].shape == (6, 20, 10)


def test_static_liquid():
    """ Tests the initial guess at the bottom of a static liquid layer """

    # Test for order l = 2
    # For the static liquid guess there is only one solution, therefore the output is NOT a tuple.
    liquid_guess = liquid_static_guess(radius_array_to_use, order_l=2)
    assert type(liquid_guess) == np.ndarray
    assert type(liquid_guess) == np.ndarray
    # Since there is no shear (or bulk) dependence for static liquid tides, then the results will be real not complex.
    #    The above is true, but in 0.3.0a5 a forced complex "asarray" was added.
    assert liquid_guess.dtype == 'complex128'
    # Static liquid solution should have 2 y values across the radius domain (10).
    assert liquid_guess.shape == (2, 10)

    # Test for order l = 3
    liquid_guess = liquid_static_guess(radius_array_to_use, order_l=3)
    assert type(liquid_guess) == np.ndarray
    assert type(liquid_guess) == np.ndarray
    # Since there is no shear (or bulk) dependence for static liquid tides, then the results will be real not complex.
    #    The above is true, but in 0.3.0a5 a forced complex "asarray" was added.
    assert liquid_guess.dtype == 'complex128'
    # Static liquid solution should have 2 y values across the radius domain (10).
    assert liquid_guess.shape == (2, 10)


def test_dynamic_liquid():
    """ Tests the initial guess at the bottom of a dynamic liquid layer """

    # Test for order l = 2
    liquid_guess = liquid_dynamic_guess(radius_array_to_use, bulk_array, density_array, frequency, order_l=2)
    assert type(liquid_guess) == tuple
    assert len(liquid_guess) == 2
    for sn in range(2):
        assert type(liquid_guess[sn]) == np.ndarray
        # Since there is no shear dependence for dynamic liquid tides, then the results will be real not complex unless
        #    there is bulk dependence
        #    The above is true, but in 0.3.0a5 a forced complex "asarray" was added.
        assert liquid_guess[sn].dtype == 'complex128'
        # Dynamic liquid solution should have 4 y values across the radius domain (10).
        assert liquid_guess[sn].shape == (4, 10)

    # Test for order l = 3
    liquid_guess = liquid_dynamic_guess(radius_array_to_use, bulk_array, density_array, frequency, order_l=3)
    assert type(liquid_guess) == tuple
    assert len(liquid_guess) == 2
    for sn in range(2):
        assert type(liquid_guess[sn]) == np.ndarray
        # Since there is no shear dependence for dynamic liquid tides, then the results will be real not complex unless
        #    there is bulk dependence
        #    The above is true, but in 0.3.0a5 a forced complex "asarray" was added.
        assert liquid_guess[sn].dtype == 'complex128'
        # Dynamic liquid solution should have 4 y values across the radius domain (10).
        assert liquid_guess[sn].shape == (4, 10)

    # Test for an array in the frequency variable
    freq_domain = np.linspace(-1., 1., 20)
    rad_mtx, freq_mtx = np.meshgrid(radius_array_to_use, freq_domain)
    bulk_mtx, _ = np.meshgrid(bulk_array, freq_domain)
    den_mtx, _ = np.meshgrid(density_array, freq_domain)
    liquid_guess = liquid_dynamic_guess(rad_mtx, bulk_mtx, den_mtx, freq_mtx, order_l=2)
    assert type(liquid_guess) == tuple
    assert len(liquid_guess) == 2
    for sn in range(2):
        assert type(liquid_guess[sn]) == np.ndarray
        # Since there is no shear dependence for static liquid tides, then the results will be real not complex unless
        #    there is bulk dependence
        #    The above is true, but in 0.3.0a5 a forced complex "asarray" was added.
        assert liquid_guess[sn].dtype == 'complex128'
        # Dynamic liquid solution should have 4 y values across the radius domain (10)
        #     and again across the new freq domain (20)
        assert liquid_guess[sn].shape == (4, 20, 10)
