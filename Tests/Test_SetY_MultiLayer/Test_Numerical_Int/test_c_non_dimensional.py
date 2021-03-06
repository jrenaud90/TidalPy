""" Tests for calculating the multilayer solution using a non-dimensional solver.
"""

import numpy as np

import TidalPy
from TidalPy.constants import G
from TidalPy.tides.multilayer.nondimensional import (non_dimensionalize_physicals,
                                                     re_dimensionalize_physicals, re_dimensionalize_radial_func)
from TidalPy.utilities.types import float_eps

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
bulk_array = 10.e10 * np.ones(10, dtype=np.complex128)
radius_array_to_use = radius_array[1:]
frequency = 2. * np.pi / (86400. * 1.)


def test_non_dimensionalize_physicals():
    # Create non-dimensionalized properties.
    R_nd = 1.123e6
    rho_nd = 1030.
    result = non_dimensionalize_physicals(
        radius_array[1:], gravity_array, density_array, shear_array, bulk_array,
        frequency, R_nd, rho_nd
        )
    # Perform type checks
    assert type(result) == tuple
    radius_prime, gravity_prime, density_prime, shear_modulus_prime, bulk_modulus_prime, frequency_prime, \
    newton_g_prime = result

    for arr in [radius_prime, gravity_prime, density_prime]:
        assert type(arr) == np.ndarray
        assert arr.dtype in [float, np.float, np.float64]
        assert arr.shape == (10,)

    for arr in [shear_modulus_prime, bulk_modulus_prime]:
        assert type(arr) == np.ndarray
        assert arr.dtype in [complex, np.complex, np.complex128]
        assert arr.shape == (10,)

    assert type(frequency_prime) in [float, np.float, np.float64]
    assert type(newton_g_prime) in [float, np.float, np.float64]

    # Try to redimensionalize them
    result = \
        re_dimensionalize_physicals(
            radius_prime, gravity_prime, density_prime, shear_modulus_prime, bulk_modulus_prime,
            frequency_prime, R_nd, rho_nd
            )
    # Perform type checks
    assert type(result) == tuple
    radius_dbl_prime, gravity_dbl_prime, density_dbl_prime, shear_modulus_dbl_prime, bulk_modulus_dbl_prime, \
    frequency_dbl_prime = result

    # Perform type checks
    for arr in [radius_dbl_prime, gravity_dbl_prime, density_dbl_prime]:
        assert type(arr) == np.ndarray
        assert arr.dtype in [float, np.float, np.float64]
        assert arr.shape == (10,)

    for arr in [shear_modulus_dbl_prime, bulk_modulus_dbl_prime]:
        assert type(arr) == np.ndarray
        assert arr.dtype in [complex, np.complex, np.complex128]
        assert arr.shape == (10,)

    assert type(frequency_dbl_prime) in [float, np.float, np.float64]

    # Compare results. If both functions worked correctly the double prime parameters should be equal to the initial input.
    np.testing.assert_allclose(radius_array[1:], radius_dbl_prime)
    np.testing.assert_allclose(gravity_array, gravity_dbl_prime)
    np.testing.assert_allclose(density_array, density_dbl_prime)
    np.testing.assert_allclose(shear_array, shear_modulus_dbl_prime)
    np.testing.assert_allclose(bulk_array, bulk_modulus_dbl_prime)

    assert np.abs(frequency_dbl_prime - frequency) < float_eps


def test_non_dimensionalize_tidaly():
    # TODO TESTS - this is not actually testing that it is preforming the conversion correctly. Just a unit and process check.

    # Make a random tidal y array
    tidal_y_prime = np.stack(
        (
            np.linspace(0.1, 10.0, 10),
            np.linspace(0.8, 10.0, 10),
            np.linspace(11.0, 10.0, 10),
            np.linspace(0.1, 10.0, 10),
            np.linspace(83.2, 10.0, 10),
            np.linspace(0.1, 10.0, 10),
            )
        )

    # Attempt to re-dimensionalize it
    R_nd = 1.123e6
    rho_nd = 1030.
    result = re_dimensionalize_radial_func(tidal_y_prime, R_nd, rho_nd)

    assert type(result) == np.ndarray
    assert result.shape == (6, 10)
    assert result.dtype in [float, np.float, np.float64]
