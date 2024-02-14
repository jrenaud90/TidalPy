import pytest

import numpy as np

import TidalPy
TidalPy.test_mode = True

from scipy.constants import G as G_

from TidalPy.utilities.dimensions.nondimensional import non_dimensionalize_physicals, redimensionalize_physicals


def test_non_dimensionalize_physicals():

    frequency     = 1.0e-3
    mean_radius   = 1.0e6
    bulk_density  = 5500.
    radius_array  = np.linspace(0., mean_radius, 10)
    density_array = np.linspace(7000, 3500, 10)
    bulk_array    = 100.0e9 * np.ones_like(radius_array)
    gravity_array = np.linspace(0.1, 7.0, 10)
    shear_array   = (50.0e9 + 20.0e3j) * np.ones(radius_array.size, dtype=np.complex128)

    # Make copies of the arrays so we can keep the original values for comparison.
    radius_array_out  = radius_array.copy()
    density_array_out = density_array.copy()
    bulk_array_out    = bulk_array.copy()
    gravity_array_out = gravity_array.copy()
    shear_array_out   = shear_array.copy()    

    # non-dimensionalize the arrays and get some non-dim constants.
    output_freq_nondim, G_nondim = non_dimensionalize_physicals(
        frequency, mean_radius, bulk_density,
        radius_array_out, density_array_out, gravity_array_out, bulk_array_out, shear_array_out)
    
    # Build conversions for comparison
    second2_conversion = 1. / (np.pi * G_ * bulk_density)
    second_conversion  = np.sqrt(second2_conversion)
    length_conversion  = mean_radius
    density_conversion = bulk_density
    mass_conversion    = bulk_density * mean_radius**3
    pascal_conversion  = mass_conversion / (length_conversion * second2_conversion)

    # Check conversion is correct
    assert np.allclose(radius_array / length_conversion, radius_array_out)
    assert np.allclose(density_array / density_conversion, density_array_out)
    assert np.allclose(gravity_array / (length_conversion / second2_conversion), gravity_array_out)
    assert np.allclose(bulk_array / pascal_conversion, bulk_array_out)
    assert np.allclose(shear_array / pascal_conversion, shear_array_out)

    assert np.isclose(G_nondim, G_ / (length_conversion**3 / (mass_conversion * second2_conversion)))
    assert np.isclose(output_freq_nondim, frequency / (1. / second_conversion))

    # Check that the redimensionalize function works too.
    # Use the output from the non-dim as inputs to the redim.
    output_freq_redim, G_redim = redimensionalize_physicals(
        output_freq_nondim, mean_radius, bulk_density,
        radius_array_out, density_array_out, gravity_array_out, bulk_array_out, shear_array_out)

    # These should all match the original again.
    assert np.allclose(radius_array, radius_array_out)
    assert np.allclose(density_array, density_array_out)
    assert np.allclose(gravity_array, gravity_array_out)
    assert np.allclose(bulk_array, bulk_array_out)
    assert np.allclose(shear_array, shear_array_out)
    assert np.isclose(G_redim, G_)
    assert np.isclose(output_freq_redim, frequency)
