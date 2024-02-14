""" Tests for spherical helper function to calculate the volume of spherical voxels. """

import numpy as np

import TidalPy


from TidalPy.constants import G
from TidalPy.utilities.spherical_helper.mass import calculate_mass_gravity_arrays

radius = np.linspace(0.01e3, 6378.1e3, 50)
layer_0 = radius < 3483.e3
layer_1 = radius >= 3483.e3
l1_n = len(radius[layer_1])
density = np.zeros_like(radius)
density[layer_0] = 10000.
density[layer_1] = 4500.


def test_calculate_mass_gravity_arrays():
    """ Test the calculation of volume, mass, and gravity from radius and density arrays. """

    volumes, masses, gravities = calculate_mass_gravity_arrays(radius, density)

    # Check array sizes
    assert len(volumes) == len(radius)
    assert len(masses) == len(radius)
    assert len(gravities) == len(radius)

    # Check values
    real_volume = (4. / 3.) * np.pi * radius[-1]**3
    np.testing.assert_almost_equal(np.abs(np.sum(volumes) - real_volume)/real_volume, 0.)
    mass = np.sum(masses)
    real_gravity_surf = G * mass / radius[-1]**2
    np.testing.assert_almost_equal(np.abs(gravities[-1] - real_gravity_surf)/real_gravity_surf, 0.)