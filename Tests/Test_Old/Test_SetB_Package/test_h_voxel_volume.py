""" Tests for spherical helper function to calculate the volume of spherical voxels. """

import numpy as np

import TidalPy


from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes, calculate_voxel_volumes_npy

planet_radius = 6300.0e3


def test_voxel_volume_numba():
    """ Test the voxel volume calculation using the numba version of the func.
    Check that the results match expectations. """

    radius_array = np.linspace(0., planet_radius, 50)
    longitude_array_deg = np.linspace(0., 360., 20)
    colatitude_array_deg = np.linspace(0., 180., 25)
    longitude_array = np.radians(longitude_array_deg)
    colatitude_array = np.radians(colatitude_array_deg)

    voxel_volumes = calculate_voxel_volumes(radius_array, longitude_array, colatitude_array)

    # Check shape
    assert len(voxel_volumes.shape) == 3
    assert voxel_volumes.shape[0] == len(radius_array)
    assert voxel_volumes.shape[1] == len(longitude_array)
    assert voxel_volumes.shape[2] == len(colatitude_array)

    # See how results compare to expectations
    real_total_volume = (4. / 3.) * np.pi * planet_radius**3
    voxel_total_volume = np.sum(voxel_volumes)
    percent_diff = np.abs(real_total_volume - voxel_total_volume) / real_total_volume

    np.testing.assert_almost_equal(0., percent_diff, decimal=1)


def test_voxel_volume_numpy():
    """ Test the voxel volume calculation using the numpy version of the func.
    Check that the results match expectations. """

    radius_array = np.linspace(0., planet_radius, 20)
    longitude_array_deg = np.linspace(0., 360., 20)
    colatitude_array_deg = np.linspace(0., 180., 20)
    longitude_array = np.radians(longitude_array_deg)
    colatitude_array = np.radians(colatitude_array_deg)

    voxel_volumes = calculate_voxel_volumes_npy(radius_array, longitude_array, colatitude_array)

    # Check shape
    assert len(voxel_volumes.shape) == 3
    assert voxel_volumes.shape[0] == len(radius_array)
    assert voxel_volumes.shape[1] == len(longitude_array)
    assert voxel_volumes.shape[2] == len(colatitude_array)

    # See how results compare to expectations
    real_total_volume = (4. / 3.) * np.pi * planet_radius**3
    voxel_total_volume = np.sum(voxel_volumes)
    percent_diff = np.abs(real_total_volume - voxel_total_volume) / real_total_volume

    np.testing.assert_almost_equal(0., percent_diff, decimal=1)


def test_voxel_volume_numba_higherN():
    """ See if using a higher N will result in better results. """

    radius_array = np.linspace(0., planet_radius, 400)
    longitude_array_deg = np.linspace(0., 360., 100)
    colatitude_array_deg = np.linspace(0., 180., 100)
    longitude_array = np.radians(longitude_array_deg)
    colatitude_array = np.radians(colatitude_array_deg)

    voxel_volumes = calculate_voxel_volumes(radius_array, longitude_array, colatitude_array)

    # Check shape
    assert len(voxel_volumes.shape) == 3
    assert voxel_volumes.shape[0] == len(radius_array)
    assert voxel_volumes.shape[1] == len(longitude_array)
    assert voxel_volumes.shape[2] == len(colatitude_array)

    # See how results compare to expectations
    real_total_volume = (4. / 3.) * np.pi * planet_radius**3
    voxel_total_volume = np.sum(voxel_volumes)
    percent_diff = np.abs(real_total_volume - voxel_total_volume) / real_total_volume

    np.testing.assert_almost_equal(0., percent_diff, decimal=2)
