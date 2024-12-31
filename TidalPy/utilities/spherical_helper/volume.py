""" Functions to help calculate the volume of voxels in a spherical geometry
"""
import numpy as np

from TidalPy.utilities.performance import njit


def calculate_voxel_volumes_npy(radius_array: np.ndarray, longitude_array: np.ndarray, colatitude_array: np.ndarray):
    """ Calculate the volume of all voxels within a sphere assuming spherical geometry.

    This function uses numpy functions that are not currently supported by numba.

    Longitude and Colatitude must be provided in radians.

    Parameters
    ----------
    radius_array : np.ndarray
        Radius array, does not have to be equally spaced [m]
        Shape: N_r
    longitude_array : np.ndarray
        Longitude array, does not have to be equally spaced [rad]
        Shape: N_long
    colatitude_array : np.ndarray
        Colatitude array, does not have to be equally spaced [rad]
        Shape: N_colat

    Returns
    -------
    voxel_volumes : np.ndarray
        Volume for each voxel assuming spherical geometry [m3]
        Shape: (N_r, N_long, N_colat)

    """

    # Find the separation between members of each array.
    # Note we do not assume members are equally spaced (e.g., dr near the base may not equal dr near the top)
    #    This is to allow for different sized (numerical size) layers within a planet.
    dr_array = np.diff(radius_array)
    dlong_array = np.diff(longitude_array)
    dcolat_array = np.diff(colatitude_array)

    # The difference arrays are not as long as the input arrays. To counteract this we will assume the starting
    #    difference is the same as the first calculated difference. This should be fine unless there is a layer
    #    structure with a very small N (likewise if longitude or latitude are small N).
    dr_array = np.concatenate(([dr_array[0]], dr_array))
    dlong_array = np.concatenate(([dlong_array[0]], dlong_array))
    dcolat_array = np.concatenate(([dcolat_array[0]], dcolat_array))

    # Need to make a meshgrid to increase the dimensions of the values.
    radius_matrix, long_matrix, colat_matrix = \
        np.meshgrid(radius_array, longitude_array, colatitude_array, indexing='ij')

    # Same for the differences
    dr_matrix, dlong_matrix, dcolat_matrix = \
        np.meshgrid(dr_array, dlong_array, dcolat_array, indexing='ij')

    # Solve for the voxel volume
    voxel_volumes = radius_matrix**2 * np.sin(colat_matrix) * dr_matrix * dlong_matrix * dcolat_matrix

    return voxel_volumes


@njit(cacheable=True)
def calculate_voxel_volumes(radius_array: np.ndarray, longitude_array: np.ndarray, colatitude_array: np.ndarray):
    """ Calculate the volume of all voxels within a sphere assuming spherical geometry.

    This function uses loops which are supported by numba. Currently, this function is about 160% faster than the
        current numpy implementation above. Setting parallel to False will lead to a drop in performance of about 120%.

    Longitude and Colatitude must be provided in radians.

    Parameters
    ----------
    radius_array : np.ndarray
        Radius array, does not have to be equally spaced [m]
        Shape: N_r
    longitude_array : np.ndarray
        Longitude array, does not have to be equally spaced [rad]
        Shape: N_long
    colatitude_array : np.ndarray
        Colatitude array, does not have to be equally spaced [rad]
        Shape: N_colat

    Returns
    -------
    voxel_volumes : np.ndarray
        Volume for each voxel assuming spherical geometry [m3]
        Shape: (N_r, N_long, N_colat)

    """

    # Find the separation between members of each array.
    # Note we do not assume members are equally spaced (e.g., dr near the base may not equal dr near the top)
    #    This is to allow for different sized (numerical size) layers within a planet.
    dr_array = np.diff(radius_array)
    dlong_array = np.diff(longitude_array)
    dcolat_array = np.diff(colatitude_array)

    # The difference arrays are not as long as the input arrays. To counteract this we will assume the starting
    #    difference is the same as the first calculated difference. This should be fine unless there is a layer
    #    structure with a very small N (likewise if longitude or latitude are small N).
    dr_array = np.concatenate((np.asarray([dr_array[0]], dtype=np.float64), dr_array))
    dlong_array = np.concatenate((np.asarray([dlong_array[0]], dtype=np.float64), dlong_array))
    dcolat_array = np.concatenate((np.asarray([dcolat_array[0]], dtype=np.float64), dcolat_array))

    # Get size of arrays
    n_r = len(radius_array)
    n_long = len(longitude_array)
    n_colat = len(colatitude_array)

    # Initialize empty array for voxel volumes
    voxel_volumes = np.empty((n_r, n_long, n_colat), dtype=np.float64)

    for i_r in range(n_r):
        r = radius_array[i_r]
        dr = dr_array[i_r]
        for i_long in range(n_long):
            dlong = dlong_array[i_long]
            for i_colat in range(n_colat):
                colat = colatitude_array[i_colat]
                dcolat = dcolat_array[i_colat]

                # Calculate this voxel's volume.
                voxel_volumes[i_r, i_long, i_colat] = r**2 * np.sin(colat) * dr * dlong * dcolat

    return voxel_volumes
