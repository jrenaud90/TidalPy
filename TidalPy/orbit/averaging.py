""" Functions used during orbit averaging. """

from typing import List

import numpy as np

from ..utilities.performance import njit, nbList

@njit
def orbit_average(orbital_period: float, time_domain: np.ndarray, array: np.ndarray, scale_by_period: bool = True):
    """ Calculates the orbit averaged values for the provided array. 
    This function assumes that the provided array is only a function of time (one dimensional).

    Parameters
    ----------
    orbital_period : float
        Orbital period [must be the same units as `time_domain`]
    time_domain : np.ndarray
        Time domain array (1D array) [must be the same units as `orbital_period`]
    array : np.ndarray
        Array that will be orbit averaged.
        This function assumes that the provided array is only a function of time (shape should equal shape of time_domain).

    Returns
    -------
    array_averaged : np.ndarray
        The orbit averaged array.
        This array will be one less dimension than the input array (since the time domain information has collapsed).

    Raises
    ------
    ValueError
        No points in the time domain found within orbital period.
    ValueError
        Very few points in the time domain found within orbital period.

    """

    # Limit the time domain to only include points that are within one orbital period.
    time_slice = time_domain <= orbital_period
    time_domain_sliced = time_domain[time_slice]
    n_time = len(time_domain_sliced)
    if n_time == 0:
        raise ValueError('No points in the time domain found within orbital period.')
    elif n_time <= 3:
        raise ValueError('Very few points in the time domain found within orbital period.')

    # The following uses a method similar to `np.trapz`
    diff = np.diff(time_domain_sliced)
    if scale_by_period:
        diff = diff / orbital_period
    
    # Array is looped over time, radius, longitude, and colatitude
    array_averaged = 0.
    for t_i in range(n_time):
        # t_i + 1 does not exist at the end of the time domain so once there, break
        if t_i == n_time - 1:
            t_ip1 = -1
            break
        else:
            t_ip1 = t_i + 1
        diff_at_time = diff[t_i]
        diff_over_2 = diff_at_time / 2.
        # Add this time step's contribution
        array_averaged += diff_over_2 * (array[t_i] + array[t_ip1])

    return array_averaged


@njit
def orbit_average_3d(orbital_period: float, time_domain: np.ndarray, array: np.ndarray,
                     n_radius: int, n_longitude: int, n_colatitude: int, scale_by_period: bool = True):
    """ Calculates the orbit averaged values for the provided array. 
    This function assumes that arrays include 3D spatial information: radius, longitude, and colatitude.

    Parameters
    ----------
    orbital_period : float
        Orbital period [must be the same units as `time_domain`]
    time_domain : np.ndarray
        Time domain array (1D array) [must be the same units as `orbital_period`]
    array : np.ndarray
        Array that will be orbit averaged.
        This function assumes this array is structured with the shape (radius, longitude, colatitude, time).
    n_radius : int
        Number of radius points.
    n_longitude : int
        Number of longitude points.
    n_colatitude : int
        Number of colatitude points.

    Returns
    -------
    array_averaged : np.ndarray
        The orbit averaged array.
        This array will be one less dimension than the input array (since the time domain information has collapsed).
    
    Raises
    ------
    ValueError
        No points in the time domain found within orbital period.
    ValueError
        Very few points in the time domain found within orbital period.

    """

    # Limit the time domain to only include points that are within one orbital period.
    time_slice = time_domain <= orbital_period
    time_domain_sliced = time_domain[time_slice]
    n_time = len(time_domain_sliced)
    if n_time == 0:
        raise ValueError('No points in the time domain found within orbital period.')
    elif n_time <= 3:
        raise ValueError('Very few points in the time domain found within orbital period.')
    
    # The following uses a method similar to `np.trapz`
    n_time = len(time_domain_sliced)
    diff = np.diff(time_domain_sliced)
    if scale_by_period:
        diff = diff / orbital_period
    
    # Array is looped over time, radius, longitude, and colatitude
    array_averaged = np.zeros((n_radius, n_longitude, n_colatitude), dtype=array.dtype)

    for t_i in range(n_time):
        # t_i + 1 does not exist at the end of the time domain so once there, break
        if t_i == n_time - 1:
            t_ip1 = -1
            break
        else:
            t_ip1 = t_i + 1
        diff_at_time = diff[t_i]
        diff_over_2 = diff_at_time / 2.
        for r_i in range(n_radius):
            for l_i in range(n_longitude):
                for c_i in range(n_colatitude):
                    # Add this time step's contribution
                    array_averaged[r_i, l_i, c_i] += diff_over_2 * (array[r_i, l_i, c_i, t_i] + array[r_i, l_i, c_i, t_ip1])

    return array_averaged


@njit
def orbit_average_3d_multiarray(orbital_period: float, time_domain: np.ndarray, array_list: List[np.ndarray],
                                n_radius: int, n_longitude: int, n_colatitude: int, scale_by_period: bool = True):
    """ Calculates the orbit averaged values for the provided arrays. 
    This function assumes that arrays include 3D spatial information: radius, longitude, and colatitude.

    Parameters
    ----------
    orbital_period : float
        Orbital period [must be the same units as `time_domain`]
    time_domain : np.ndarray
        Time domain array (1D array) [must be the same units as `orbital_period`]
    array_list : List[np.ndarray]
        List of arrays that will be orbit averaged. All arrays must have the same dimensions and dtypes.
        This function assumes these arrays are structured with the shape (radius, longitude, colatitude, time).
    radius_N : int
        Number of radius points.
    longitude_N : int
        Number of longitude points.
    colatitude_N : int
        Number of colatitude points.
    scale_by_period : bool, optional
      If True, then the result will be scaled by 1/orbital_period, by default True

    Returns
    -------
    arrays_averaged : List[np.ndarray]
        The orbit averaged arrays stored in the same order as `array_list`.
        Each array will be one less dimension than the input array (since the time domain information has collapsed).

    Raises
    ------
    ValueError
        No points in the time domain found within orbital period.
    ValueError
        Very few points in the time domain found within orbital period.
    TypeError
        Arrays must have the same dtype.
    TypeError
        Arrays must have the same shape.
    """

    # Limit the time domain to only include points that are within one orbital period.
    time_slice = time_domain <= orbital_period
    time_domain_sliced = time_domain[time_slice]
    n_time = len(time_domain_sliced)
    if n_time == 0:
        raise ValueError('No points in the time domain found within orbital period.')
    elif n_time <= 3:
        raise ValueError('Very few points in the time domain found within orbital period.')
    
    # The following uses a method similar to `np.trapz`
    n_time = len(time_domain_sliced)
    diff = np.diff(time_domain_sliced)
    if scale_by_period:
        diff = diff / orbital_period

    arrays_averaged = nbList()
    dtype_old = array_list[0].dtype
    shape_old = array_list[0].shape
    for array in array_list:
        dtype = array.dtype
        shape = array.shape
        if not dtype is dtype_old:
            raise TypeError('Arrays must have the same dtype.')
        if not shape is shape_old:
            raise TypeError('Arrays must have the same shape.')
        arrays_averaged.append(np.zeros((n_radius, n_longitude, n_colatitude), dtype=dtype))
        
    # Arrays are looped over time, radius, longitude, and colatitude
    for t_i in range(n_time):
        # t_i + 1 does not exist at the end of the time domain so once there, break
        if t_i == n_time - 1:
            t_ip1 = -1
            break
        else:
            t_ip1 = t_i + 1
        diff_at_time = diff[t_i]
        diff_over_2 = diff_at_time / 2.
        for r_i in range(n_radius):
            for l_i in range(n_longitude):
                for c_i in range(n_colatitude):
                    # Loop through all provided arrays and add this time step's contribution
                    for array, array_avg in zip(array_list, arrays_averaged):
                        array_avg[r_i, l_i, c_i] += diff_over_2 * (array[r_i, l_i, c_i, t_i] + array[r_i, l_i, c_i, t_ip1])

    return arrays_averaged

@njit
def orbit_average_4d_multiarray(orbital_period: float, time_domain: np.ndarray, array_list: List[np.ndarray],
                                n_outer: int, n_radius: int, n_longitude: int, n_colatitude: int, scale_by_period: bool = True):
    """ Calculates the orbit averaged values for the provided arrays. 
    This function assumes that arrays include 3D spatial information: radius, longitude, and colatitude.

    Parameters
    ----------
    orbital_period : float
        Orbital period [must be the same units as `time_domain`]
    time_domain : np.ndarray
        Time domain array (1D array) [must be the same units as `orbital_period`]
    array_list : List[np.ndarray]
        List of arrays that will be orbit averaged. All arrays must have the same dimensions and dtypes.
        This function assumes these arrays are structured with the shape (radius, longitude, colatitude, time).
    radius_N : int
        Number of radius points.
    longitude_N : int
        Number of longitude points.
    colatitude_N : int
        Number of colatitude points.
    scale_by_period : bool, optional
      If True, then the result will be scaled by 1/orbital_period, by default True

    Returns
    -------
    arrays_averaged : List[np.ndarray]
        The orbit averaged arrays stored in the same order as `array_list`.
        Each array will be one less dimension than the input array (since the time domain information has collapsed).

    Raises
    ------
    ValueError
        No points in the time domain found within orbital period.
    ValueError
        Very few points in the time domain found within orbital period.
    TypeError
        Arrays must have the same dtype.
    TypeError
        Arrays must have the same shape.
    """

    # Limit the time domain to only include points that are within one orbital period.
    time_slice = time_domain <= orbital_period
    time_domain_sliced = time_domain[time_slice]
    n_time = len(time_domain_sliced)
    if n_time == 0:
        raise ValueError('No points in the time domain found within orbital period.')
    elif n_time <= 3:
        raise ValueError('Very few points in the time domain found within orbital period.')
    
    # The following uses a method similar to `np.trapz`
    n_time = len(time_domain_sliced)
    diff = np.diff(time_domain_sliced)
    if scale_by_period:
        diff = diff / orbital_period

    arrays_averaged = nbList()
    dtype_old = array_list[0].dtype
    shape_old = array_list[0].shape
    for array in array_list:
        dtype = array.dtype
        shape = array.shape
        if not dtype is dtype_old:
            raise TypeError('Arrays must have the same dtype.')
        if not shape is shape_old:
            raise TypeError('Arrays must have the same shape.')
        arrays_averaged.append(np.zeros((n_outer, n_radius, n_longitude, n_colatitude), dtype=dtype))
        
    # Arrays are looped over time, radius, longitude, and colatitude
    for t_i in range(n_time):
        # t_i + 1 does not exist at the end of the time domain so once there, break
        if t_i == n_time - 1:
            t_ip1 = -1
            break
        else:
            t_ip1 = t_i + 1
        diff_at_time = diff[t_i]
        diff_over_2 = diff_at_time / 2.
        for o_i in range(n_outer):
            for r_i in range(n_radius):
                for l_i in range(n_longitude):
                    for c_i in range(n_colatitude):
                        # Loop through all provided arrays and add this time step's contribution
                        for array, array_avg in zip(array_list, arrays_averaged):
                            array_avg[o_i, r_i, l_i, c_i] += diff_over_2 * (array[o_i, r_i, l_i, c_i, t_i] + array[o_i, r_i, l_i, c_i, t_ip1])

    return arrays_averaged