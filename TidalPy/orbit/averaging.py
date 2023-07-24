""" Functions used during orbit averaging. """

from typing import List

import numpy as np

from ..utilities.performance import njit, nbList

def orbit_average(orbital_period: float, time_domain: np.ndarray, array_list: List[np.ndarray],
                  radius_N: int, longitude_N: int, colatitude_N: int):

    # Find relevant slice of time domain.
    time_slice = time_domain <= orbital_period
    time_domain_sliced = time_domain[time_slice]

    # The following uses a method similar to `np.trapz`
    N_time = len(time_domain_sliced)
    diff = np.diff(time_domain_sliced)

    arrays_averaged = nbList()
    for array in array_list:
        dtype = array.dtype
        arrays_averaged.append(np.zeros((radius_N, longitude_N, colatitude_N), dtype=dtype))

    # Arrays are looped over radius, longitude, and colatitude
    for t_i in range(N_time):
        # t_i + 1 does not exist at the end of the time domain so once there, break
        if t_i == N_time - 1:
            t_ip1 = -1
            break
        else:
            t_ip1 = t_i + 1
        diff_at_time = diff[t_i] / orbital_period
        diff_over_2 = diff_at_time / 2.
        for r_i in range(radius_N):
            for l_i in range(longitude_N):
                for c_i in range(colatitude_N):
                    # Loop through all provided arrays
                    for array, array_avg in zip(array_list, arrays_averaged):
                        array_avg[r_i, l_i, c_i] += diff_over_2 * (array[r_i, l_i, c_i, t_i] + array[r_i, l_i, c_i, t_ip1])

    return arrays_averaged
