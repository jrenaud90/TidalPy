import numpy as np
from numba import njit
from typing import Union

@njit
def find_nearest(array: np.ndarray, value: Union[int, float]):
    """ Returns the index of the value closest to a provided value in an array

    :param array: <ndarray> Array to search
    :param value: <float or int> value to search for
    :return:      <int> index of nearest value
    """

    array = np.asarray(array)
    return (np.abs(array - value)).argmin()
