from typing import Union

import numpy as np

from TidalPy.exceptions import BadArrayShape
from ..performance import njit
from ..types import FloatArray


@njit
def find_nearest(array: np.ndarray, value: Union[int, float]):
    """ Returns the index of the value closest to a provided value in an array

    :param array: <ndarray> Array to search
    :param value: <float or int> value to search for
    :return:      <int> index of nearest value
    """

    array = np.asarray(array)
    return (np.abs(array - value)).argmin()


def value_cleanup(value):
    if type(value) != np.ndarray:
        value = np.asarray([value])

    return value


def match_array(array_to_be_matched: FloatArray, *reference_arrays):
    reference_array = None
    # Pick the reference array that is not None and preferentially one that is not (1,).
    for ref_array in reference_arrays:
        if ref_array is None:
            continue
        if type(ref_array) != np.ndarray:
            continue
        reference_array = ref_array
        if ref_array.shape != (1,):
            break

    # Covert input to array
    if type(array_to_be_matched) != np.ndarray:
        array_to_be_matched = np.asarray([array_to_be_matched])

    if reference_array is not None:

        # Check to see if length is 1 and if so set it equal to reference array length
        if array_to_be_matched.shape == (1,) and reference_array.shape != (1,):
            array_to_be_matched = array_to_be_matched * np.ones_like(reference_array)
        else:
            if reference_array.shape != (1,):
                if array_to_be_matched.shape != reference_array.shape:
                    raise BadArrayShape

    return array_to_be_matched


def neg_array_for_log_plot(array_with_negatives: np.ndarray):
    """ Converts one numpy array into two where both new arrays only have positive values. Useful for log-plotting.

    Parameters
    ----------
    array_with_negatives : np.ndarray
        Numpy array that may have negatives or zeros

    Returns
    -------
    array_positive : np.ndarray
        Numpy array with original array's positives
    array_negative : np.ndarray
        Numpy array with original array's negatives set to positive
    """

    assert type(array_with_negatives) == np.ndarray

    array_positive = array_with_negatives.copy()
    array_negative = array_with_negatives.copy()

    array_positive[array_positive <= 0.] = np.nan
    array_negative[array_negative >= 0.] = np.nan
    array_negative[array_negative < 0.] *= -1.

    return array_positive, array_negative
