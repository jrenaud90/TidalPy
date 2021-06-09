from typing import Dict

import numpy as np

from ..performance.numba import njit
from ..types import FloatArray, NumericalType
from ...exceptions import BadArrayShape


def normalize_dict(dict_of_values: Dict[str, np.ndarray], pass_negatives: bool = False,
                   new_max: float = 1.0, new_min: float = 0.0):
    """ Normalizes values provided in a name separated dictionary to the specified range.

    Parameters
    ----------
    dict_of_values : Dict[str, float]
        Dictionary of reference keys pointing to the to-be-normalized values.
    pass_negatives : bool = False
        If true then any values that are negative will be excluded from the normalization.
    new_max : float = 1.0
        The upper-most of the post-normalized values
    new_min : float = 0.0
        The lower-most of the post-normalized values

    Returns
    -------
    dict_of_normalized_values : Dict[str, float]
        Dictionary of reference keys pointing to the post-normalized values.

    """

    if pass_negatives:
        max_ = 1.0e-10
    else:
        max_ = -1.0e100
    min_ = 1.0e100

    for ref_name, value in dict_of_values.items():
        if np.max(value) > max_:
            max_ = np.max(value)
        if min_ > np.min(value):
            if pass_negatives and np.min(value) < 0.:
                # Negative values may indicate an integration problem and the user may want to exclude them from the
                # normalization.
                if min_ > np.min(value[value >= 0]):
                    min_ = np.min(value[value >= 0])
            else:
                min_ = np.min(value)

    if min_ == max_:
        min_ = 0.9 * max_

    new_dict = {ref: np.zeros_like(value) for ref, value in dict_of_values.items()}
    slope = (new_max - new_min) / (max_ - min_)
    intercept = new_max - slope * max_
    for ref_name, value in dict_of_values.items():
        if pass_negatives:
            new_dict[ref_name][value >= 0.] = slope * value[value >= 0.] + intercept
            new_dict[ref_name][value < 0.] = value[value < 0.]
        else:
            new_dict[ref_name] = slope * value + intercept
    return new_dict


@njit()
def find_nearest(array: np.ndarray, value: NumericalType):
    """ Returns the index of the value closest to a provided value in a numpy array.

    Parameters
    ----------
    array : np.ndarray
        Numpy array to search in
    value : NumericalType
        Value to search for within array

    Returns
    -------
    index: int
        Index of nearest value
    """

    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index


def value_np_cleanup(value):
    if type(value) in [float, int]:
        value = np.asarray(value)

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
