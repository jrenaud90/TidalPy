from typing import Dict, TYPE_CHECKING

import numpy as np

from TidalPy.utilities.types import float_eps
from TidalPy.utilities.performance.numba import njit

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray, NumericalType


def normalize_dict(
    dict_of_values: Dict[str, 'FloatArray'], pass_negatives: bool = False,
    new_max: float = 1.0, new_min: float = 0.0
    ):
    """ Normalizes values provided in a name separated dictionary to the specified range.

    Parameters
    ----------
    dict_of_values : Dict[str, FloatArray]
        Dictionary of reference keys pointing to the to-be-normalized values.
    pass_negatives : bool = False
        If true then any values that are negative will be excluded from the normalization.
    new_max : float = 1.0
        The upper-most of the post-normalized values
    new_min : float = 0.0
        The lower-most of the post-normalized values

    Returns
    -------
    dict_of_normalized_values : Dict[str, FloatArray]
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
            if pass_negatives and np.min(value) < -float_eps:
                # Negative values may indicate an integration problem and the user may want to exclude them from the
                # normalization.
                if min_ > np.min(value[value >= float_eps]):
                    min_ = np.min(value[value >= float_eps])
            else:
                min_ = np.min(value)

    if min_ == max_:
        min_ = 0.9 * max_

    new_dict = {ref: np.zeros_like(value) for ref, value in dict_of_values.items()}
    slope = (new_max - new_min) / (max_ - min_)
    intercept = new_max - slope * max_
    for ref_name, value in dict_of_values.items():
        if pass_negatives:
            new_dict[ref_name][value >= float_eps] = slope * value[value >= float_eps] + intercept
            new_dict[ref_name][value < -float_eps] = value[value < -float_eps]
        else:
            new_dict[ref_name] = slope * value + intercept
    return new_dict


@njit(cacheable=True)
def find_nearest(array: np.ndarray, value: 'NumericalType'):
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

@njit(cacheable=True)
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

    # TODO: Numba currently does not support fancy indexing with 2+D arrays.
    #    for now we will flatten and reshape the array.
    org_shape = array_with_negatives.shape

    array_positive = array_with_negatives.copy().flatten()
    array_negative = array_with_negatives.copy().flatten()

    array_positive[array_positive <= 0.] = np.nan
    array_negative[array_negative >= 0.] = np.nan
    array_negative = np.abs(array_negative)

    # Reshape
    array_positive = array_positive.reshape(org_shape)
    array_negative = array_negative.reshape(org_shape)

    return array_positive, array_negative
