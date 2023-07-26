from typing import Tuple, TYPE_CHECKING

import numpy as np

from TidalPy.exceptions import BadArrayShape

if TYPE_CHECKING:
    from TidalPy.utilities.types import NumArray


def reshape_help(
    value: 'NumArray', comparison_shape: tuple, call_locale=None,
    force_into_new_shape: bool = False, force_ints_to_floats: bool = True
    ) -> \
        Tuple[Tuple[int, ...], np.ndarray]:
    """ Attempts to reshape value into an array that matches the shape of new_shape. Raises an error if it is unable
    to do so.

    Parameters
    ----------
    value : NumArray
        A number or array to be reshaped
    comparison_shape : tuple
        New shape to try to force the value's shape into
    call_locale :
        class that will be used to generate an error message should it be needed.
    force_into_new_shape : bool = False
        If True, the new value will be put into an array of new_shape instead of just a zero-dimensional array.
    force_ints_to_floats : bool = True
        Even if value appears to be an integer, force it to be a float.

    Returns
    -------
    new_shape : Tuple[int, ...]
        Flag if the array appears to be a new shape.
    new_array : np.ndarray
        New array based on value in the shape of new_shape.
    """

    new_shape = False

    if type(value) != np.ndarray:
        scalar = value
    else:
        if value.shape == tuple():
            # 0-D array
            scalar = value
        elif value.shape == (1,):
            scalar = value[0]
        else:
            # It has multiple values. Perhaps they are all the same?
            if np.all(value == value[0]):
                scalar = value[0]
            else:
                # Can not reshape.
                scalar = None
                if comparison_shape is not None:
                    error_text = f'Could not reshape a value to the global shape. Please provide a new value in the shape of {comparison_shape}, or reset global shape.'
                    if call_locale is not None:
                        error_text += f' Encountered at {call_locale}'
                    raise BadArrayShape(error_text)

    if scalar is None:
        dtype = None
    else:
        if type(scalar) in [complex, np.complex128]:
            dtype = np.complex128
        elif type(scalar) in [int, np.int64, np.int32]:
            if force_ints_to_floats:
                dtype = np.float64
            else:
                dtype = np.int64
        else:
            dtype = np.float64

    if comparison_shape is not None:
        if comparison_shape == tuple():
            new_array = np.asarray(scalar, dtype=dtype)
        else:
            if force_into_new_shape:
                new_array = scalar * np.ones(comparison_shape, dtype=dtype)
            else:
                new_array = np.asarray(scalar, dtype=dtype)
    else:
        if scalar is None:
            # Already an array, perhaps of a new class (or one whose global shape has been reset). Leave array alone.
            new_shape = value.shape
            new_array = value
        else:
            new_array = np.asarray(scalar, dtype=dtype)

    return new_shape, new_array
