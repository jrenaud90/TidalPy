import numpy as np

from ...exceptions import BadArrayShape
from ...types import NumArray

def reshape_help(value: NumArray, new_shape: tuple, call_locale = None) -> np.ndarray:
    """ Attempts to reshape value into an array that matches the shape of new_shape. Raises an error if it is unable
    to do so.

    Parameters
    ----------
    value : NumArray
        A number or array to be reshaped
    new_shape : tuple
        New shape to try to force the value's shape into
    call_locale :
        class that will be used to generate an error message should it be needed.

    Returns
    -------
    new_array : np.ndarray
        New array based on value in the shape of new_shape.
    """

    if type(value) != np.ndarray:
        if new_shape is None:
            # No shape provided, just make the float an array
            new_array = np.asarray(value)
        else:
            new_array = value * np.ones(new_shape)
    else:
        # It is an array.
        if new_shape is not None:
            # Is it a 0-D or 1-D array
            if value.shape == tuple():
                new_array = value * np.ones(new_shape)
            elif value.shape == (1,):
                new_array = value[0] * np.ones(new_shape)
            else:
                # It has multiple values. Perhaps they are all the same?
                if np.all(value == value[0]):
                    new_array = value[0] * np.ones(new_shape)
                else:
                    # Can not reshape.
                    error_text = 'Could not reshape a value to the global shape.'
                    if call_locale is not None:
                        error_text += f' Encountered at {call_locale}'
                    raise BadArrayShape(error_text)
        else:
            # Can't do much with the info provided. Simply pass the new value back.
            new_array = value

    return new_array