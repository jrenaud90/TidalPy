""" Function to collapse tidal mode frequency dictionary """

from typing import Dict, Tuple

import numpy as np
from numba.core import types
from numba import typed

from TidalPy.utilities.performance import njit

MIN_MODE_FREQ = 10. * np.finfo(np.float64).eps

CollapseModeValueType = types.Tuple(types=[types.float64,
                                           types.int64,
                                           types.unicode_type])


@njit(cacheable=True)
def collapse_modes(tidal_modes: Dict[str, float], unique_frequencies: bool = False) -> \
        Dict[str, Tuple[float, int, str]]:
    """ Takes a (numba) dictionary of tidal modes and returns a dictionary of unique modes or frequencies

    Parameters
    ----------
    tidal_modes : Dict[str, float]
        Dictionary of tidal modes stored like {'n-o': 0.005} where the number is in rad s-1.
        If numba is being used then this must be passed as a numba TypedDict of the form:
             Dict.empty(
                key_type=types.unicode_type,
                value_type=types.float64,
            )
    unique_frequencies : bool = False
        If True, then the absolute value of the mode will be used instead of the mode itself for comparison.

    Returns
    -------
    collapsed_modes : Dict[str, Tuple[float, int, str]]
        Dictionary of collapsed modes (or frequencies). Along with the value, the number of instances it appeared,
            and a string of which modes were not unique is also provided.
        Example:
            {
            '2n' : (0.005, 1, ''),            # Mode was unique.
            'n-o': (0.110, 3, '2n-2o, 3n-3o') # Two other modes were found to have a similar value.
            }
    """
    # Create copy of input dict
    tidal_modes_copy = typed.Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
    )
    for mode_name, mode_value in tidal_modes.items():
        # Convert to frequency if requested by the user
        if unique_frequencies:
            mode_value = abs(mode_value)

        tidal_modes_copy[mode_name] = mode_value

    # Create storage for results
    collapsed_modes = typed.Dict.empty(
        key_type=types.unicode_type,
        value_type=CollapseModeValueType,
    )

    # Create storage for checks
    modes_to_skip = typed.List.empty_list(
        item_type=types.unicode_type
    )

    # Loop through the provided tidal modes twice, merging like-valued modes.
    for mode_outer, value_outer in tidal_modes.items():

        if mode_outer in modes_to_skip:
            # Skip modes that have already been merged into a previous mode
            continue

        # Convert to frequency if requested by the user
        if unique_frequencies:
            value_outer = abs(value_outer)

        # If this mode has not been skipped, then it should be added into the final results.
        # Remove this mode from the copied dict so that it is not looked at in the inner loop.
        del tidal_modes_copy[mode_outer]

        # Create a reference string so that the user knows which modes were merged.
        mode_str = ''
        mode_count = 1

        # Loop through the modes again to find similar values
        for mode_inner, value_inner in tidal_modes_copy.items():

            # Check if this mode's value is similar (to machine precision) to the outer value
            if abs(value_outer - value_inner) < MIN_MODE_FREQ:

                # These two modes are very similar in value, merge them together
                modes_to_skip.append(mode_inner)
                mode_count += 1

                # Update the reference string to indicate that this mode was skipped.
                if mode_str == '':
                    mode_str = mode_inner
                else:
                    mode_str += ', ' + mode_inner

        # If any modes were merged into the outer one, update the reference to reflect this
        collapsed_modes[mode_outer] = (value_outer, mode_count, mode_str)

    return collapsed_modes
