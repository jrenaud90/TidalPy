""" Function to collapse tidal mode frequency dictionary """

from numba.core import types
from numba.typed import Dict
from numba.typed import List
from ...utilities.performance import njit

import numpy as np
MIN_MODE_FREQ = 10. * np.finfo(np.float64).eps

@njit(cacheable=True)
def collapse_modes(tidal_modes):

    # Create storage for results
    collapsed_modes = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
        )
    collapsed_mode_reference = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.unicode_type,
        )

    # Create storage for checks
    modes_to_skip = List.empty_list(
        item_type=types.unicode_type
        )
    modes_done = List.empty_list(
        item_type=types.unicode_type
        )

    # Loop through the provided tidal modes twice, merging like-valued modes.
    for mode_outer, value_outer in tidal_modes.items():

        if mode_outer in modes_to_skip:
            # Skip modes that have already been merged into a previous mode
            continue

        # If this mode has not been skipped, then it should be added into the final results.
        collapsed_modes[mode_outer] = value_outer
        modes_done.append(mode_outer)

        # Create a reference string so that the user knows which modes were merged.
        mode_str = ''

        # Loop through the modes again to find similar values
        for mode_inner, value_inner in tidal_modes.items():

            if mode_inner in modes_done:
                # Skip modes that have already been compiled into the final results
                continue

            # Check if this mode's value is similar (to machine precision) to the outer value
            if abs(value_outer - value_inner) < MIN_MODE_FREQ:

                # These two modes are very similar in value, merge them together
                modes_to_skip.append(mode_inner)

                # Update the reference string to indicate that this mode was skipped.
                if mode_str == '':
                    mode_str = mode_inner
                else:
                    mode_str += ', ' + mode_inner

        # If any modes were merged into the outer one, update the reference to reflect this
        if mode_str != '':
            collapsed_mode_reference[mode_outer] = mode_str

    return collapsed_modes, collapsed_mode_reference