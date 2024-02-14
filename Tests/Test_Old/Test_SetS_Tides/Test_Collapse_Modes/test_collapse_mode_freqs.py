import TidalPy

import numpy as np
from TidalPy.tides.modes.collapse_modes import collapse_modes
from numba.core import types
from numba.typed import Dict
from numba.typed.typeddict import Dict as TypedDict
from numba.typed import List


def test_collapse_mode():
    """ Test the collapse tidal mode function with the frequency switch off """

    # This function requires numba typed dict as input
    # The Dict.empty() constructs a typed dictionary.
    # The key and value typed must be explicitly declared.
    test_modes = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
        )

    for key, val in {
        'n'  : 0.2231,
        '2n' : 0.2121,
        '3n' : 0.11,
        '4n' : 0.1,
        '5n' : -0.2231,
        '6n' : 0.3232231,
        '7n' : 0.222221231,
        '8n' : 0.22310000000000000000,
        '9n' : 0.72231,
        '10n': 0.382231,
        '11n': 0.942231,
        '12n': 0.2231,
        '13n': 0.1235622231,
        }.items():
        test_modes[key] = val

    collapsed_modes = collapse_modes(test_modes, unique_frequencies=False)

    # Check types
    assert type(collapsed_modes) in [dict, TypedDict]
    for mode_name, mode_data in collapsed_modes.items():
        assert type(mode_data) == tuple or isinstance(mode_data, types.Tuple)
        assert len(mode_data) == 3
        mode_value, mode_multiplier, mode_ref_str = mode_data
        assert type(mode_value) in (float, np.float64, types.float64)
        assert type(mode_multiplier) in (int, np.int64, types.int64)
        assert type(mode_ref_str) == str

        if mode_ref_str != '':
            num_modes = len(mode_ref_str.split(', ')) + 1
            assert mode_multiplier == num_modes

        if mode_name == 'n':
            assert mode_multiplier == 3
            assert '8n' in mode_ref_str
            assert '12n' in mode_ref_str

def test_collapse_mode_unique_freqs():
    """ Test the collapse tidal mode function with the frequency switch on """

    # This function requires numba typed dict as input
    # The Dict.empty() constructs a typed dictionary.
    # The key and value typed must be explicitly declared.
    test_modes = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
        )

    for key, val in {
        'n'  : 0.2231,
        '2n' : 0.2121,
        '3n' : 0.11,
        '4n' : 0.1,
        '5n' : -0.2231,
        '6n' : 0.3232231,
        '7n' : 0.222221231,
        '8n' : 0.22310000000000000000,
        '9n' : 0.72231,
        '10n': 0.382231,
        '11n': 0.942231,
        '12n': 0.2231,
        '13n': 0.1235622231,
        }.items():
        test_modes[key] = val

    collapsed_modes = collapse_modes(test_modes, unique_frequencies=True)

    # Check types
    assert type(collapsed_modes) in [dict, TypedDict]
    for mode_name, mode_data in collapsed_modes.items():
        assert type(mode_data) == tuple or isinstance(mode_data, types.Tuple)
        assert len(mode_data) == 3
        mode_value, mode_multiplier, mode_ref_str = mode_data
        assert type(mode_value) in (float, np.float64, types.float64)
        assert type(mode_multiplier) in (int, np.int64, types.int64)
        assert type(mode_ref_str) == str

        if mode_ref_str != '':
            num_modes = len(mode_ref_str.split(', ')) + 1
            assert mode_multiplier == num_modes

        if mode_name == 'n':
            assert mode_multiplier == 4
            assert '8n' in mode_ref_str
            assert '12n' in mode_ref_str
            assert '5n' in mode_ref_str
