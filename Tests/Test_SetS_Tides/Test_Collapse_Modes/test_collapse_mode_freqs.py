import TidalPy
TidalPy.test_mode()

from TidalPy.tides.modes.collapse_modes import collapse_modes
from numba.core import types
from numba.typed import Dict
from numba.typed import List


def test_collapse_mode():
    """ Test the collapse mode frequency function """

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

    collapsed_modes, collapsed_mode_reference = collapse_modes(test_modes)

    # Check types
    assert type(collapsed_modes) is Dict
    assert type(collapsed_mode_reference) is Dict

    # There should be 1 mode that has two values
    assert len(collapsed_mode_reference) == 1
    assert len(collapsed_modes) == (len(test_modes) - 2)
    assert 'n' in collapsed_mode_reference
    assert type(collapsed_mode_reference['n']) == str
    assert '8n' in collapsed_mode_reference['n']
    assert '12n' in collapsed_mode_reference['n']


