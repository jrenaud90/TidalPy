import numpy as np

import TidalPy


from TidalPy.utilities.conversions.timing import convert_to_hms

def test_convert_hms():
    """ Test the convert_to_hms function """

    result = convert_to_hms(60.)
    assert type(result) == tuple

    day, hour, minute, second = result
    assert type(day) in [int, np.int32]
    assert type(hour) in [int, np.int32]
    assert type(minute) in [int, np.int32]
    assert type(second) in [float, np.float64]

    assert day == 0
    assert hour == 0
    assert minute == 1
    assert second == 0.

    day, hour, minute, second = convert_to_hms(198000.)
    assert day == 2
    assert hour == 7
    assert minute == 0
    assert second == 0.

    day, hour, minute, second = convert_to_hms(198025.5)
    assert day == 2
    assert hour == 7
    assert minute == 0
    assert second == 25.5
