"""Tests for the integer-keyed maps in ``TidalPy.Utilities.lookups``."""
from math import isclose

import pytest

from TidalPy.Utilities.lookups import (
    IntMap1, IntMap2, IntMap3, IntMap4,
    IntMap1Complex, IntMap2Complex, IntMap3Complex, IntMap4Complex,
)

# (map class, key size, value type, values set by item, by a negative key, and by .set)
_CASES = [
    (IntMap1, 1, float, (70.0, 167.8, 45.4)),
    (IntMap2, 2, float, (70.0, 167.8, 45.4)),
    (IntMap3, 3, float, (70.0, 167.8, 45.4)),
    (IntMap4, 4, float, (70.0, 167.8, 45.4)),
    (IntMap1Complex, 1, complex, (70.0 - 4.0j, 167.8 + 78j, 45.4 - 9.4j)),
    (IntMap2Complex, 2, complex, (70.0 - 4.0j, 167.8 + 78j, 45.4 - 9.4j)),
    (IntMap3Complex, 3, complex, (70.0 - 4.0j, 167.8 + 78j, 45.4 - 9.4j)),
    (IntMap4Complex, 4, complex, (70.0 - 4.0j, 167.8 + 78j, 45.4 - 9.4j)),
]


def _assert_value(got, expected):
    if isinstance(expected, complex):
        assert isclose(got.real, expected.real)
        assert isclose(got.imag, expected.imag)
    else:
        assert isclose(got, expected)


@pytest.mark.parametrize(
    "map_class, key_size, value_type, values",
    _CASES,
    ids=[case[0].__name__ for case in _CASES])
def test_intmap(map_class, key_size, value_type, values):
    """Setting, getting, negative keys, missing keys, clear, reserve, set/get, len, and iteration work."""
    item_value, negative_key_value, set_value = values
    test_map = map_class()
    test_key = [index + 2 for index in range(key_size)]
    # The varied key slot is the second one when there is one.
    varied_slot = 0 if key_size == 1 else 1

    test_map[tuple(test_key)] = item_value
    assert test_map.size() == 1
    _assert_value(test_map[tuple(test_key)], item_value)

    test_key[varied_slot] = -2
    test_map[tuple(test_key)] = negative_key_value
    _assert_value(test_map[tuple(test_key)], negative_key_value)

    for index in range(10):
        test_key[varied_slot] = index
        test_map[tuple(test_key)] = 25.2
        _assert_value(test_map[tuple(test_key)], value_type(25.2))
    assert test_map.size() == 11

    test_key[0] = 1999
    with pytest.raises(KeyError):
        test_map[tuple(test_key)]

    test_map.clear()
    assert test_map.size() == 0
    test_map.reserve(10)
    assert test_map.size() == 0
    test_key = [index + 2 for index in range(key_size)]
    with pytest.raises(KeyError):
        test_map[tuple(test_key)]

    test_map.set(tuple(test_key), set_value)
    _assert_value(test_map[tuple(test_key)], set_value)
    _assert_value(test_map.get(tuple(test_key)), set_value)

    assert len(test_map) == test_map.size()

    for key_tuple, result in test_map:
        assert isinstance(key_tuple, tuple)
        assert len(key_tuple) == key_size
        assert isinstance(result, value_type)
