import pytest
from math import isclose

from TidalPy.utilities.lookups import (
    IntMap1, IntMap2, IntMap3, IntMap4,
    IntMap1Complex, IntMap2Complex, IntMap3Complex, IntMap4Complex,
)

class_dict = {
    'IntMap1': IntMap1,
    'IntMap2': IntMap2,
    'IntMap3': IntMap3,
    'IntMap4': IntMap4,
    'IntMap1Complex': IntMap1Complex,
    'IntMap2Complex': IntMap2Complex,
    'IntMap3Complex': IntMap3Complex,
    'IntMap4Complex': IntMap4Complex,
}


@pytest.mark.parametrize('use_complex', (True, False))
@pytest.mark.parametrize('key_size', (1, 2, 3, 4))
def test_intmap(use_complex, key_size):
    if use_complex:
        test_map = class_dict[f'IntMap{key_size}Complex']()
    else:
        test_map = class_dict[f'IntMap{key_size}']()

    test_key = [i + 2 for i in range(key_size)]    

    # Try setting
    if use_complex:
        test_map[tuple(test_key)] = 70.0 - 4.0j
    else:
        test_map[tuple(test_key)] = 70.0
    assert test_map.size() == 1

    # Try getting
    if use_complex:
        assert isclose(test_map[tuple(test_key)].real, 70.0)
        assert isclose(test_map[tuple(test_key)].imag, -4.0)
    else:
        assert isclose(test_map[tuple(test_key)], 70.0)

    # Try using negatives
    if key_size == 1:
        test_key[0] = -2
    else:
        test_key[1] = -2
    
    if use_complex:
        test_map[tuple(test_key)] = 167.8 + 78j
    else:
        test_map[tuple(test_key)] = 167.8
    if use_complex:
        assert isclose(test_map[tuple(test_key)].real, 167.8)
        assert isclose(test_map[tuple(test_key)].imag, 78.0)
    else:
        assert isclose(test_map[tuple(test_key)], 167.8)

    # More tests
    for i in range(10):
        if key_size == 1:
            test_key[0] = i
        else:
            test_key[1] = i
        test_map[tuple(test_key)] = 25.2

        if use_complex:
            assert isclose(test_map[tuple(test_key)].real, 25.2)
            assert isclose(test_map[tuple(test_key)].imag, 0.0)
    
    assert test_map.size() == 11

    # Check unset value
    with pytest.raises(KeyError):
        test_key[0] = 1999
        result = test_map[tuple(test_key)]

    # Check other methods
    test_map.clear()
    assert test_map.size() == 0
    test_map.reserve(10)
    assert test_map.size() == 0
    test_key = [i + 2 for i in range(key_size)]
    with pytest.raises(KeyError):
        result = test_map[tuple(test_key)]

    if use_complex:
        test_map.set(tuple(test_key), 45.4 - 9.4j)
        assert isclose(test_map[tuple(test_key)].real, 45.4)
        assert isclose(test_map.get(tuple(test_key)).real, 45.4)
        assert isclose(test_map[tuple(test_key)].imag, -9.4)
        assert isclose(test_map.get(tuple(test_key)).imag, -9.4)
    else:
        test_map.set(tuple(test_key), 45.4)
        assert isclose(test_map[tuple(test_key)], 45.4)
        assert isclose(test_map.get(tuple(test_key)), 45.4)
    
    # Check that __len__ is working
    assert len(test_map) == test_map.size()

    # Check that the iterator is working
    for key_tuple, result in test_map:
        assert isinstance(key_tuple, tuple)
        assert len(key_tuple) == key_size
        if use_complex:
            assert isinstance(result, complex)
        else:
            assert isinstance(result, float)
