import pytest
from math import isclose

from TidalPy.utilities.lookups import IntMap3

def test_intmap3():
    test_map = IntMap3()

    # Try setting
    test_map[(1,2,3)] = 70.0
    assert test_map.size() == 1

    # Try getting
    assert isclose(test_map[(1,2,3)], 70.0)

    # More tests
    for i in range(10):
        test_map[(1,i,3)] = 25.2
    
    assert test_map.size() == 10
    assert isclose(test_map[(1,0,3)], 25.2)
    assert isclose(test_map[(1,2,3)], 25.2)
    assert isclose(test_map[(1,9,3)], 25.2)

    # Check unset value
    with pytest.raises(KeyError):
        result = test_map[(2,2,3)]

    # Check other methods
    test_map.clear()
    assert test_map.size() == 0
    test_map.reserve(10)
    assert test_map.size() == 0
    with pytest.raises(KeyError):
        result = test_map[(1,2,3)]

    test_map.set(1, 1, 2, 45.4)
    assert isclose(test_map[(1,1,2)], 45.4)
    assert isclose(test_map.get(1,1,2), 45.4)
    assert test_map.size() == 1

test_intmap3()
