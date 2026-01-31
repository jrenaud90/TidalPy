import pytest
from math import isclose

from TidalPy.utilities.lookups import IntMap3, IntMap4, IntMap3Complex, IntMap4Complex


@pytest.mark.parametrize('use_complex', (True, False))
def test_intmap3(use_complex):
    if use_complex:
        test_map = IntMap3Complex()
    else:
        test_map = IntMap3()

    # Try setting
    if use_complex:
        test_map[(1,2,3)] = 70.0 - 4.0j
    else:
        test_map[(1,2,3)] = 70.0
    assert test_map.size() == 1

    # Try getting
    if use_complex:
        assert isclose(test_map[(1,2,3)].real, 70.0)
        assert isclose(test_map[(1,2,3)].imag, -4.0)
    else:
        assert isclose(test_map[(1,2,3)], 70.0)

    # Try using negatives
    if use_complex:
        test_map[(1,-2,3)] = 167.8 + 78j
    else:
        test_map[(1,-2,3)] = 167.8
    if use_complex:
        assert isclose(test_map[(1,-2,3)].real, 167.8)
        assert isclose(test_map[(1,-2,3)].imag, 78.0)
    else:
        assert isclose(test_map[(1,-2,3)], 167.8)

    # More tests
    for i in range(10):
        test_map[(1,i,3)] = 25.2
    
    assert test_map.size() == 11
    if use_complex:
        assert isclose(test_map[(1,0,3)].real, 25.2)
        assert isclose(test_map[(1,2,3)].real, 25.2)
        assert isclose(test_map[(1,9,3)].real, 25.2)
        assert isclose(test_map[(1,0,3)].imag, 0.0)
        assert isclose(test_map[(1,2,3)].imag, 0.0)
        assert isclose(test_map[(1,9,3)].imag, 0.0)
    else:
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

    if use_complex:
        test_map.set(1, 1, 2, 45.4 - 9.4j)
        assert isclose(test_map[(1,1,2)].real, 45.4)
        assert isclose(test_map.get(1,1,2).real, 45.4)
        assert isclose(test_map[(1,1,2)].imag, -9.4)
        assert isclose(test_map.get(1,1,2).imag, -9.4)
    else:
        test_map.set(1, 1, 2, 45.4)
        assert isclose(test_map[(1,1,2)], 45.4)
        assert isclose(test_map.get(1,1,2), 45.4)


@pytest.mark.parametrize('use_complex', (True, False))
def test_intmap4(use_complex):
    if use_complex:
        test_map = IntMap4Complex()
    else:
        test_map = IntMap4()

    # Try setting
    if use_complex:
        test_map[(1,2,3,5)] = 70.0 - 4.0j
    else:
        test_map[(1,2,3,5)] = 70.0
    assert test_map.size() == 1

    # Try getting
    if use_complex:
        assert isclose(test_map[(1,2,3,5)].real, 70.0)
        assert isclose(test_map[(1,2,3,5)].imag, -4.0)
    else:
        assert isclose(test_map[(1,2,3,5)], 70.0)

    # Try using negatives
    if use_complex:
        test_map[(1,-2,3,-5)] = 167.8 + 78j
    else:
        test_map[(1,-2,3,-5)] = 167.8
    if use_complex:
        assert isclose(test_map[(1,-2,3,-5)].real, 167.8)
        assert isclose(test_map[(1,-2,3,-5)].imag, 78.0)
    else:
        assert isclose(test_map[(1,-2,3,-5)], 167.8)

    # More tests
    for i in range(10):
        test_map[(1,i,3,5)] = 25.2
    
    assert test_map.size() == 11
    if use_complex:
        assert isclose(test_map[(1,0,3,5)].real, 25.2)
        assert isclose(test_map[(1,2,3,5)].real, 25.2)
        assert isclose(test_map[(1,9,3,5)].real, 25.2)
        assert isclose(test_map[(1,0,3,5)].imag, 0.0)
        assert isclose(test_map[(1,2,3,5)].imag, 0.0)
        assert isclose(test_map[(1,9,3,5)].imag, 0.0)
    else:
        assert isclose(test_map[(1,0,3,5)], 25.2)
        assert isclose(test_map[(1,2,3,5)], 25.2)
        assert isclose(test_map[(1,9,3,5)], 25.2)

    # Check unset value
    with pytest.raises(KeyError):
        result = test_map[(2,2,3,5)]

    # Check other methods
    test_map.clear()
    assert test_map.size() == 0
    test_map.reserve(10)
    assert test_map.size() == 0
    with pytest.raises(KeyError):
        result = test_map[(1,2,3,5)]

    if use_complex:
        test_map.set(1, 1, 2, 5, 45.4 - 9.4j)
        assert isclose(test_map[(1,1,2,5)].real, 45.4)
        assert isclose(test_map.get(1,1,2,5).real, 45.4)
        assert isclose(test_map[(1,1,2,5)].imag, -9.4)
        assert isclose(test_map.get(1,1,2,5).imag, -9.4)
    else:
        test_map.set(1, 1, 2, 5, 45.4)
        assert isclose(test_map[(1,1,2,5)], 45.4)
        assert isclose(test_map.get(1,1,2,5), 45.4)

    assert test_map.size() == 1