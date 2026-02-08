import pytest
import math

def test_mode_map():
    """Tests `ModeMap` a Cython wrapper of `c_ModeMap`."""
    from TidalPy.Tides_x.potential import ModeMap, test_mode_map

    mode_map, total_size = test_mode_map()
    assert isinstance(mode_map, ModeMap)
    assert mode_map.size() == total_size
    assert len(mode_map) == total_size

    # Test Getters
    for l in range(2, 4):
        for m in range(0, l + 1):
            for p in range(0, l + 1):
                for q in range(-1, 2):
                    result_1 = mode_map[(l, m, p, q)]
                    result_2 = mode_map.get((l, m, p, q))
                    assert result_1 == result_2
                    assert isinstance(result_1, tuple)
                    assert isinstance(result_1[0], float)
                    assert math.isclose(result_1[0], (l - 2 * p + q) * 10.0 - m * 5.0)
                    assert isinstance(result_1[1], float)
                    assert math.isclose(result_1[1], l + p + q + m - 4.5)
                    assert isinstance(result_1[2], int)
                    assert result_1[2] == (l - 2 * p + q)
                    assert isinstance(result_1[3], int)
                    assert result_1[3] == -m
    
    # Test Iters
    for (l, m, p, q), result in mode_map:
        assert isinstance(result, tuple)
        assert isinstance(result[0], float)
        assert math.isclose(result[0], (l - 2 * p + q) * 10.0 - m * 5.0)
        assert isinstance(result[1], float)
        assert math.isclose(result[1], l + p + q + m - 4.5)
        assert isinstance(result[2], int)
        assert result[2] == (l - 2 * p + q)
        assert isinstance(result[3], int)
        assert result[3] == -m
    
    # Test setter
    mode_map.set((100, 100, 10, 10), (1.5, 2.5, 2, 1))
    assert mode_map[(100, 100, 10, 10)] == (1.5, 2.5, 2, 1)
    mode_map[(200, 200, 10, 10)] = (2.5, 3.5, 4, 5)
    assert mode_map[(200, 200, 10, 10)] == (2.5, 3.5, 4, 5)
