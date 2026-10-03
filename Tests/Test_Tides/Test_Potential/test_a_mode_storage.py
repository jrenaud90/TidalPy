"""``ModeMap``: setting, getting, and iterating tidal mode entries."""
import math

from TidalPy.Tides.potential import ModeMap

_KEYS = [(l, m, p, q) for l in range(2, 4) for m in range(0, l + 1) for p in range(0, l + 1) for q in range(-1, 2)]


def _entry(l, m, p, q):
    return (float((l - 2 * p + q) * 10.0 - m * 5.0), float(l + p + q + m - 4.5), l - 2 * p + q, -m)


def _check_entry(key, result):
    l, m, p, q = key
    assert isinstance(result, tuple)
    assert isinstance(result[0], float)
    assert math.isclose(result[0], (l - 2 * p + q) * 10.0 - m * 5.0)
    assert isinstance(result[1], float)
    assert math.isclose(result[1], l + p + q + m - 4.5)
    assert isinstance(result[2], int)
    assert result[2] == (l - 2 * p + q)
    assert isinstance(result[3], int)
    assert result[3] == -m


def test_mode_map():
    """Entries round-trip through item access, get, iteration, and set."""
    mode_map = ModeMap()
    for key in _KEYS:
        mode_map[key] = _entry(*key)
    assert isinstance(mode_map, ModeMap)
    assert mode_map.size() == len(_KEYS)
    assert len(mode_map) == len(_KEYS)

    for key in _KEYS:
        result = mode_map[key]
        assert result == mode_map.get(key)
        _check_entry(key, result)

    for key, result in mode_map:
        _check_entry(key, result)

    mode_map.set((100, 100, 10, 10), (1.5, 2.5, 2, 1))
    assert mode_map[(100, 100, 10, 10)] == (1.5, 2.5, 2, 1)
    mode_map[(200, 200, 10, 10)] = (2.5, 3.5, 4, 5)
    assert mode_map[(200, 200, 10, 10)] == (2.5, 3.5, 4, 5)
