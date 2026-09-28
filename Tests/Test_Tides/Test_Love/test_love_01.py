"""``LoveNumbers``: storage, equality, iteration, to_dict, repr, and access through a PhysicsLayer."""
import pytest

from TidalPy.Tides.love.love import LoveNumbers

_COMPLEX = (0.3 - 0.01j, 0.6 - 0.02j, 0.1 - 0.005j)
_DICT_KEYS = ("love_number_k_re", "love_number_k_im", "love_number_h_re", "love_number_h_im",
              "love_number_l_re", "love_number_l_im")


def _make(values):
    k, h, l = values
    return LoveNumbers(k=k, h=h, l=l)


@pytest.fixture
def physics_layer():
    physics = pytest.importorskip("TidalPy.Structures.layers.physics")
    return physics.PhysicsLayer(
        "test",
        0,
        0.0,
        1e6,
        1e20,
        love_number_k=_COMPLEX[0],
        love_number_h=_COMPLEX[1],
        love_number_l=_COMPLEX[2])


_VALUE_CASES = pytest.mark.parametrize("values", [
    (0 + 0j, 0 + 0j, 0 + 0j),
    (0.3, 0.6, 0.1),
    _COMPLEX,
    (1.0 + 2.0j, 3.0 + 4.0j, 5.0 + 6.0j),
    (1.0j, 2.0j, 3.0j),
], ids=["zero", "real", "complex", "complex_large", "imaginary"])


@_VALUE_CASES
def test_love_numbers_store_values(values):
    """k, h, l are stored as complex."""
    love = _make(values)
    for found, expected in zip((love.k, love.h, love.l), values):
        assert isinstance(found, complex)
        assert found == pytest.approx(complex(expected))


@_VALUE_CASES
def test_love_numbers_iterate_and_unpack(values):
    """Iteration and unpacking yield (k, h, l) in order."""
    love = _make(values)
    as_list = list(love)
    assert len(as_list) == 3
    for found, expected in zip(as_list, values):
        assert found == pytest.approx(complex(expected))
    k, h, l = love
    assert (k, h, l) == pytest.approx(tuple(complex(value) for value in values))


@_VALUE_CASES
def test_love_numbers_to_dict(values):
    """to_dict holds the real and imaginary part of each Love number."""
    as_dict = _make(values).to_dict()
    for key in _DICT_KEYS:
        assert key in as_dict, f"Missing key: {key}"
    for name, value in zip("khl", values):
        assert as_dict[f"love_number_{name}_re"] == pytest.approx(complex(value).real)
        assert as_dict[f"love_number_{name}_im"] == pytest.approx(complex(value).imag)
    if not any(values):
        for key in as_dict:
            assert as_dict[key] == pytest.approx(0.0), f"{key} should be 0.0"


@pytest.mark.parametrize("left, right, expected", [
    pytest.param(_COMPLEX, None, True, id="self"),
    pytest.param(_COMPLEX, _COMPLEX, True, id="copy"),
    pytest.param((0.3, 0.6, 0.1), (0.4, 0.6, 0.1), False, id="different"),
])
def test_love_numbers_equality(left, right, expected):
    """Equality compares values; an object equals itself."""
    first = _make(left)
    second = first if right is None else _make(right)
    assert (first == second) is expected
    assert (first != second) is not expected


def test_love_numbers_eq_not_implemented_for_non_love():
    assert _make((0j, 0j, 0j)).__eq__(42) is NotImplemented


def test_love_numbers_repr():
    """repr names the class and each Love number."""
    text = repr(_make((0.3 - 0.01j, 0.6, 0.0)))
    assert "LoveNumbers" in text
    assert "k=" in text
    assert "h=" in text
    assert "l=" in text


def test_love_numbers_from_physics_layer(physics_layer):
    """PhysicsLayer.love_numbers holds the layer's Love numbers."""
    love = physics_layer.love_numbers
    assert isinstance(love, LoveNumbers)
    assert love.k == pytest.approx(_COMPLEX[0])
    assert love.h == pytest.approx(_COMPLEX[1])
    assert love.l == pytest.approx(_COMPLEX[2])


def test_love_numbers_from_physics_layer_to_dict(physics_layer):
    """PhysicsLayer.love_numbers.to_dict matches the layer's config dict."""
    as_dict = physics_layer.love_numbers.to_dict()
    config = physics_layer.get_config_dict()
    for key in _DICT_KEYS:
        assert as_dict[key] == pytest.approx(config[key]), f"Mismatch for {key}"
