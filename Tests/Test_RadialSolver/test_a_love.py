"""LoveNumbers storage, quality factor, and phase lag."""
import math

import pytest

from TidalPy.RadialSolver.love import LoveNumbers

COMPONENTS = ('k', 'h', 'l')


def test_love_numbers_initialization():
    """Complex k, h, and l are stored and returned unchanged."""
    love_k = 1 + 2j
    love_h = 3 + 4j
    love_l = 5 + 6j
    love = LoveNumbers(love_k, love_h, love_l)
    assert love.k == love_k
    assert love.h == love_h
    assert love.l == love_l


@pytest.mark.parametrize('component', COMPONENTS)
@pytest.mark.parametrize('value, expected', (
    (3 + 4j, -abs(3 + 4j) / 4.0),
    (5 + 0j, math.inf),
), ids=('normal', 'zero_imag'))
def test_love_numbers_Q(value, expected, component):
    """Q = -|x| / Im(x), infinite when Im(x) = 0."""
    love = LoveNumbers(value, value, value)
    assert math.isclose(getattr(love, f'Q_{component}'), expected)


@pytest.mark.parametrize('component', COMPONENTS)
@pytest.mark.parametrize('value, expected', (
    (3 + 4j, math.atan(-4.0 / 3.0)),
    (5 + 0j, 0.0),
    (0 + 4j, math.pi / 2.0),
), ids=('normal', 'zero_imag', 'zero_real'))
def test_love_numbers_lag(value, expected, component):
    """lag = atan(-Im(x) / Re(x)), 0 when Im(x) = 0 and pi/2 when Re(x) = 0."""
    love = LoveNumbers(value, value, value)
    assert math.isclose(getattr(love, f'lag_{component}'), expected)
