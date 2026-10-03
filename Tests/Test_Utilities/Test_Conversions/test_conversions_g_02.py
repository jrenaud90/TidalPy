"""Tests for the gravitational constant option and input checks of the Kepler conversions."""
import math

import pytest

from TidalPy.constants import G
from TidalPy.Utilities.conversions import orbital_motion2semi_a, semi_a2orbital_motion

_HOST_MASS = 1.989e30
_TARGET_MASS = 5.972e24
_ORBITAL_MOTION = 2.0 * math.pi / (365.25 * 86400.0)
_SEMI_MAJOR_AXIS = 1.495978707e11


def _expected_semi_major_axis(gravitational_constant):
    return (gravitational_constant * (_HOST_MASS + _TARGET_MASS) / _ORBITAL_MOTION ** 2) ** (1.0 / 3.0)


def _expected_orbital_motion(gravitational_constant):
    return math.sqrt(gravitational_constant * (_HOST_MASS + _TARGET_MASS) / _SEMI_MAJOR_AXIS ** 3)


@pytest.mark.parametrize(
    "g_kwargs, gravitational_constant",
    [
        pytest.param({}, G, id="default"),
        pytest.param({"G_to_use": 6.674e-11}, 6.674e-11, id="6.674e-11"),
        pytest.param({"G_to_use": 2.0 * G}, 2.0 * G, id="2G"),
        pytest.param({"G_to_use": 1}, 1, id="1"),
    ])
def test_G_used_in_conversions(g_kwargs, gravitational_constant):
    """Both conversions use the package G by default and an explicit ``G_to_use`` otherwise."""
    assert math.isclose(
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS, **g_kwargs),
        _expected_semi_major_axis(gravitational_constant),
        rel_tol=1e-14)
    assert math.isclose(
        semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, _TARGET_MASS, **g_kwargs),
        _expected_orbital_motion(gravitational_constant),
        rel_tol=1e-14)


def test_none_G_uses_the_config():
    """``G_to_use=None`` gives the same result as the default and as the config's G."""
    assert orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS, G_to_use=None) == \
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS)
    assert semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, _TARGET_MASS, G_to_use=None) == \
        semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, _TARGET_MASS)
    assert math.isclose(
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS, G_to_use=None),
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS, G_to_use=G),
        rel_tol=1e-15)


def test_explicit_G_changes_the_result():
    """Doubling G scales the semi-major axis by the cube root of two."""
    default_axis = orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS)
    doubled_axis = orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, G_to_use=2.0 * G)
    assert math.isclose(doubled_axis / default_axis, 2.0 ** (1.0 / 3.0), rel_tol=1e-14)


@pytest.mark.parametrize("gravitational_constant", [0.0, -G, math.nan])
def test_non_positive_G_raises(gravitational_constant):
    """A non-positive or NaN ``G_to_use`` raises ValueError in both conversions."""
    with pytest.raises(ValueError, match="G_to_use"):
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, G_to_use=gravitational_constant)
    with pytest.raises(ValueError, match="G_to_use"):
        semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, G_to_use=gravitational_constant)


@pytest.mark.parametrize(
    "convert, value, match",
    [
        (orbital_motion2semi_a, 0.0, "Orbital motion"),
        (orbital_motion2semi_a, -1.0e-5, "Orbital motion"),
        (orbital_motion2semi_a, math.nan, "Orbital motion"),
        (semi_a2orbital_motion, 0.0, "Semi-major axis"),
        (semi_a2orbital_motion, -1.0e8, "Semi-major axis"),
        (semi_a2orbital_motion, math.nan, "Semi-major axis"),
    ],
    ids=["motion_zero", "motion_negative", "motion_nan", "axis_zero", "axis_negative", "axis_nan"])
def test_non_positive_orbital_element_raises(convert, value, match):
    """A non-positive or NaN orbital motion or semi-major axis raises ValueError."""
    with pytest.raises(ValueError, match=match):
        convert(value, _HOST_MASS)
