"""Gravitational constant and input checks of the Kepler conversions in ``TidalPy.Utilities_x.conversions``.

``G_to_use`` defaults to ``None``, which reads the TidalPy config's gravitational constant at call time; an explicit
value overrides it. A non-positive orbital motion or semi-major axis, or a non-positive ``G_to_use``, raises
``ValueError``.
"""

import math

import pytest

from TidalPy.constants import G
from TidalPy.Utilities_x.conversions import orbital_motion2semi_a, semi_a2orbital_motion

_HOST_MASS = 1.989e30
_TARGET_MASS = 5.972e24
_ORBITAL_MOTION = 2.0 * math.pi / (365.25 * 86400.0)
_SEMI_MAJOR_AXIS = 1.495978707e11


def _expected_semi_major_axis(gravitational_constant):
    return (gravitational_constant * (_HOST_MASS + _TARGET_MASS) / _ORBITAL_MOTION ** 2) ** (1.0 / 3.0)


def _expected_orbital_motion(gravitational_constant):
    return math.sqrt(gravitational_constant * (_HOST_MASS + _TARGET_MASS) / _SEMI_MAJOR_AXIS ** 3)


# =====================================================================================================================
# Gravitational constant
# =====================================================================================================================
def test_default_G_is_the_package_constant():
    assert math.isclose(
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS), _expected_semi_major_axis(G),
        rel_tol=1e-14)
    assert math.isclose(
        semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, _TARGET_MASS), _expected_orbital_motion(G),
        rel_tol=1e-14)


def test_none_G_uses_the_config():
    assert orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS, G_to_use=None) == \
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS)
    assert semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, _TARGET_MASS, G_to_use=None) == \
        semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, _TARGET_MASS)
    assert math.isclose(
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS, G_to_use=None),
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS, G_to_use=G),
        rel_tol=1e-15)


@pytest.mark.parametrize("gravitational_constant", [6.674e-11, 2.0 * G, 1])
def test_explicit_G_overrides_the_config(gravitational_constant):
    assert math.isclose(
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, _TARGET_MASS, G_to_use=gravitational_constant),
        _expected_semi_major_axis(gravitational_constant),
        rel_tol=1e-14)
    assert math.isclose(
        semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, _TARGET_MASS, G_to_use=gravitational_constant),
        _expected_orbital_motion(gravitational_constant),
        rel_tol=1e-14)


def test_explicit_G_changes_the_result():
    default_axis = orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS)
    doubled_axis = orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, G_to_use=2.0 * G)
    assert math.isclose(doubled_axis / default_axis, 2.0 ** (1.0 / 3.0), rel_tol=1e-14)


@pytest.mark.parametrize("gravitational_constant", [0.0, -G, math.nan])
def test_non_positive_G_raises(gravitational_constant):
    with pytest.raises(ValueError, match="G_to_use"):
        orbital_motion2semi_a(_ORBITAL_MOTION, _HOST_MASS, G_to_use=gravitational_constant)
    with pytest.raises(ValueError, match="G_to_use"):
        semi_a2orbital_motion(_SEMI_MAJOR_AXIS, _HOST_MASS, G_to_use=gravitational_constant)


# =====================================================================================================================
# Orbital elements
# =====================================================================================================================
@pytest.mark.parametrize("orbital_motion", [0.0, -1.0e-5, math.nan])
def test_non_positive_orbital_motion_raises(orbital_motion):
    """A zero frequency would give an infinite orbit and a negative one no orbit."""
    with pytest.raises(ValueError, match="Orbital motion"):
        orbital_motion2semi_a(orbital_motion, _HOST_MASS)


@pytest.mark.parametrize("semi_major_axis", [0.0, -1.0e8, math.nan])
def test_non_positive_semi_major_axis_raises(semi_major_axis):
    with pytest.raises(ValueError, match="Semi-major axis"):
        semi_a2orbital_motion(semi_major_axis, _HOST_MASS)
