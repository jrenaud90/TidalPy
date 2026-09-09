"""Tests for ``TidalPy.Utilities_x.dimensions.nondimensional`` (non-dimensionalization scales).

Ported from the classic ``TidalPy.utilities.dimensions`` tests; the two implementations must agree
while both exist.
"""

from math import isclose, isnan

from TidalPy.constants import G, pi
from TidalPy.Utilities_x.dimensions import NonDimensionalScalesClass, build_nondimensional_scales

_FREQUENCY = 1.0e-3
_MEAN_RADIUS = 1.0e6
_BULK_DENSITY = 5500.0


def test_non_dimensionalize_structure_initializes_nan():
    test_struct = NonDimensionalScalesClass()
    assert isnan(test_struct.second2_conversion)
    assert isnan(test_struct.second_conversion)
    assert isnan(test_struct.length_conversion)
    assert isnan(test_struct.length3_conversion)
    assert isnan(test_struct.density_conversion)
    assert isnan(test_struct.mass_conversion)
    assert isnan(test_struct.pascal_conversion)


def test_build_nondimensional_scales():
    scales = build_nondimensional_scales(_FREQUENCY, _MEAN_RADIUS, _BULK_DENSITY)

    second2 = 1.0 / (pi * G * _BULK_DENSITY)
    assert isclose(scales.second2_conversion, second2)
    assert isclose(scales.second_conversion, second2 ** 0.5)
    assert isclose(scales.length_conversion, _MEAN_RADIUS)
    assert isclose(scales.length3_conversion, _MEAN_RADIUS ** 3)
    assert isclose(scales.density_conversion, _BULK_DENSITY)
    assert isclose(scales.mass_conversion, _BULK_DENSITY * _MEAN_RADIUS ** 3)
    assert isclose(scales.pascal_conversion,
                   (_BULK_DENSITY * _MEAN_RADIUS ** 3) / (_MEAN_RADIUS * second2))


def test_matches_classic_implementation():
    """The ported scales equal the classic implementation's while both backends exist."""
    from TidalPy.utilities.dimensions.nondimensional import (
        build_nondimensional_scales as build_classic)

    ported = build_nondimensional_scales(_FREQUENCY, _MEAN_RADIUS, _BULK_DENSITY)
    classic = build_classic(_FREQUENCY, _MEAN_RADIUS, _BULK_DENSITY)
    for name in ("second2_conversion", "second_conversion", "length_conversion",
                 "length3_conversion", "density_conversion", "mass_conversion",
                 "pascal_conversion"):
        assert isclose(getattr(ported, name), getattr(classic, name), rel_tol=1e-15), name
