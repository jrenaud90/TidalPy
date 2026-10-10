"""Tests for ``TidalPy.Utilities.dimensions.nondimensional`` (non-dimensionalization scales).

Classic scales were frozen from the 0.8.0 snapshot 8b8e0b12 into ``frozen/test_nondimensional_01.npz``.
"""
from math import isclose, isnan
from pathlib import Path

import numpy as np
import pytest

from TidalPy.constants import G, pi
from TidalPy.Utilities.dimensions import NonDimensionalScalesClass, build_nondimensional_scales

_FREQUENCY = 1.0e-3
_MEAN_RADIUS = 1.0e6
_BULK_DENSITY = 5500.0
_SECOND2 = 1.0 / (pi * G * _BULK_DENSITY)

# Classic build_nondimensional_scales(_FREQUENCY, _MEAN_RADIUS, _BULK_DENSITY), one scalar per scale name.
_FROZEN_PATH = Path(__file__).parent / "frozen" / "test_nondimensional_01.npz"
with np.load(_FROZEN_PATH, allow_pickle=False) as _frozen_file:
    _CLASSIC_SCALES = {key: float(_frozen_file[key]) for key in _frozen_file.files}

_EXPECTED_SCALES = {
    "second2_conversion": _SECOND2,
    "second_conversion": _SECOND2 ** 0.5,
    "length_conversion": _MEAN_RADIUS,
    "length3_conversion": _MEAN_RADIUS ** 3,
    "density_conversion": _BULK_DENSITY,
    "mass_conversion": _BULK_DENSITY * _MEAN_RADIUS ** 3,
    "pascal_conversion": (_BULK_DENSITY * _MEAN_RADIUS ** 3) / (_MEAN_RADIUS * _SECOND2),
}


@pytest.mark.parametrize("name", list(_EXPECTED_SCALES))
def test_non_dimensionalize_structure_initializes_nan(name):
    """A default constructed scales structure holds NaN."""
    assert isnan(getattr(NonDimensionalScalesClass(), name))


@pytest.mark.parametrize("name", list(_EXPECTED_SCALES))
def test_build_nondimensional_scales(name):
    """Each built scale matches its analytic value and the frozen classic value."""
    scales = build_nondimensional_scales(_MEAN_RADIUS, _BULK_DENSITY)
    assert isclose(getattr(scales, name), _EXPECTED_SCALES[name])
    assert isclose(getattr(scales, name), _CLASSIC_SCALES[name], rel_tol=1e-15)


def test_builder_takes_no_frequency():
    """The builder rejects the classic, unused frequency argument."""
    with pytest.raises(TypeError):
        build_nondimensional_scales(_FREQUENCY, _MEAN_RADIUS, _BULK_DENSITY)
