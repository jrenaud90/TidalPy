import pytest

import numpy as np

from math import isnan, isclose, sqrt

from TidalPy.constants import G, PI_DBL
from TidalPy.utilities.dimensions.nondimensional import NonDimensionalScalesClass, build_nondimensional_scales


frequency     = 1.0e-3
mean_radius   = 1.0e6
bulk_density  = 5500.

def test_non_dimensionalize_structure():
    """ Test that the non-dimensionalize structure initializes correctly. """

    test_struct = NonDimensionalScalesClass()
    
    # All conversions should initialize as nans
    assert isnan(test_struct.second2_conversion)
    assert isnan(test_struct.second_conversion)
    assert isnan(test_struct.length_conversion)
    assert isnan(test_struct.length3_conversion)
    assert isnan(test_struct.density_conversion)
    assert isnan(test_struct.mass_conversion)
    assert isnan(test_struct.pascal_conversion)


def test_build_nondimensional_scales():
    """ Test building a non-dimensionalize structure with real inputs. """

    # Build
    test_struct = build_nondimensional_scales(frequency, mean_radius, bulk_density)

    # Check values were set correctly
    assert not isnan(test_struct.second2_conversion)
    assert not isnan(test_struct.second_conversion)
    assert not isnan(test_struct.length_conversion)
    assert not isnan(test_struct.length3_conversion)
    assert not isnan(test_struct.density_conversion)
    assert not isnan(test_struct.mass_conversion)
    assert not isnan(test_struct.pascal_conversion)

    # Check that they are the correct values.
    second2_conv = 1. / (PI_DBL * G * bulk_density)
    assert isclose(test_struct.second2_conversion, second2_conv)
    assert isclose(test_struct.second_conversion, np.sqrt(second2_conv))
    assert isclose(test_struct.length_conversion, mean_radius)
    assert isclose(test_struct.length3_conversion, mean_radius**3)
    assert isclose(test_struct.density_conversion, bulk_density)
    assert isclose(test_struct.mass_conversion, bulk_density * mean_radius**3)
    assert isclose(test_struct.pascal_conversion, (bulk_density * mean_radius**3) / (mean_radius * second2_conv))
