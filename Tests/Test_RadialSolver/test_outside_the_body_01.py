"""A solution's profile getters read NaN above the surface on both the shooting and the propagation-matrix path."""
import math

import pytest

from TidalPy.RadialSolver.helpers import homogeneous_love_numbers

RADIUS = 1.0e6
DENSITY = 3000.0
SHEAR = 5.0e10 + 1.0e8j
FREQUENCY = 1.0e-5


@pytest.mark.parametrize("method", ["shooting", "pm"])
def test_the_moduli_are_nan_above_the_surface(method):
    # The propagation matrix needs a static, incompressible layer.
    solution = homogeneous_love_numbers(
        RADIUS, DENSITY, SHEAR, FREQUENCY, layer_is_static=True, layer_is_incompressible=(method == "pm"),
        love_method=method)
    inside = solution.get_complex_shear_modulus(0.5 * RADIUS)
    assert inside.real == pytest.approx(SHEAR.real)
    outside = solution.get_complex_shear_modulus(1.5 * RADIUS)
    assert math.isnan(outside.real) and math.isnan(outside.imag)
