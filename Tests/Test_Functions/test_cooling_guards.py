"""The convection model's thin-layer guards agree at the minimum layer thickness.

A layer exactly at `MIN_THICKNESS` once kept its Rayleigh number (a `>=` guard) while its Nusselt number and boundary
layer thickness were set to their thin-layer values (`>` guards). It is now too thin for every output.
"""
import pytest

from TidalPy.constants import MIN_THICKNESS
from TidalPy.cooling.cooling_models import convection

# Convecting silicate mantle inputs [K, Pa s, W m-1 K-1, m2 s-1, K-1, m s-2, kg m-3].
DELTA_TEMP = 1000.0
THERMAL_CONDUCTIVITY = 4.0
THERMAL_DIFFUSIVITY = 1.0e-6
THERMAL_EXPANSION = 3.0e-5
GRAVITY = 9.8
DENSITY = 3300.0


def _convection(layer_thickness, viscosity):
    return convection(
        DELTA_TEMP,
        viscosity,
        THERMAL_CONDUCTIVITY,
        THERMAL_DIFFUSIVITY,
        THERMAL_EXPANSION,
        layer_thickness,
        GRAVITY,
        DENSITY,
        1.0,
        1.0 / 3.0,
        1100.0)


def test_guards_agree_at_the_minimum_thickness():
    # A viscosity low enough that the layer would convect hard if its thickness were accepted.
    cooling_flux, boundary_layer_thickness, rayleigh, nusselt = _convection(MIN_THICKNESS, 1.0e-6)
    assert rayleigh == 0.0
    assert nusselt == 2.0
    assert boundary_layer_thickness == MIN_THICKNESS


def test_layer_above_the_minimum_thickness_convects():
    cooling_flux, boundary_layer_thickness, rayleigh, nusselt = _convection(2.0 * MIN_THICKNESS, 1.0e-6)
    assert rayleigh > 0.0
    assert nusselt > 2.0
    assert boundary_layer_thickness == pytest.approx(2.0 * MIN_THICKNESS / nusselt)
