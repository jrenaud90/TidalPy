"""The convection model's Nusselt convention: Nu = 1 is conduction, so a sub-critical or rigid layer conducts."""
import numpy as np
import pytest

import TidalPy
from TidalPy import Cooling

# A 100 km silicate layer with a 500 K drop, MKS.
_DELTA_TEMP = 500.0
_THICKNESS = 1.0e5
_GRAVITY = 3.0
_DENSITY = 3500.0
_CONDUCTIVITY = 3.75
_DIFFUSIVITY = _CONDUCTIVITY / (3500.0 * 1200.0)
_EXPANSION = 5.2e-5
_CRITICAL_RAYLEIGH = 1100.0


def _convective(delta_temp, viscosity):
    return Cooling.convective(
        delta_temp, _THICKNESS, _GRAVITY, _DENSITY, viscosity, _CONDUCTIVITY, _DIFFUSIVITY, _EXPANSION)


def _rayleigh(viscosity):
    return _EXPANSION * _DENSITY * _GRAVITY * _DELTA_TEMP * _THICKNESS ** 3 / (viscosity * _DIFFUSIVITY)


def test_the_default_nusselt_floor_is_conduction():
    assert TidalPy.config["numerical"]["minimum_nusselt"] == 1.0


@pytest.mark.parametrize("viscosity", [1.0e22, 1.0e26, np.inf], ids=["subcritical", "stiff", "rigid"])
def test_a_layer_that_does_not_convect_reports_the_conductive_flux(viscosity):
    """Below the critical Rayleigh number the convection model gives what the conduction model gives."""
    assert _rayleigh(viscosity) < _CRITICAL_RAYLEIGH
    convecting = _convective(_DELTA_TEMP, viscosity)
    conducting = Cooling.conductive(_DELTA_TEMP, _THICKNESS, _CONDUCTIVITY)
    assert convecting.nusselt == 1.0
    assert convecting.cooling_flux == pytest.approx(conducting.cooling_flux, rel=1.0e-12)
    assert convecting.boundary_layer_thickness == pytest.approx(conducting.boundary_layer_thickness, rel=1.0e-12)


def test_a_convecting_layer_carries_nusselt_times_the_conductive_flux():
    """q = Nu k dT / D, with Nu the flux over that of conduction across the whole layer, and a boundary layer D / Nu."""
    viscosity = 1.0e18
    result = _convective(_DELTA_TEMP, viscosity)
    expected_nusselt = (_rayleigh(viscosity) / _CRITICAL_RAYLEIGH) ** (1.0 / 3.0)
    assert result.nusselt == pytest.approx(expected_nusselt, rel=1.0e-12)
    conductive_flux = _CONDUCTIVITY * _DELTA_TEMP / _THICKNESS
    assert result.cooling_flux == pytest.approx(result.nusselt * conductive_flux, rel=1.0e-12)
    assert result.boundary_layer_thickness == pytest.approx(_THICKNESS / result.nusselt, rel=1.0e-12)


def test_the_flux_is_continuous_at_the_onset_of_convection():
    """Just above the critical Rayleigh number the flux leaves the conductive value without a jump."""
    onset_viscosity = _rayleigh(1.0) / _CRITICAL_RAYLEIGH
    below = _convective(_DELTA_TEMP, onset_viscosity * (1.0 + 1.0e-6))
    above = _convective(_DELTA_TEMP, onset_viscosity * (1.0 - 1.0e-6))
    assert below.nusselt == 1.0
    assert above.nusselt > 1.0
    assert above.cooling_flux == pytest.approx(below.cooling_flux, rel=1.0e-5)


def test_a_degenerate_layer_sits_on_the_conductive_floor():
    """No contrast gives no flux and the whole layer as its boundary layer."""
    result = _convective(0.0, 1.0e18)
    assert result.cooling_flux == 0.0
    assert result.rayleigh == 0.0
    assert result.nusselt == 1.0
    assert result.boundary_layer_thickness == pytest.approx(_THICKNESS)
