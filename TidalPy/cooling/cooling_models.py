from typing import Tuple

import numpy as np

from ..utilities.performance.numba import njit
from ..utilities.types import FloatArray, float_eps

MIN_VISCOSITY = 1.
MIN_THICKNESS = 50.
CoolingOutputType = Tuple[FloatArray, FloatArray, FloatArray, FloatArray]


@njit(cacheable=True)
def off(delta_temp: FloatArray, layer_thickness: float) -> CoolingOutputType:
    """ No cooling - Max boundary layer thickness

    !TPY_args live: self.thickness

    Parameters
    ----------
    delta_temp : float
        Difference in temperature between top and mid (or bottom) of layer [K]
        Not used in this model.
    layer_thickness : float
        Thickness of layer [m]

    Returns
    -------
    cooling_flux : FloatArray
        Heat flux leaving the layer [W m-2]
    boundary_layer_thickness : FloatArray
        Thickness of boundary layer (if any) [m]
    rayleigh : FloatArray
        Rayleigh number (only non-zero for convection)
    nusselt : FloatArray
        Nusselt number (only != 1 for convection)
    """

    shape = 0. * delta_temp

    boundary_layer_thickness = (0.5 * layer_thickness) + shape
    cooling_flux = 0. + shape

    # Values to mimic convection that is turned off
    rayleigh = 0. + shape
    nusselt = 1. + shape

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt

@njit(cacheable=True)
def convection(delta_temp: FloatArray,
               viscosity: FloatArray, thermal_conductivity: float, thermal_diffusivity: float, thermal_expansion: float,
               layer_thickness: float, gravity: float, density: float,
               convection_alpha: float, convection_beta: float, critical_rayleigh: float) -> CoolingOutputType:
    """ Calculates cooling by a parameterized convection model via the Rayleigh number

    !TPY_args live: self.viscosity, self.thermal_conductivity, self.thermal_diffusivity, self.thermal_expansion, self.thickness, self.gravity, self.density_bulk
    !TPY_args const: convection_alpha, convection_beta, critical_rayleigh

    Parameters
    ----------
    delta_temp : FloatArray
        Difference in temperature between top and mid (or bottom) of layer [K]
    viscosity : FloatArray
        Viscosity of convecting layer [Pa s]
    thermal_conductivity : float
        Thermal conductivity of layer material [W K-1 m-1]
    thermal_diffusivity : float
        Thermal diffusivity of layer [m2 s-1]
    thermal_expansion : float
        Thermal expansion of layer [K-1]
    layer_thickness : float
        Thickness of layer [m]
    gravity : float
        Surface gravity of layer [m s-2]
    density : float
        Bulk density of layer [kg m-3]
    convection_alpha : float
        Convection scaling term that is ~1. Nusselt number = alpha * (rayleigh / rayleigh_crit) ^ beta
    convection_beta : float
        Nusselt number = alpha * (rayleigh / rayleigh_crit) ^ beta
    critical_rayleigh : float
        Nusselt number = alpha * (rayleigh / rayleigh_crit) ^ beta

    Returns
    -------
    cooling_flux : FloatArray
        Heat flux leaving the layer [W m-2]
    boundary_layer_thickness : FloatArray
        Thickness of boundary layer (if any) [m]
    rayleigh : FloatArray
        Rayleigh number (only non-zero for convection)
    nusselt : FloatArray
        Nusselt number (only != 1 for convection)

    """

    shape = 0. * (delta_temp + viscosity)

    delta_temp_shape = delta_temp + shape
    layer_thickness_shape = layer_thickness + shape

    # Perform main calculations
    rate_heat_loss = thermal_diffusivity / layer_thickness
    parcel_rise_rate = thermal_expansion * density * gravity * delta_temp * layer_thickness**2 / viscosity

    # Calculate the Rayleigh number and then check for over/under shoots
    rayleigh = parcel_rise_rate / rate_heat_loss
    rayleigh = (delta_temp_shape > float_eps) * rayleigh
    rayleigh = (layer_thickness_shape >= MIN_THICKNESS) * rayleigh

    # Calculate the Nusselt number and then check for over/under shoots (minimum nusselt = 2.)
    nusselt = convection_alpha * (rayleigh / critical_rayleigh)**convection_beta
    nusselt = (delta_temp_shape > float_eps) * nusselt + \
              (delta_temp_shape <= float_eps) * 2.
    nusselt = (layer_thickness_shape > MIN_THICKNESS) * nusselt + \
              (layer_thickness_shape <= MIN_THICKNESS) * 2.
    nusselt = (nusselt > 2.) * nusselt + \
              (nusselt <= 2.) * 2.

    # Calculate the thermal boundary layer thickness and then check for over/under shoots
    boundary_layer_thickness = layer_thickness / nusselt
    boundary_layer_thickness = (delta_temp_shape > float_eps) * boundary_layer_thickness + \
                               (delta_temp_shape <= float_eps) * 1.
    boundary_layer_thickness = (layer_thickness_shape > MIN_THICKNESS) * boundary_layer_thickness + \
                               (layer_thickness_shape <= MIN_THICKNESS) * layer_thickness

    # Calculate the parameterized heat flux through the convective layer
    cooling_flux = thermal_conductivity * delta_temp / boundary_layer_thickness

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt

@njit(cacheable=True)
def conduction(delta_temp: FloatArray,
               thermal_conductivity: float, layer_thickness: float) -> CoolingOutputType:
    """ Calculates cooling by conduction through a sub-layer half the thickness of the layer - NonArrays Only

    Half layer thickness is used based on the assumption that delta_temp is the average (central) temperature minus
        the surface temperature.

    !TPY_args live: self.thermal_conductivity, self.thickness

    Parameters
    ----------
    delta_temp : FloatArray
        Difference in temperature between top and mid (or bottom) of layer [K]
    thermal_conductivity : float
        Thermal conductivity of layer material [W K-1 m-1]
    layer_thickness : float
        Thickness of layer [m]

    Returns
    -------
    cooling_flux : FloatArray
        Heat flux leaving the layer [W m-2]
    boundary_layer_thickness : FloatArray
        Thickness of boundary layer (if any) [m]
    rayleigh : FloatArray
        Rayleigh number (only non-zero for convection)
    nusselt : FloatArray
        Nusselt number (only != 1 for convection)

    """

    shape = 0. * delta_temp

    boundary_layer_thickness = layer_thickness + shape

    cooling_flux = thermal_conductivity * delta_temp / boundary_layer_thickness

    # Values to mimic convection that is turned off
    rayleigh = 0. + shape
    nusselt = 1. + shape

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt

# Put New Models Below Here!
