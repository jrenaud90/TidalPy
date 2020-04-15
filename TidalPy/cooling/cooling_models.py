from typing import Tuple

import numpy as np

from ..performance import njit
from ..types import float_min


MIN_THICKNESS = 50.
CoolingOutputType = Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]

@njit
def off(delta_temp: np.ndarray, layer_thickness: float) -> CoolingOutputType:
    """ No cooling - Max boundary layer thickness

    !TPY_args live: self.thickness

    Parameters
    ----------
    delta_temp : np.ndarray
        Difference in temperature between top and mid (or bottom) of layer [K]
    layer_thickness : float
        Thickness of layer [m]

    Returns
    -------
    cooling_flux : np.ndarray
        Heat flux leaving the layer [W m-2]
    boundary_layer_thickness : np.ndarray
        Thickness of boundary layer (if any) [m]
    rayleigh : np.ndarray
        Rayleigh number (only non-zero for convection)
    nusselt : np.ndarray
        Nusselt number (only != 1 for convection)

    """

    boundary_layer_thickness = (0.5 * layer_thickness) * np.ones_like(delta_temp)
    cooling_flux = np.zeros_like(delta_temp)

    # Values to mimic convection that is turned off
    rayleigh = np.zeros_like(cooling_flux)
    nusselt = np.zeros_like(cooling_flux)

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt


@njit
def convection(delta_temp: np.ndarray,
               viscosity: np.ndarray, thermal_conductivity: float, thermal_diffusivity: float, thermal_expansion: float,
               layer_thickness: float, gravity: float, density: float,
               convection_alpha: float, convection_beta: float, critical_rayleigh: float) -> CoolingOutputType:
    """ Calculates cooling by a parameterized convection model via the Rayleigh number

    !TPY_args live: self.viscosity, self.thermal_conductivity, self.thermal_diffusivity, self.thermal_expansion, self.thickness, self.gravity, self.density_bulk
    !TPY_args const: gravity, density, convection_alpha, convection_beta, critical_rayleigh

    Parameters
    ----------
    delta_temp : np.ndarray
        Difference in temperature between top and mid (or bottom) of layer [K]
    viscosity : np.ndarray
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
    cooling_flux : np.ndarray
        Heat flux leaving the layer [W m-2]
    boundary_layer_thickness : np.ndarray
        Thickness of boundary layer (if any) [m]
    rayleigh : np.ndarray
        Rayleigh number (only non-zero for convection)
    nusselt : np.ndarray
        Nusselt number (only != 1 for convection)

    """

    rate_heat_loss = thermal_diffusivity / layer_thickness
    delta_temp[delta_temp < 0.] = 0.
    viscosity[viscosity < float_min] = float_min
    
    # TODO: Right now negative delta-temps do not lead to backwards cooling.
    #  Convection should def shut off, but there could still be backwards conduction. How to handle this and keep
    #  the boundary layer still smaller (doesn't make sense for it to explode to 50% layer thickness between two time
    #  steps).
    parcel_rise_rate = thermal_expansion * density * gravity * delta_temp * layer_thickness**2 / viscosity
    rayleigh = parcel_rise_rate / rate_heat_loss
    nusselt = convection_alpha * (rayleigh / critical_rayleigh)**convection_beta
    # Technically convection should shut off at nusselt == 1, but the 2 here forces the max boundary layer thickness
    #  to be 50% the layer thickness.
    nusselt[nusselt < 2.] = 2.

    # Thickness of boundary layer must be between MIN_THICKNESS and 50% of layer_thickness
    boundary_layer_thickness = layer_thickness / nusselt
    boundary_layer_thickness[boundary_layer_thickness < MIN_THICKNESS] = MIN_THICKNESS

    cooling_flux = thermal_conductivity * delta_temp / boundary_layer_thickness

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt


@njit
def conduction(delta_temp: np.ndarray,
               thermal_conductivity: float, layer_thickness: float) -> CoolingOutputType:
    """ Calculates cooling by conduction through a sub-layer half the thickness of the layer.

    Half layer thickness is used based on the assumption that delta_temp is the average (central) temperature minus
        the surface temperature.

    !TPY_args live: self.thermal_conductivity, self.thickness
    !TPY_args const: gravity, density, convection_alpha, convection_beta, critical_rayleigh

    Parameters
    ----------
    delta_temp : np.ndarray
        Difference in temperature between top and mid (or bottom) of layer [K]
    thermal_conductivity : float
        Thermal conductivity of layer material [W K-1 m-1]
    layer_thickness : float
        Thickness of layer [m]

    Returns
    -------
    cooling_flux : np.ndarray
        Heat flux leaving the layer [W m-2]
    boundary_layer_thickness : np.ndarray
        Thickness of boundary layer (if any) [m]
    rayleigh : np.ndarray
        Rayleigh number (only non-zero for convection)
    nusselt : np.ndarray
        Nusselt number (only != 1 for convection)

    """

    boundary_layer_thickness = (0.5 * layer_thickness) * np.ones_like(delta_temp)

    cooling_flux = thermal_conductivity * delta_temp / boundary_layer_thickness

    # Values to mimic convection that is turned off
    rayleigh = np.zeros_like(cooling_flux)
    nusselt = np.zeros_like(cooling_flux)

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt

# Put New Models Below Here!