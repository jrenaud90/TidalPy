from typing import Tuple

import numpy as np

from ..performance import njit
from ..types import FloatArray


MIN_THICKNESS = 0.1


@njit
def off(delta_temp: FloatArray, layer_thickness: float, thermal_conductivity: float) -> \
        Tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """ Sets cooling to zero

    Parameters
    ----------
    delta_temp : FloatArray
        Difference in temperature between top and mid (or bottom) of layer [K]
    layer_thickness : float
        Thickness of layer [m]
    thermal_conductivity : float
        Thermal conductivity of layer material [W K-1 m-1]

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

    boundary_layer_thickness = (layer_thickness / 2.0) * np.ones_like(delta_temp)
    cooling_flux = np.zeros_like(delta_temp)

    # Values to mimic convection that is turned off
    rayleigh = np.zeros_like(delta_temp)
    nusselt = np.zeros_like(delta_temp)

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt


@njit
def convection(delta_temp: FloatArray, layer_thickness: float, thermal_conductivity: float,
               viscosity: FloatArray, thermal_diffusivity: float, thermal_expansion: float,
               gravity: float, density: float,
               convection_alpha: float, convection_beta: float, critical_rayleigh: float) -> \
        Tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """ Calculates cooling by a parameterized convection model via the Rayleigh number

    !TPY_args live: self.layer.viscosity, self.layer.thermal_diffusivity, self.layer.thermal_expansion
    !TPY_args const: gravity, density, convection_alpha, convection_beta, critical_rayleigh

    Parameters
    ----------
    delta_temp : FloatArray
        Difference in temperature between top and mid (or bottom) of layer [K]
    layer_thickness : float
        Thickness of layer [m]
    thermal_conductivity : float
        Thermal conductivity of layer material [W K-1 m-1]
    viscosity : FloatArray
        Viscosity of convecting layer [Pa s]
    thermal_diffusivity : float
        Thermal diffusivity of layer [m2 s-1]
    thermal_expansion : float
        Thermal expansion of layer [K-1]
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

    rate_heat_loss = thermal_diffusivity / layer_thickness
    parcel_rise_rate = thermal_expansion * density * gravity * delta_temp * layer_thickness**2 / viscosity

    rayleigh = parcel_rise_rate / rate_heat_loss
    nusselt = convection_alpha * (rayleigh / critical_rayleigh)**convection_beta
    nusselt[nusselt < 1.] = 1.

    # Thickness of boundary layer must be between MIN_THICKNESS and 50% of layer_thickness
    boundary_layer_thickness = layer_thickness / nusselt
    # TODO: look at the below commented out equation - some people use this some dont.
    # boundary_layer_thickness[boundary_layer_thickness > .5 * layer_thickness] = .5 * layer_thickness
    boundary_layer_thickness[boundary_layer_thickness < MIN_THICKNESS] = MIN_THICKNESS
    # The minimum thickness sets a (not expressed) limit on the affect of the Nusselt number:
    #     Nusselt_Max = layer_thickness / (2. * MIN_THICKNESS)

    cooling_flux = thermal_conductivity * delta_temp / boundary_layer_thickness

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt


@njit
def conduction(delta_temp: FloatArray, layer_thickness: float, thermal_conductivity: float) -> \
        Tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """ Calculates cooling by conduction through the full thickness of the layer

    Parameters
    ----------
    delta_temp : FloatArray
        Difference in temperature between top and mid (or bottom) of layer [K]
    layer_thickness : float
        Thickness of layer [m]
    thermal_conductivity : float
        Thermal conductivity of layer material [W K-1 m-1]

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
    """ Calculates cooling by conduction through half of the layer thickness


    :return: Cooling flux, Boundary thickness, Rayleigh, Nusselt
    """

    # TODO: Should this be the full thickness or half the thickness?
    boundary_layer_thickness = layer_thickness * np.ones_like(delta_temp)

    cooling_flux = thermal_conductivity * delta_temp / boundary_layer_thickness

    # Values to mimic convection that is turned off
    rayleigh = np.zeros_like(delta_temp)
    nusselt = np.zeros_like(delta_temp)

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt

# Put New Models Below Here!
