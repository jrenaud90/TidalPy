from typing import Tuple

import numpy as np

from ..performance import njit


MIN_THICKNESS = 0.1


@njit
def off(delta_temp: np.ndarray, layer_thickness: float, thermal_conductivity: float) -> \
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """ Sets cooling to zero


    :return: Cooling flux, Boundary thickness, Rayleigh, Nusselt
    """

    boundary_layer_thickness = (layer_thickness / 2.0) * np.ones_like(delta_temp)
    cooling_flux = np.zeros_like(delta_temp)

    # Values to mimic convection that is turned off
    rayleigh = np.zeros_like(delta_temp)
    nusselt = np.zeros_like(delta_temp)

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt


@njit
def convection(delta_temp: np.ndarray, layer_thickness: float, thermal_conductivity: float,
               viscosity: np.ndarray, thermal_diffusivity: float, thermal_expansion: float,
               gravity: float, density: float,
               convection_alpha: float, convection_beta: float, critical_rayleigh: float) -> \
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """ Calculates cooling by a parameterized convection model via the Rayleigh number

    !TPY_args live: self.layer.viscosity, self.layer.thermal_diffusivity, self.layer.thermal_expansion
    !TPY_args const: gravity, density, convection_alpha, convection_beta, critical_rayleigh

    :return: Cooling flux, Boundary thickness, Rayleigh, Nusselt
    """

    rate_heat_loss = thermal_diffusivity / layer_thickness
    parcel_rise_rate = thermal_expansion * density * gravity * delta_temp * layer_thickness**2 / viscosity

    rayleigh = parcel_rise_rate / rate_heat_loss
    nusselt = convection_alpha * (rayleigh / critical_rayleigh)**convection_beta
    nusselt[nusselt < 1.] = 1.

    # Thickness of boundary layer must be between MIN_THICKNESS and 50% of layer_thickness
    boundary_layer_thickness = layer_thickness / nusselt
    # boundary_layer_thickness[boundary_layer_thickness > .5 * layer_thickness] = .5 * layer_thickness
    boundary_layer_thickness[boundary_layer_thickness < MIN_THICKNESS] = MIN_THICKNESS
    # The minimum thickness sets a (not expressed) limit on the affect of the Nusselt number:
    #     Nusselt_Max = layer_thickness / (2. * MIN_THICKNESS)

    cooling_flux = thermal_conductivity * delta_temp / boundary_layer_thickness

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt


@njit
def conduction(delta_temp: np.ndarray, layer_thickness: float, thermal_conductivity: float) -> \
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """ Calculates cooling by conduction through half of the layer thickness


    :return: Cooling flux, Boundary thickness, Rayleigh, Nusselt
    """

    boundary_layer_thickness = (layer_thickness / 2.0) * np.ones_like(delta_temp)

    cooling_flux = thermal_conductivity * delta_temp / boundary_layer_thickness

    # Values to mimic convection that is turned off
    rayleigh = np.zeros_like(delta_temp)
    nusselt = np.zeros_like(delta_temp)

    return cooling_flux, boundary_layer_thickness, rayleigh, nusselt

# Put New Models Below Here!
