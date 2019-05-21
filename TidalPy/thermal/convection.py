from numba import njit
from ..types import FloatArray
import numpy as np
from typing import Tuple

# Following are using for over and undershoots
MIN_THICKNESS = 0.1

@njit
def calc_convection(viscosity: np.ndarray, delta_temp: np.ndarray,
                    layer_thickness: float, gravity: float, density: float,
                    thermal_diffusivity: float, thermal_expansion: float,
                    convection_alpha: float, convection_beta: float, critical_rayleigh: float) -> \
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Calculates the parameterized convection via the Rayleigh number

    :return: Boundary thickness, Rayleigh, Nusselt
    """

    rate_heat_loss = thermal_diffusivity / layer_thickness
    parcel_rise_rate = thermal_expansion * density * gravity * delta_temp * layer_thickness**2 / viscosity

    rayleigh = parcel_rise_rate / rate_heat_loss
    nusselt = convection_alpha * (rayleigh / critical_rayleigh)**convection_beta
    nusselt[nusselt < 1.] = 1.

    # Thickness of boundary layer must be between MIN_THICKNESS and 50% of layer_thickness
    boundary_layer_thickness = layer_thickness / (2. * nusselt)
    boundary_layer_thickness[boundary_layer_thickness < MIN_THICKNESS] = MIN_THICKNESS
    # The minimum thickness sets a (not expressed) limit on the affect of the Nusselt number:
    #     Nusselt_Max = layer_thickness / (2. * MIN_THICKNESS)

    return boundary_layer_thickness, rayleigh, nusselt


