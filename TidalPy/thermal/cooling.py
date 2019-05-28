from numba import njit

from TidalPy.exceptions import ParameterMissingError
from TidalPy.utilities.classes import ModelHolder, LayerModel
from ..types import FloatArray
import numpy as np
from ..structures.layers import LayerType
from typing import Tuple

# Following are using for over and undershoots
MIN_THICKNESS = 0.1

cooling_param_defaults = {
    'ice': {
        'use_convection': True,
        'convection_alpha': 1.,
        'convection_beta': 1./3.,
        'critical_rayleigh': 1600.
    },
    'rock': {
        'use_convection': True,
        'convection_alpha': 1.,
        'convection_beta': 1./3.,
        'critical_rayleigh': 1100.
    },
    'iron': {
        'use_convection': False,
    }
}


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



class Cooling(LayerModel):

    default_config = cooling_param_defaults
    config_key = 'cooling'

    def __init__(self, layer: LayerType):

        super().__init__(layer=layer, function_searcher=None, automate=True)

        if self.config['use_convection']:
            self.calculate = self._convection
            self.inputs = (self.config['convection_alpha'], self.config['convection_beta'],
                           self.config['critical_rayleigh'])
        else:
            self.calculate = self._conduction
            self.inputs = None

    def _conduction(self):
        """ Calculates the cooling flux in a layer via conduction

        :return: <Tuple[ndarray's]> Heat Flux, Boundary Thickness, Rayleigh, Nusselt
        """

        temperature = self.layer.temperature
        delta_temp = temperature - self.layer.temperature_surf
        cooling_thickness = self.layer.thickness

        cooling_flux = self.layer.thermal_conductivity * delta_temp / cooling_thickness

        # Numbers that indicate that convection is off
        rayleigh = np.zeros_like(temperature)
        nusselt = np.ones_like(temperature)

        return cooling_flux, cooling_thickness, rayleigh, nusselt

    def _convection(self):
        """ Calculates the cooling flux in a layer via convection

        :return: <Tuple[ndarray's]> Heat Flux, Boundary Thickness, Rayleigh, Nusselt
        """

        temperature = self.layer.temperature
        delta_temp = temperature - self.layer.temperature_surf
        cooling_thickness = self.layer.thickness

        # Calculate Convection
        # Wherever delta_temp is negative convection is not happening. To avoid numerical errors we will zero out those
        delta_temp_for_calc = np.copy(delta_temp)
        delta_temp_for_calc[delta_temp_for_calc < 0.] = 0.

        cooling_thickness, rayleigh, nusselt = \
            calc_convection(self.layer.viscosity, delta_temp_for_calc, cooling_thickness, self.layer.gravity,
                            self.layer.density, self.layer.thermal_diffusivity, self.layer.thermal_expansion,
                            *self.inputs)

        # The cooling flux should use the real delta_temp as it could be negative (conducting only)
        cooling_flux = self.layer.thermal_conductivity * delta_temp / cooling_thickness

        # Numbers that indicate that convection is off
        rayleigh = np.zeros_like(temperature)
        nusselt = np.ones_like(temperature)

        return cooling_flux, cooling_thickness, rayleigh, nusselt