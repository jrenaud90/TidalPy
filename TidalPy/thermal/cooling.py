from numba import njit

from TidalPy.exceptions import ParameterMissingError
from TidalPy.utilities.classes import ModelHolder
from ..types import FloatArray
import numpy as np
from typing import Tuple

# Following are using for over and undershoots
MIN_THICKNESS = 0.1

cooling_param_defaults = {
    'ice': {
        'use_convection': True,
        'convection_alpha': 1.,
        'convection_beta': 1./3.,
        'thermal_conductivity': 2.3,
        'thermal_diffusivity': 2.3 / (1000. * 2000.),
        'thermal_expansion': 5.0e-5,
        'critical_rayleigh': 1600.
    },
    'rock': {
        'use_convection': True,
        'convection_alpha': 1.,
        'convection_beta': 1./3.,
        'thermal_conductivity': 3.75,
        'thermal_diffusivity': 3.75 / (3250. * 1260.),
        'thermal_expansion': 5.2e-5,
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



class Cooling(ModelHolder):

    def __init__(self, layer_type: str, layer_thickness: float, cooling_config: dict = None,
                 gravity: float = None, density: float = None):

        super().__init__(user_config=cooling_config, default_config=cooling_param_defaults[layer_type], automate=False)

        self.layer_thickness = layer_thickness
        self.thermal_conductivity = self.config['thermal_conductivity']
        if self.config['use_convection']:
            self.calculate = self._convection
            if gravity is None:
                raise ParameterMissingError
            if density is None:
                raise ParameterMissingError
            self.inputs = (gravity, density, self.config['thermal_diffusivity'], self.config['thermal_expansion'],
                           self.config['convection_alpha'], self.config['convection_beta'],
                           self.config['critical_rayleigh'])
        else:
            self.inputs = None

    def _conduction(self, temperature: np.ndarray, surface_temperature: FloatArray, layer_thickness: float = None):
        """ Calculates the cooling flux in a layer via conduction

        :param temperature:         <ndarray> Temperature within the layer
        :param surface_temperature: <FloatArray> Temperature at the top of the layer
        :param layer_thickness:     <float> (Optional) Total thickness of the layer in question. Provided in __init__ so it
                                        should only be provided here if the layer is actively growing/shrinking
        :return:                    <Tuple[ndarray's]> Heat Flux, Boundary Thickness, Rayleigh, Nusselt
        """

        delta_temp = temperature - surface_temperature
        if layer_thickness is None:
            cooling_thickness = self.layer_thickness
        else:
            cooling_thickness = layer_thickness

        cooling_flux = self.thermal_conductivity * delta_temp / cooling_thickness

        # Numbers that indicate that convection is off
        rayleigh = np.zeros_like(temperature)
        nusselt = np.ones_like(temperature)

        return cooling_flux, cooling_thickness, rayleigh, nusselt

    def _convection(self, temperature: np.ndarray, surface_temperature: FloatArray, viscosity: np.ndarray,
                    layer_thickness: float = None):
        """ Calculates the cooling flux in a layer via convection

        :param temperature:         <ndarray> Temperature within the layer
        :param surface_temperature: <FloatArray> Temperature at the top of the layer
        :param viscosity:           <ndarray> Viscosity of the layer
        :param layer_thickness:     <float> (Optional) Total thickness of the layer in question. Provided in __init__ so it
                                        should only be provided here if the layer is actively growing/shrinking
        :return:                    <Tuple[ndarray's]> Heat Flux, Boundary Thickness, Rayleigh, Nusselt
        """

        if layer_thickness is None:
            layer_thickness = self.layer_thickness

        delta_temp = temperature - surface_temperature

        # Calculate Convection
        # Wherever delta_temp is negative convection is not happening. To avoid numerical errors we will zero out those
        delta_temp_for_calc = np.copy(delta_temp)
        delta_temp_for_calc[delta_temp_for_calc < 0.] = 0.

        cooling_thickness, rayleigh, nusselt = calc_convection(viscosity, delta_temp_for_calc,
                                                               layer_thickness, *self.inputs)

        # The cooling flux should use the real delta_temp as it could be negative (conducting only)
        cooling_flux = self.thermal_conductivity * delta_temp / cooling_thickness

        # Numbers that indicate that convection is off
        rayleigh = np.zeros_like(temperature)
        nusselt = np.ones_like(temperature)

        return cooling_flux, cooling_thickness, rayleigh, nusselt