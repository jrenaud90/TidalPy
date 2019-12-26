import copy
from typing import TYPE_CHECKING, List, Dict

import numpy as np

from ..constants import G
from . import calculate_tides, calculate_tides_withmodes, calc_tidal_susceptibility, calc_tidal_susceptibility_reduced
from ..exceptions import (AttributeNotSetError, ImplementationException, ParameterMissingError, UnknownModelError,
                          BadValueError, AttributeChangeRequiresReINIT, ImproperAttributeHandling, ParameterValueError,
                          OuterscopeAttributeSetError, ConfigAttributeChangeError, MissingArgumentError)

from ..dynamics import mode_types
from ..utilities.classes import WorldConfigHolder, LayerConfigHolder
from ..types import ArrayNone
from ..utilities.model import LayerModelHolder
from .calculation import calc_tidal_susceptibility
from .love_1d import complex_love_general, effective_rigidity_general

if TYPE_CHECKING:
    from ..structures import TidalWorld
    from ..structures import ThermalLayer


class Tides(WorldConfigHolder):

    """ Tides class - Holder for all tidal heating and torque calculation functions

    Tides class stores model parameters and methods for heating and torque which are functions of
        (T, P, melt_frac, w, e, I)
    """

    world_config_key = 'tides'

    def __init__(self, world: 'TidalWorld', store_config_in_world: bool = True):

        super().__init__(world, store_config_in_world=store_config_in_world)

        # Model setup
        self._orbital_truncation_lvl = self.config['orbital_truncation_level']
        self._tidal_order_lvl = self.config['max_tidal_order_l']

        # Ensure the tidal order and orbital truncation levels make sense
        if self._tidal_order_lvl > 3:
            raise ImplementationException(f'Tidal order {self._tidal_order_lvl} has not been implemented yet.')
        if self._orbital_truncation_lvl % 2 != 0:
            raise ParameterValueError('Orbital truncation level must be an even integer.')
        if self._orbital_truncation_lvl <= 2:
            raise ParameterValueError('Orbital truncation level must be greater than or equal to 2.')
        if self._orbital_truncation_lvl not in [2, 4, 6]:
            raise ImplementationException(f'Orbital truncation level of {self._orbital_truncation_lvl} is not '
                                          f'currently supported.')

        # Setup mode calculator
        self._tidal_order_num_list = range(2, self._tidal_order_lvl+1)
        self.mode_calculators = [mode_types[self.use_nsr][order_l][self._orbital_truncation_lvl]
                                 for order_l in self._tidal_order_num_list]

        # State properties
        self._modes = None
        self._frequencies = None
        self._unique_tidal_frequencies = None
        self._mode_names = None
        self._heating_subterms = None
        self._ztorque_subterms = None
        self._tidal_susceptibility = None
        self._tidal_susceptibility_reduced = None
        self._tidal_heating = None

        # Parameters initialized from configuration but can be changed later
        self.quality_factor = self.config['quality_factor']

    def orbit_update(self):
        """ Should be called whenever the orbit is updated
        """

        self.calculate_tidal_susceptibility()
        self.calculate_modes()


    def calculate_tidal_susceptibility(self):
        """ Calculate the tidal susceptibility.

        Returns
        -------
        tidal_susceptibility : np.ndarray
            Tidal Susceptibility [N m]
        """

        self._tidal_susceptibility = \
            calc_tidal_susceptibility(self.tidal_host.mass, self.world.radius, self.semi_major_axis)

        return self.tidal_susceptibility

    def calculate_tidal_susceptibility_reduced(self):
        """ Calculate the reduced tidal susceptibility.

        Returns
        -------
        tidal_susceptibility_reduced : np.ndarray
            Reduced Tidal Susceptibility [kg m8 s-2]

        """

        self._tidal_susceptibility_reduced = calc_tidal_susceptibility_reduced(self.tidal_host.mass, self.world.radius)

        return self._tidal_susceptibility_reduced

    def calculate_modes(self):
        """ Calculate tidal heating and torques based on the layer's current state.

        Relies on the layer.rheology class' love_numbers

        Returns
        -------
        modes
        tidal_heating : np.ndarray
        tidal_ztorque : np.ndarray
        """

        unique_tidal_freqs = dict()
        mode_names_by_orderl = list()
        modes_by_orderl = list()
        freqs_by_orderl = list()
        heating_subterms_by_orderl = list()
        ztorque_subterms_by_orderl = list()

        for mode_calculator in self.mode_calculators:
            mode_names, modes, freqs, heating_subterms, ztorque_subterms = \
                mode_calculator(self.orbital_freq, self.spin_freq, self.eccentricity, self.inclination)

            mode_names_by_orderl.append(mode_names)
            modes_by_orderl.append(modes)
            freqs_by_orderl.append(freqs)
            heating_subterms_by_orderl.append(heating_subterms)
            ztorque_subterms_by_orderl.append(ztorque_subterms)

            for mode_name, freq in zip(mode_names, freqs):
                if mode_name not in unique_tidal_freqs:
                    unique_tidal_freqs[mode_name] = freq

        self._modes = modes_by_orderl
        self._frequencies = freqs_by_orderl
        self._unique_tidal_frequencies = unique_tidal_freqs
        self._mode_names = mode_names_by_orderl
        self._heating_subterms = heating_subterms_by_orderl
        self._ztorque_subterms = ztorque_subterms_by_orderl

        return self.modes, self.heating_subterms, self.ztorque_subterms

    def calculate_tides_for_layer(self, layer: 'ThermalLayer'):
        """ Calculate heating and torque for a thermal layer.

        Assumes that tidal modes and frequencies have already been calculated and updated within the function.

        Parameters
        ----------
        layer : ThermalLayer
            Layer in which you want to calculate tides in.

        Returns
        -------
        tidal_heating : np.ndarray
            Tidal heating [W]
        tidal_ztorque : np.ndarray
            Tidal ztorque [N m]
        """

        tidal_heating_terms = list()
        tidal_ztorque_terms = list()



    def calculate_from_complex_func(self, complex_compliance_func, layer: 'ThermalLayer' = None,
                                    compliance: np.ndarray = None, viscosity: np.ndarray = None,
                                    complex_compliance_inputs: tuple = None, tidal_scale: float = 1.):

        if complex_compliance_inputs is None:
            complex_compliance_inputs = tuple()

        # First calculate the modes and subterms
        self.calculate_modes()

        # Determine the viscosity and compliances
        gravity, radius, density_bulk = self.world.gravity, self.world.radius, self.world.density_bulk
        if layer is not None:
            viscosity = layer.viscosity
            compliance = layer.compliance

            # If layer is provided use its geometry
            gravity, radius, density_bulk = layer.gravity, layer.radius, layer.density_bulk
            tidal_scale = layer.tidal_scale
        else:
            if compliance is None or viscosity is None:
                raise MissingArgumentError('Viscosity and compliance must be provided if layer is not.')

        shear = 1. / compliance

        tidal_heating, tidal_ztorque, love_numbers, tidal_frequencies, cached_complex_comps = \
                calculate_tides_withmodes(shear, viscosity, gravity, radius, density_bulk, self.tidal_susceptibility,
                                          complex_compliance_func, complex_compliance_inputs, self.semi_major_axis,
                                          self.mode_names, self.modes, self.frequencies, self.heating_subterms,
                                          self.ztorque_subterms, tidal_scale, self.use_nsr,
                                          self.orbital_truncation_lvl, self.tidal_order_lvl)

        return tidal_heating, tidal_ztorque

    # Configuration properties
    @property
    def orbital_truncation_lvl(self) -> int:
        return self._orbital_truncation_lvl

    @orbital_truncation_lvl.setter
    def orbital_truncation_lvl(self, value):
        raise ConfigAttributeChangeError

    @property
    def tidal_order_lvl(self) -> int:
        return self._tidal_order_lvl

    @tidal_order_lvl.setter
    def tidal_order_lvl(self, value):
        raise ConfigAttributeChangeError

    # State properties
    @property
    def modes(self) -> List[List[np.ndarray]]:
        return self._modes

    @modes.setter
    def modes(self, value):
        raise ImproperAttributeHandling

    @property
    def frequencies(self) -> List[List[np.ndarray]]:
        return self._frequencies

    @frequencies.setter
    def frequencies(self, value):
        raise ImproperAttributeHandling

    @property
    def unique_tidal_frequencies(self) -> Dict[str, np.ndarray]:
        return self._unique_tidal_frequencies

    @unique_tidal_frequencies.setter
    def unique_tidal_frequencies(self, value):
        raise ImproperAttributeHandling

    @property
    def mode_names(self) -> List[List[str]]:
        return self._mode_names

    @mode_names.setter
    def mode_names(self, value):
        raise ImproperAttributeHandling

    @property
    def heating_subterms(self) -> List[List[np.ndarray]]:
        return self._heating_subterms

    @heating_subterms.setter
    def heating_subterms(self, value):
        raise ImproperAttributeHandling

    @property
    def ztorque_subterms(self) -> List[List[np.ndarray]]:
        return self._ztorque_subterms

    @ztorque_subterms.setter
    def ztorque_subterms(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_susceptibility(self) -> np.ndarray:
        return self._tidal_susceptibility

    @tidal_susceptibility.setter
    def tidal_susceptibility(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_susceptibility_reduced(self) -> np.ndarray:
        return self._tidal_susceptibility_reduced

    @tidal_susceptibility_reduced.setter
    def tidal_susceptibility_reduced(self, value):
        raise ImproperAttributeHandling

    # Outer scope properties
    @property
    def semi_major_axis(self):
        return self.world.semi_major_axis

    @semi_major_axis.setter
    def semi_major_axis(self, value):
        raise OuterscopeAttributeSetError

    @property
    def tidal_host(self):
        return self.world.tidal_host

    @tidal_host.setter
    def tidal_host(self, value):
        raise OuterscopeAttributeSetError

    @property
    def use_nsr(self):
        return not self.world.is_spin_sync

    @use_nsr.setter
    def use_nsr(self, value):
        raise OuterscopeAttributeSetError

    @property
    def eccentricity(self):
        return self.world.eccentricity

    @eccentricity.setter
    def eccentricity(self, value):
        raise OuterscopeAttributeSetError

    @property
    def inclination(self):
        return self.world.inclination

    @inclination.setter
    def inclination(self, value):
        raise OuterscopeAttributeSetError

    @property
    def spin_freq(self):
        return self.world.spin_freq

    @spin_freq.setter
    def spin_freq(self, value):
        raise OuterscopeAttributeSetError

    @property
    def orbital_freq(self):
        return self.world.orbital_freq

    @orbital_freq.setter
    def orbital_freq(self, value):
        raise OuterscopeAttributeSetError



#
#
#     # No defaults or config information needed in this version of TidalPy.tides
#     #    I am keeping tides a child class of configholder in case there is ever a reason to add in config info - it
#     #    will be easy to do so.
#     default_config = None
#     layer_config_key = None
#
#     def __init__(self, world: TidalWorld, store_config_in_world: bool = True):
#
#         super().__init__(world, store_config_in_world)
#
#         # State properties
#         self._tidal_susceptibility = None
#         self._tidal_susceptibility_reduced = None
#
#         # Pull out configurations
#         self._use_nsr = self.config['use_nsr']
#         self._max_tidal_order_l = self.config['max_tidal_order_l']
#         self._orbital_truncation = self.config['orbital_truncation']
#
#         # Ensure the tidal order and orbital truncation levels make sense
#         if self.max_tidal_order_l > 3:
#             raise ImplementationException(f'Tidal order {self.max_tidal_order_l} has not been implemented yet.')
#         if self.orbital_truncation % 2 != 0:
#             raise ParameterValueError('Orbital truncation level must be an even integer.')
#         if self.orbital_truncation <= 2:
#             raise ParameterValueError('Orbital truncation level must be greater than or equal to 2.')
#         if self.orbital_truncation not in [2, 4, 6]:
#             raise ImplementationException(f'Orbital truncation level of {self.orbital_truncation} is not currently '
#                                           f'supported.')
#         self.order_l_list = range(2, self.max_tidal_order_l + 1)
#
#         # Find the correct mode functions for heating and torque
#         self.tidal_mode_calculators = [mode_types[self.use_nsr][l][self.orbital_truncation] for l in self.order_l_list]
#
#     def calc_tidal_susceptibility_reduced(self) -> float:
#         """ Reduced version of tidal susceptibility that does not take into account the orbital seperation between the
#         target body and host.
#
#         Returns
#         -------
#         tidal_susceptibility_reduced: float
#             Reduced tidal susceptibility [N m]
#         """
#
#         self._tidal_susceptibility_reduced = (3. / 2.) * G * self.host.mass**2 * self.world.radius**5
#
#         return self._tidal_susceptibility_reduced
#
#     def calc_tidal_susceptibility(self) -> np.ndarray:
#         """ Tidal susceptibility is the relative strength of tides based on a planet's size, its host's mass, and their
#             orbital seperation.
#
#         Tidal susceptibility is equal to tidal_susceptibility_reduced / a^6
#
#         Returns
#         -------
#         tidal_susceptibility : np.ndarray
#             Tidal susceptibility of this planet.
#         """
#
#         if self.tidal_susceptibility_reduced is None:
#             self.calc_tidal_susceptibility_reduced()
#
#         self._tidal_susceptibility = self.tidal_susceptibility_reduced / self.semi_major_axis**6
#
#         return self._tidal_susceptibility
#
#     def calc_tidal_modes(self):
#         """ Calculates the tidal modes, frequencies, and heating & torque coefficients based on the planet's spin-rate
#             and orbital motion.
#
#         Returns
#         -------
#
#
#         """
#
#         mode_names = list()
#         modes = list()
#         freqs = list()
#         heating_subterms = list()
#
#         tidal_modes = list()
#         tidal
#
#         mode_names, modes, freqs, heating_subterms, ztorque_subterms
#
#
#     # Outerscope properties
#     @property
#     def host(self):
#         if self.world.host is None:
#             raise AttributeNotSetError
#         return self.world.host
#
#     @host.setter
#     def host(self, value):
#         raise ImproperAttributeHandling
#
#     @property
#     def semi_major_axis(self):
#         if self.world.semi_major_axis is None:
#             raise AttributeNotSetError
#         return self.world.semi_major_axis
#
#     @semi_major_axis.setter
#     def semi_major_axis(self, value):
#         raise ImproperAttributeHandling
#
#     # State properties
#     @property
#     def use_nsr(self):
#         return self._use_nsr
#
#     @use_nsr.setter
#     def use_nsr(self, value):
#         raise AttributeChangeRequiresReINIT
#
#     @property
#     def max_tidal_order_l(self):
#         return self._max_tidal_order_l
#
#     @max_tidal_order_l.setter
#     def max_tidal_order_l(self, value):
#         raise AttributeChangeRequiresReINIT
#
#     @property
#     def orbital_truncation(self):
#         return self._orbital_truncation
#
#     @orbital_truncation.setter
#     def orbital_truncation(self, value) -> float:
#         raise AttributeChangeRequiresReINIT
#
#     @property
#     def tidal_susceptibility_reduced(self):
#         return self._tidal_susceptibility_reduced
#
#     @tidal_susceptibility_reduced.setter
#     def tidal_susceptibility_reduced(self, value):
#         raise ImproperAttributeHandling
#
#     @property
#     def tidal_susceptibility(self) -> np.ndarray:
#         return self._tidal_susceptibility
#
#     @tidal_susceptibility.setter
#     def tidal_susceptibility(self, value):
#         raise ImproperAttributeHandling
#
#
#
#     class Rheology(LayerConfigHolder):
#
# class Tides(LayerModelHolder):
#     default_config = copy.deepcopy(rheology_param_defaults)
#     config_key = 'rheology'
#
#     def __init__(self, layer):
#
#         super().__init__(layer=layer, function_searcher=None, call_reinit=True)
#
#         # Override model if layer has tidal off
#         if not self.layer.config['is_tidal']:
#             self.model = 'off'
#
#         # Pull out model information
#         self.orbital_truncation_level = self.world.config['orbital_truncation_level']
#         self.order_l = self.world.config['tidal_order_l']
#
#         # Pull out information about the planet/layer
#         self.config['planet_beta'] = self.layer.world.beta
#         self.config['quality_factor'] = self.layer.world.fixed_q
#
#         if int(self.order_l) != self.order_l or self.order_l < 2:
#             raise BadValueError('Tidal order l must be an integer >= 2')
#         if self.order_l > 3:
#             raise ImplementationException('Tidal order l > 3 has not been fully implemented in TidalPy.')
#
#         if int(self.orbital_truncation_level) != self.orbital_truncation_level or \
#                 self.orbital_truncation_level % 2 != 0 or self.orbital_truncation_level < 2:
#             raise BadValueError('Orbital truncation level must be an even integer >= 2')
#         if self.orbital_truncation_level > 6:
#             raise ImplementationException('Orbital truncation level > 6 has not been fully implemented in TidalPy.')
#
#         # Setup Love number calculator
#         self.full_calculation = True
#         if self.model == 'regular':
#             # FIXME: I don't like how the model is stored for the tidal calculations.
#             #    right now this stuff is stored in the /rheology/default.py which doesn't make a lot of sense.
#             pass
#         elif self.model == 'off':
#             self.full_calculation = False
#         elif self.model == 'multilayer':
#             # TODO: Multilayer code
#             raise ImplementationException
#         else:
#             raise UnknownModelError
#
#         # Setup complex compliance calculator
#         model_searcher = ComplianceModelSearcher(compliance_models, andrade_frequency_models, self.config,
#                                                  defaults_require_key=False)
#
#         if not self.full_calculation:
#             self.rheology_name = 'off'
#         else:
#             self.rheology_name = self.config['compliance_model']
#
#         comp_func, comp_input, comp_live_args = model_searcher.find_model(self.rheology_name)
#         # TODO: As of 0.1.0a no compliance functions have live args, so comp_live_args is never used.
#         self.compliance_func = comp_func
#         self.compliance_inputs = comp_input
#
#         # Switches
#         self.is_spin_sync = self.layer.is_spin_sync
#
#     def _calculate(self) -> Tuple[ArrayNone, ArrayNone]:
#         """ Calculates the Tidal Heating and Tidal Torque for this layer.
#
#         Returns
#         -------
#         total_heating : ArrayNone
#             Total Tidal Heating [W]
#         total_torque : ArrayNone
#             Total Tidal Torque [N m]
#
#         """
#
#         # Physical and Frequency Data
#         try:
#             shear = self.layer.shear_modulus
#             viscosity = self.layer.viscosity
#             orbital_freq = self.world.orbital_freq
#             spin_freq = self.world.spin_freq
#             eccentricity = self.world.eccentricity
#             inclination = self.world.inclination
#         except AttributeNotSetError:
#             raise ParameterMissingError
#         if shear is None:
#             raise ParameterMissingError
#         if viscosity is None:
#             raise ParameterMissingError
#         if orbital_freq is None:
#             raise ParameterMissingError
#         if spin_freq is None:
#             if self.is_spin_sync:
#                 spin_freq = orbital_freq
#             else:
#                 raise ParameterMissingError
#         if eccentricity is None:
#             eccentricity = np.asarray([0])
#         if inclination is None:
#             inclination = np.asarray([0])
#
#         if self.full_calculation:
#             # Make a call to the calculate tides function
#             tidal_heating, tidal_torque = \
#                 calculate_tides(shear, viscosity, self.layer.gravity, self.layer.radius, self.layer.density,
#                                 self.world.tidal_susceptibility, self.compliance_func, self.compliance_inputs,
#                                 self.world.semi_major_axis, eccentricity, inclination,
#                                 orbital_freq, spin_freq,
#                                 tidal_volume_fraction=self.layer.tidal_scale, use_nsr=self.is_spin_sync,
#                                 truncation=self.orbital_truncation_level, order_l=self.order_l)
#         else:
#             tidal_heating, tidal_torque = None, None
#
#         return tidal_heating, tidal_torque
#
