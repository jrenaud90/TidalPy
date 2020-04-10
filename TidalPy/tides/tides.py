import copy
from typing import TYPE_CHECKING, List, Dict

import numpy as np

from ..constants import G
from ..exceptions import (AttributeNotSetError, ImplementationException, ParameterMissingError, UnknownModelError,
                          BadValueError, AttributeChangeRequiresReINIT, ImproperAttributeHandling, ParameterValueError,
                          OuterscopeAttributeSetError, ConfigAttributeChangeError, MissingArgumentError)

from ..utilities.classes import WorldConfigHolder, LayerConfigHolder
from ..types import ArrayNone
from ..utilities.model import LayerModelHolder
from .love1d import complex_love_general, effective_rigidity_general

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
        return self.world.use_nsr

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
