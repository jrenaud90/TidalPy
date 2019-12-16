import copy
from typing import Tuple, TYPE_CHECKING

import numpy as np

from ..constants import G
from . import calculate_tides
from ..exceptions import (AttributeNotSetError, ImplementationException, ParameterMissingError, UnknownModelError,
                          BadValueError, AttributeChangeRequiresReINIT, ImproperAttributeHandling, ParameterValueError)
from TidalPy.rheology.andrade_frequency import andrade_frequency_models
from TidalPy.rheology.compliance import compliance_models
from ..rheology.compliance import ComplianceModelSearcher
from ..rheology.defaults import rheology_param_defaults
from ..dynamics import mode_types
from ..utilities.classes import WorldConfigHolder
from ..types import ArrayNone
from ..utilities.model import LayerModelHolder
from .calculation import calc_tidal_susceptibility

if TYPE_CHECKING:
    from ..structures import TidalWorld, WorldType


FAKE_FREQS = (np.asarray([0.]),)


class Tides(WorldConfigHolder):

    """ Tides class - Holder for all tidal heating and torque calculation functions

    Tides class stores model parameters and methods for heating and torque which are functions of
        (T, P, melt_frac, w, e, I)
    """

    # No defaults or config information needed in this version of TidalPy.tides
    #    I am keeping tides a child class of configholder in case there is ever a reason to add in config info - it
    #    will be easy to do so.
    default_config = None
    layer_config_key = None

    def __init__(self, world: TidalWorld, store_config_in_world: bool = True):

        super().__init__(world, store_config_in_world)

        # State properties
        self._tidal_susceptibility = None
        self._tidal_susceptibility_reduced = None

        # Pull out configurations
        self._use_nsr = self.config['use_nsr']
        self._max_tidal_order_l = self.config['max_tidal_order_l']
        self._orbital_truncation = self.config['orbital_truncation']

        # Ensure the tidal order and orbital truncation levels make sense
        if self.max_tidal_order_l > 3:
            raise ImplementationException(f'Tidal order {self.max_tidal_order_l} has not been implemented yet.')
        if self.orbital_truncation % 2 != 0:
            raise ParameterValueError('Orbital truncation level must be an even integer.')
        if self.orbital_truncation <= 2:
            raise ParameterValueError('Orbital truncation level must be greater than or equal to 2.')
        if self.orbital_truncation not in [2, 4, 6]:
            raise ImplementationException(f'Orbital truncation level of {self.orbital_truncation} is not currently '
                                          f'supported.')
        self.order_l_list = range(2, self.max_tidal_order_l + 1)

        # Find the correct mode functions for heating and torque
        self.tidal_mode_calculators = [mode_types[self.use_nsr][l][self.orbital_truncation] for l in self.order_l_list]

    def calc_tidal_susceptibility_reduced(self) -> float:
        """ Reduced version of tidal susceptibility that does not take into account the orbital seperation between the
        target body and host.

        Returns
        -------
        tidal_susceptibility_reduced: float
            Reduced tidal susceptibility [N m]
        """

        self._tidal_susceptibility_reduced = (3. / 2.) * G * self.host.mass**2 * self.world.radius**5

        return self._tidal_susceptibility_reduced

    def calc_tidal_susceptibility(self) -> np.ndarray:
        """ Tidal susceptibility is the relative strength of tides based on a planet's size, its host's mass, and their
            orbital seperation.

        Tidal susceptibility is equal to tidal_susceptibility_reduced / a^6

        Returns
        -------
        tidal_susceptibility : np.ndarray
            Tidal susceptibility of this planet.
        """

        if self.tidal_susceptibility_reduced is None:
            self.calc_tidal_susceptibility_reduced()

        self._tidal_susceptibility = self.tidal_susceptibility_reduced / self.semi_major_axis**6

        return self._tidal_susceptibility

    def calc_tidal_modes(self):
        """ Calculates the tidal modes, frequencies, and heating & torque coefficients based on the planet's spin-rate
            and orbital motion.

        Returns
        -------


        """

        mode_names = list()
        modes = list()
        freqs = list()
        heating_subterms = list()

        tidal_modes = list()
        tidal

        mode_names, modes, freqs, heating_subterms, ztorque_subterms


    # Outerscope properties
    @property
    def host(self):
        if self.world.host is None:
            raise AttributeNotSetError
        return self.world.host

    @host.setter
    def host(self, value):
        raise ImproperAttributeHandling

    @property
    def semi_major_axis(self):
        if self.world.semi_major_axis is None:
            raise AttributeNotSetError
        return self.world.semi_major_axis

    @semi_major_axis.setter
    def semi_major_axis(self, value):
        raise ImproperAttributeHandling

    # State properties
    @property
    def use_nsr(self):
        return self._use_nsr

    @use_nsr.setter
    def use_nsr(self, value):
        raise AttributeChangeRequiresReINIT

    @property
    def max_tidal_order_l(self):
        return self._max_tidal_order_l

    @max_tidal_order_l.setter
    def max_tidal_order_l(self, value):
        raise AttributeChangeRequiresReINIT

    @property
    def orbital_truncation(self):
        return self._orbital_truncation

    @orbital_truncation.setter
    def orbital_truncation(self, value) -> float:
        raise AttributeChangeRequiresReINIT

    @property
    def tidal_susceptibility_reduced(self):
        return self._tidal_susceptibility_reduced

    @tidal_susceptibility_reduced.setter
    def tidal_susceptibility_reduced(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_susceptibility(self) -> np.ndarray:
        return self._tidal_susceptibility

    @tidal_susceptibility.setter
    def tidal_susceptibility(self, value):
        raise ImproperAttributeHandling



    class Rheology(LayerConfigHolder):

class Tides(LayerModelHolder):
    default_config = copy.deepcopy(rheology_param_defaults)
    config_key = 'rheology'

    def __init__(self, layer):

        super().__init__(layer=layer, function_searcher=None, call_reinit=True)

        # Override model if layer has tidal off
        if not self.layer.config['is_tidal']:
            self.model = 'off'

        # Pull out model information
        self.orbital_truncation_level = self.world.config['orbital_truncation_level']
        self.order_l = self.world.config['tidal_order_l']

        # Pull out information about the planet/layer
        self.config['planet_beta'] = self.layer.world.beta
        self.config['quality_factor'] = self.layer.world.fixed_q

        if int(self.order_l) != self.order_l or self.order_l < 2:
            raise BadValueError('Tidal order l must be an integer >= 2')
        if self.order_l > 3:
            raise ImplementationException('Tidal order l > 3 has not been fully implemented in TidalPy.')

        if int(self.orbital_truncation_level) != self.orbital_truncation_level or \
                self.orbital_truncation_level % 2 != 0 or self.orbital_truncation_level < 2:
            raise BadValueError('Orbital truncation level must be an even integer >= 2')
        if self.orbital_truncation_level > 6:
            raise ImplementationException('Orbital truncation level > 6 has not been fully implemented in TidalPy.')

        # Setup Love number calculator
        self.full_calculation = True
        if self.model == 'regular':
            # FIXME: I don't like how the model is stored for the tidal calculations.
            #    right now this stuff is stored in the /rheology/default.py which doesn't make a lot of sense.
            pass
        elif self.model == 'off':
            self.full_calculation = False
        elif self.model == 'multilayer':
            # TODO: Multilayer code
            raise ImplementationException
        else:
            raise UnknownModelError

        # Setup complex compliance calculator
        model_searcher = ComplianceModelSearcher(compliance_models, andrade_frequency_models, self.config,
                                                 defaults_require_key=False)

        if not self.full_calculation:
            self.rheology_name = 'off'
        else:
            self.rheology_name = self.config['compliance_model']

        comp_func, comp_input, comp_live_args = model_searcher.find_model(self.rheology_name)
        # TODO: As of 0.1.0a no compliance functions have live args, so comp_live_args is never used.
        self.compliance_func = comp_func
        self.compliance_inputs = comp_input

        # Switches
        self.is_spin_sync = self.layer.is_spin_sync

    def _calculate(self) -> Tuple[ArrayNone, ArrayNone]:
        """ Calculates the Tidal Heating and Tidal Torque for this layer.

        Returns
        -------
        total_heating : ArrayNone
            Total Tidal Heating [W]
        total_torque : ArrayNone
            Total Tidal Torque [N m]

        """

        # Physical and Frequency Data
        try:
            shear = self.layer.shear_modulus
            viscosity = self.layer.viscosity
            orbital_freq = self.world.orbital_freq
            spin_freq = self.world.spin_freq
            eccentricity = self.world.eccentricity
            inclination = self.world.inclination
        except AttributeNotSetError:
            raise ParameterMissingError
        if shear is None:
            raise ParameterMissingError
        if viscosity is None:
            raise ParameterMissingError
        if orbital_freq is None:
            raise ParameterMissingError
        if spin_freq is None:
            if self.is_spin_sync:
                spin_freq = orbital_freq
            else:
                raise ParameterMissingError
        if eccentricity is None:
            eccentricity = np.asarray([0])
        if inclination is None:
            inclination = np.asarray([0])

        if self.full_calculation:
            # Make a call to the calculate tides function
            tidal_heating, tidal_torque = \
                calculate_tides(shear, viscosity, self.layer.gravity, self.layer.radius, self.layer.density,
                                self.world.tidal_susceptibility, self.compliance_func, self.compliance_inputs,
                                self.world.semi_major_axis, eccentricity, inclination,
                                orbital_freq, spin_freq,
                                tidal_volume_fraction=self.layer.tidal_scale, use_nsr=self.is_spin_sync,
                                truncation=self.orbital_truncation_level, order_l=self.order_l)
        else:
            tidal_heating, tidal_torque = None, None

        return tidal_heating, tidal_torque
