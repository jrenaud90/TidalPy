""" Tides Module
"""

import copy
from typing import TYPE_CHECKING, List, Dict

import numpy as np

from .defaults import tide_defaults
from .dissipation import calc_tidal_susceptibility, calc_tidal_susceptibility_reduced, mode_collapse
from ..constants import G
from ..exceptions import (AttributeNotSetError, ImplementationException, ParameterMissingError, UnknownModelError,
                          BadValueError, AttributeChangeRequiresReINIT, ImproperAttributeHandling, ParameterValueError,
                          OuterscopeAttributeSetError, ConfigAttributeChangeError, MissingArgumentError)

from ..types import FloatArray
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

    default_config = tide_defaults
    world_config_key = 'tides'

    def __init__(self, world: 'TidalWorld', store_config_in_world: bool = True):

        super().__init__(world, store_config_in_world=store_config_in_world)

        # Model setup
        self._eccentricity_truncation_lvl = self.config['eccentricity_truncation_lvl']
        self._tidal_order_lvl = self.config['max_tidal_order_l']
        self._use_obliquity_tides = self.config['obliquity_tides_on']

        # Ensure the tidal order and orbital truncation levels make sense
        if self._tidal_order_lvl > 5:
            raise ImplementationException(f'Tidal order {self._tidal_order_lvl} has not been implemented yet.')
        if self._eccentricity_truncation_lvl%2 != 0:
            raise ParameterValueError('Orbital truncation level must be an even integer.')
        if self._eccentricity_truncation_lvl <= 2:
            raise ParameterValueError('Orbital truncation level must be greater than or equal to 2.')
        if self._eccentricity_truncation_lvl not in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]:
            raise ImplementationException(f'Orbital truncation level of {self._eccentricity_truncation_lvl} has not '
                                          f'been implemented yet.')

        # State properties
        self._tidal_susceptibility = None
        self._tidal_susceptibility_reduced = None
        self._unique_tidal_frequencies = None
        self._tidal_terms_by_frequency = None
        self._tidal_heating_by_layer = None
        self._negative_imk_by_layer = None
        self._tidal_heating_global = None
        self._dUdM = None
        self._dUdw = None
        self._dUdO = None

    def initialize_tides(self):
        """ Initialize various tidal parameters once a tidal host is connected to the target body.
        """

        if self.tidal_host is None:
            raise AttributeNotSetError('Tidal host must be connected to target body in order to initialize tides.')

        self._tidal_susceptibility_reduced = \
            calc_tidal_susceptibility_reduced(self.tidal_host.mass, self.world.radius)

        if self.semi_major_axis is not None:
            self._tidal_susceptibility = calc_tidal_susceptibility(self.tidal_host.mass, self.world.radius,
                                                                   self.semi_major_axis)

    def orbital_change(self):
        """ Calculate tidal heating and potential derivative terms based on the current orbital state.

        This will also calculate new unique tidal frequencies which must then be digested by the rheological model
            at each planetary layer.

        See Also
        --------
        TidalPy.tides.dissipation.mode_collapse
        """

        # Calculate the various terms based on the current state
        results_by_uniquefreq = mode_collapse(self.orbital_frequency, self.spin_frequency,
                                              self.eccentricity, self.obliquity,
                                              use_obliquity=self.use_obliquity_tides,
                                              eccentricity_truncation_lvl=self.eccentricity_truncation_lvl,
                                              max_order_l=self.tidal_order_lvl)

        # Pull out unique frequencies and tidal terms. Clear old storage.
        self._unique_tidal_frequencies = dict()
        self._tidal_terms_by_frequency = dict()
        for unique_freq_signature, (unique_frequency, tidal_terms_by_orderl) in results_by_uniquefreq:
            self._unique_tidal_frequencies[unique_freq_signature] = unique_frequency
            self._tidal_terms_by_frequency[unique_freq_signature] = tidal_terms_by_orderl

    def thermal_change(self):
        """ Calculate Global Love number based on current thermal state.

        The requires a prior orbital_change() call as unique frequencies are used to calculate the complex compliances
            used to calculate the Love numbers.

        See Also
        --------
        TidalPy.tides.Tides.orbital_change
        """

        # TODO: Can any of this be pulled out into a njit'd function to improve performance?

        global_love_number_terms = list()
        tidal_heating_reduced_terms = list()
        dUdM_reduced_terms = list()
        dUdw_reduced_terms = list()
        dUdO_reduced_terms = list()

        for layer in self.world:
            # Check to see if the layer is contributing to tides
            if not layer.is_tidal:
                continue

            # The love number will be reduced by each layer's tidal volume fraction or tidal scale
            layer_tidal_scale = layer.tidal_scale

            # Pull out other parameters used in the calculations
            # TODO: These are used to calculate the effective rigidity. Should these be for the layer or for the planet
            #    as a whole?
            shear_modulus = layer.shear_modulus
            radius = layer.radius
            gravity = layer.gravity_surface
            density = layer.density_bulk

            tidal_heating_reduced_terms_for_layer = list()
            negative_imk_for_layer = list()

            # For each tidal order number calculate the complex love number and use it to make the first collapse
            #    on tidal terms.
            for tidal_order_l in range(2, self.tidal_order_lvl + 1):

                # Calculate the effective rigidity which does not use the complex compliance
                effective_rigidity = effective_rigidity_general(shear_modulus, gravity, radius, density,
                                                                order_l=tidal_order_l)

                # Pull out the already computed complex compliances for each frequency
                for unique_freq_signature, complex_compliance in layer.complex_compliance_by_frequency.items():

                    complex_love = complex_love_general(complex_compliance, shear_modulus, effective_rigidity,
                                                        order_l=tidal_order_l)

                    # Scale the Love number by the layer's contribution
                    # TODO: Should the tidal scale affect the Re(k) as well as the Im(k)?
                    neg_imk = np.imag(complex_love) * layer_tidal_scale
                    neg_imk_scaled = neg_imk * self.tidal_susceptibility

                    # The tidal potential carries one fewer dependence on the tidal host mass that is built into the
                    #    tidal susceptibility. Divide that out now.
                    neg_imk_scaled_potential = neg_imk_scaled / self.tidal_host.mass

                    # Pull out the tidal terms pre-calculated for this unique frequency. See Tides.orbital_change
                    heating_term, dUdM_term, dUdw_term, dUdO_term = self.tidal_terms_by_frequency[unique_freq_signature]

                    # Store results
                    tidal_heating_reduced_terms.append(heating_term * neg_imk_scaled)
                    dUdM_reduced_terms.append(dUdM_term * neg_imk_scaled_potential)
                    dUdw_reduced_terms.append(dUdw_term * neg_imk_scaled_potential)
                    dUdO_reduced_terms.append(dUdO_term * neg_imk_scaled_potential)

                    global_love_number_terms.append(complex_love)
                    negative_imk_for_layer.append(neg_imk_scaled)

            # The layer needs to know what the tidal heating is within it. Store this now.
            self.tidal_heating_by_layer[layer] = sum(tidal_heating_reduced_terms_for_layer)
            self.negative_imk_by_layer[layer] = sum(negative_imk_for_layer)

        # Collapse all unqiue frequencies and tidal order-l terms down into a single value (or array)
        self._tidal_heating_global = sum(tidal_heating_reduced_terms)
        self._dUdM = sum(dUdM_reduced_terms)
        self._dUdw = sum(dUdw_reduced_terms)
        self._dUdO = sum(dUdO_reduced_terms)

        return self.tidal_heating_global, self.dUdM, self.dUdw, self.dUdO

    @staticmethod
    def calc_tidal_susceptibility(host_mass: float, target_radius: float, semi_major_axis: FloatArray) -> FloatArray:
        """ Calculate the tidal susceptibility for a target object orbiting

        Wrapper for dissipation.py/calc_tidal_susceptibility

        Parameters
        ----------
        host_mass : float
            Tidal host's mass [kg]
        target_radius : float
            Target body's mean radius [m]
        semi_major_axis : FloatArray
            Orbital separation between the target and host [m]

        Returns
        -------
        tidal_susceptibility : FloatArray
            Tidal Susceptibility [N m]
        """

        tidal_susceptibility = calc_tidal_susceptibility(host_mass, target_radius, semi_major_axis)
        return tidal_susceptibility

    # Configuration properties
    @property
    def eccentricity_truncation_lvl(self) -> int:
        return self._eccentricity_truncation_lvl

    @eccentricity_truncation_lvl.setter
    def eccentricity_truncation_lvl(self, value):
        raise ConfigAttributeChangeError

    @property
    def tidal_order_lvl(self) -> int:
        return self._tidal_order_lvl

    @tidal_order_lvl.setter
    def tidal_order_lvl(self, value):
        raise ConfigAttributeChangeError

    @property
    def use_obliquity_tides(self) -> bool:
        return self._use_obliquity_tides

    @use_obliquity_tides.setter
    def use_obliquity_tides(self, value):
        raise ConfigAttributeChangeError


    # State properties
    @property
    def tidal_susceptibility_reduced(self) -> FloatArray:
        return self._tidal_susceptibility_reduced

    @tidal_susceptibility_reduced.setter
    def tidal_susceptibility_reduced(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_susceptibility(self) -> FloatArray:
        return self._tidal_susceptibility

    @tidal_susceptibility.setter
    def tidal_susceptibility(self, value):
        raise ImproperAttributeHandling

    @property
    def unique_tidal_frequencies(self) -> Dict[str, FloatArray]:
        return self._unique_tidal_frequencies

    @unique_tidal_frequencies.setter
    def unique_tidal_frequencies(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_terms_by_frequency(self) -> Dict[str, Dict[int, FloatArray]]:
        return self._tidal_terms_by_frequency

    @tidal_terms_by_frequency.setter
    def tidal_terms_by_frequency(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_heating_by_layer(self) -> Dict[ThermalLayer, FloatArray]:
        return self._tidal_heating_by_layer

    @tidal_heating_by_layer.setter
    def tidal_heating_by_layer(self, value):
        raise ImproperAttributeHandling

    @property
    def negative_imk_by_layer(self) -> Dict[ThermalLayer, FloatArray]:
        return self._negative_imk_by_layer

    @negative_imk_by_layer.setter
    def negative_imk_by_layer(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_heating_global(self) -> FloatArray:
        return self._tidal_heating_global

    @tidal_heating_global.setter
    def tidal_heating_global(self, value):
        raise ImproperAttributeHandling

    @property
    def dUdM(self) -> FloatArray:
        return self._dUdM

    @dUdM.setter
    def dUdM(self, value):
        raise ImproperAttributeHandling

    @property
    def dUdw(self) -> FloatArray:
        return self._dUdw

    @dUdw.setter
    def dUdw(self, value):
        raise ImproperAttributeHandling

    @property
    def dUdO(self) -> FloatArray:
        return self._dUdO

    @dUdO.setter
    def dUdO(self, value):
        raise ImproperAttributeHandling


    # Outer-scope Properties
    @property
    def semi_major_axis(self) -> FloatArray:
        return self.world.semi_major_axis

    @semi_major_axis.setter
    def semi_major_axis(self, value):
        raise OuterscopeAttributeSetError

    @property
    def orbital_frequency(self) -> FloatArray:
        return self.world.orbital_frequency

    @orbital_frequency.setter
    def orbital_frequency(self, value):
        raise OuterscopeAttributeSetError

    @property
    def spin_frequency(self) -> FloatArray:
        return self.world.spin_frequency

    @spin_frequency.setter
    def spin_frequency(self, value):
        raise OuterscopeAttributeSetError

    @property
    def eccentricity(self) -> FloatArray:
        return self.world.eccentricity

    @eccentricity.setter
    def eccentricity(self, value):
        raise OuterscopeAttributeSetError

    @property
    def obliquity(self):
        if self.use_obliquity_tides:
            return self.world.obliquity
        else:
            return np.zeros_like(self.eccentricity)

    @obliquity.setter
    def obliquity(self, value):
        raise OuterscopeAttributeSetError

    @property
    def tidal_host(self):
        return self.world.tidal_host

    @tidal_host.setter
    def tidal_host(self, value):
        raise OuterscopeAttributeSetError