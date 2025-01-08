from typing import List, TYPE_CHECKING, Union

import numpy as np

import TidalPy
from TidalPy.logger import get_logger
from TidalPy.exceptions import (AttributeNotSetError, ConfigPropertyChangeError, IncorrectAttributeType,
                                IncorrectMethodToSetStateProperty, OuterscopePropertySetError, ParameterMissingError,
                                UnknownModelError)
from TidalPy.utilities.classes.model import LayerModelHolder

from . import known_model_const_args, known_model_live_args, known_models

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray, NoneType
    from TidalPy.structures.layers import PhysicalLayerType

log = get_logger("TidalPy")


class Radiogenics(LayerModelHolder):
    """ Radiogenic Model Class - Child of LayerModelHolder Class

    Radiogenic model provides the functionality to calculate a layer's heating due to radioactive isotopes based on
        user provided parameters related to convection and conduction.
    """

    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = 'radiogenics'

    def __init__(
        self, layer: 'PhysicalLayerType', model_name: str = None, store_config_in_layer: bool = True,
        initialize: bool = True
        ):
        """ Constructor for the Radiogenics model

        Parameters
        ----------
        layer : PhysicalLayerType
            The layer instance which radiogenics are calculated.
        model_name : str = None
            The user-provided radiogenic model name.
        store_config_in_layer: bool = True
            Flag that determines if the final radiogenic configuration dictionary should be copied into the
            `layer.config` dictionary.
        initialize : bool = True
            Determines if initial reinit should be performed on the model (loading in data from its `self.config`).
        """

        # Set initialize to False so that the functions can be built at the end of this __init__ once a few
        #    more parameters are loaded into the class's config.
        super().__init__(layer, model_name, store_config_in_layer, initialize=False)

        # State properties
        self._heating = None

        # Configuration properties
        self._isos_name = None
        self._isos_hpr = None
        self._isos_halflife = None
        self._isos_massfrac = None
        self._isos_concentration = None
        self._radiogenic_layer_mass_frac = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        """ Reinitialization for the Radiogenic model
        Model will look at the user-provided configurations and pull out model information including constants

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this method has been called (additional steps may be
                preformed during the first reinit call).
        """

        # The model builder only looks for 'ref_time' if the user provided the reference time as 'reference_time'
        #    add the 'ref_time' key pointing to the same value.
        if 'ref_time' not in self.config:
            if 'reference_time' in self.config:
                self.config['ref_time'] = self.config['reference_time']
            else:
                self.config['ref_time'] = None

        # User can reduce the amount of the layer's mass used for radiogenic calculations if they have a reason to
        #   think that a portion of the layer contains no radiogenics (like a molten portion).
        # TODO: It would be nice to have this as a dynamic variable so that if a layer is melting it may be
        #    concentrating radiogenics into a specific portion of the layer
        self._radiogenic_layer_mass_frac = self.config['radiogenic_layer_mass_fraction']

        if self.model == 'isotope':
            # Reset the isotope lists back to empty
            self._isos_name = list()
            self._isos_hpr = list()
            self._isos_halflife = list()
            self._isos_massfrac = list()
            self._isos_concentration = list()

            # Isotopes may be given as a dictionary of individual isotopes or as a string pointing to one of the
            #  pre-built TidalPy isotope lists.
            isotopes = self.config['isotopes']
            if type(isotopes) == str:
                if isotopes.lower() not in TidalPy.config['physics']['radiogenics']['known_isotope_data']:
                    raise UnknownModelError
                iso_datas = TidalPy.config['physics']['radiogenics']['known_isotope_data'][isotopes]
            else:
                iso_datas = isotopes

            # Different Isotope data sources may have their own reference time - extract that information and
            #    store it in the model's main configuration dictionary.
            new_ref_time = None
            if 'ref_time' in iso_datas:
                new_ref_time = iso_datas['ref_time']
            elif 'reference_time' in iso_datas:
                new_ref_time = iso_datas['reference_time']

            # If no reference time was found use the model config's
            if new_ref_time is None:
                new_ref_time = self.config['ref_time']
            else:
                log.debug(f'Overriding default reference time with value provided by isotope data in {self}.')

            if new_ref_time is None:
                # If there is still no reference time then fall back on a default, but warn the user.
                log.warning(f'No reference time provided for radiogenics, using ref_time = 0 Myr in {self}.')
                self.config['ref_time'] = 0.
            else:
                self.config['ref_time'] = new_ref_time

            # Build isotope inputs and store them in radiogenics.config - these will end up in the self.inputs used in
            #    self.calculate
            for isotope, iso_data in iso_datas.items():
                if isotope in ['ref_time', 'reference_time']:
                    # Skip non-isotope keys
                    continue

                # For each isotope, extract the needed info
                self._isos_name.append(isotope)
                try:
                    self._isos_hpr.append(iso_data['hpr'])
                    self._isos_halflife.append(iso_data['half_life'])
                    self._isos_massfrac.append(iso_data['element_concentration'])
                    self._isos_concentration.append(iso_data['iso_mass_fraction'])
                except KeyError:
                    raise ParameterMissingError(f'One or more parameters are missing for isotope {isotope} in {self}.')

            # Add the isotope data to the config so that the argument builder can find them
            self.config['iso_massfracs_of_isotope'] = self.isos_massfrac
            self.config['iso_element_concentrations'] = self.isos_concentration
            self.config['iso_halflives'] = self.isos_halflife
            self.config['iso_heat_production'] = self.isos_hpr

        # Call the parent's reinit right away
        super().reinit(initial_init)

    def clear_state(self):
        """ Clear the Radiogenic model's state. """

        super().clear_state()

        self._heating = None

    def _calculate(self) -> 'FloatArray':
        """ Calculates the radiogenic heating of layer in which the radiogenic class is installed.

        Returns
        -------
        radiogenic_heating : FloatArray
            Radiogenic heating [W]
        """

        radiogenic_heating = self.func_array(*self.live_inputs, *self.inputs)
        self._heating = radiogenic_heating

        return radiogenic_heating

    def _calculate_debug(self) -> np.ndarray:

        # Attribute checks
        if self.time is None:
            raise AttributeNotSetError
        if self.mass is None:
            raise AttributeNotSetError
        if type(self.time) != np.ndarray:
            raise IncorrectAttributeType
        if type(self.mass) != float:
            raise IncorrectAttributeType

        # Value checks
        if self.model == 'isotope':
            if abs(self.time[-1] / max(self.isos_halflife)) > 1.0e4:
                log.warning(
                    'Time is much larger than radiogenic half-life - Check units of the time array and half lives.'
                    )
        elif self.model == 'fixed':
            if abs(self.time[-1] / max(self.config['average_half_life'])) > 1.0e4:
                log.warning(
                    'Time is much larger than the fixed-average radiogenic half-life - '
                    'Check units of the time array and fixed half-life.'
                    )

        # Calculate and perform more value checks
        radiogenic_heating = self._calculate()
        if np.any(radiogenic_heating < 0.):
            log.warning(f'Negative radiogenic heating encountered at time:\n{self.time[radiogenic_heating < 0.]}')
        if np.any(radiogenic_heating > 1.e23):
            log.warning(f'Very large radiogenic heating encountered at time:\n{self.time[radiogenic_heating > 1.e23]}')

        return radiogenic_heating

    # # Configuration properties
    @property
    def isos_name(self) -> Union['NoneType', List[str]]:
        """ List of isotope names used in radiogenic calculations (for isotope model only) """
        return self._isos_name

    @isos_name.setter
    def isos_name(self, value):
        raise ConfigPropertyChangeError

    @property
    def isos_hpr(self) -> Union['NoneType', List[float]]:
        """ List of isotope heat production rates [W kg-1] (for isotope model only) """
        return self._isos_hpr

    @isos_hpr.setter
    def isos_hpr(self, value):
        raise ConfigPropertyChangeError

    @property
    def isos_halflife(self) -> Union['NoneType', List[float]]:
        """ List of isotope half lives [Myr] (for isotope model only) """
        return self._isos_halflife

    @isos_halflife.setter
    def isos_halflife(self, value):
        raise ConfigPropertyChangeError

    @property
    def isos_massfrac(self) -> Union['NoneType', List[float]]:
        """ List of isotope mass fractions [kg kg-1] (for isotope model only)

        Notes
        -----
        .. This is the mass fraction of the specific radiogenic isotope in a block of all isotopes of the particular
            element.

        See Also
        --------
        Radiogenics.isos_massfrac
        """
        return self._isos_massfrac

    @isos_massfrac.setter
    def isos_massfrac(self, value):
        raise ConfigPropertyChangeError

    @property
    def isos_concentration(self) -> Union['NoneType', List[float]]:
        """ List of element concentration [kg kg-1] (for isotope model only)

        Notes
        -----
        .. This is the mass concentration of the isotope's element (e.g., Uranium) in a block of generic layer
            material. aka: how many kg of uranium is in a kg of silicate rock. This is not the amount of the specific
            radiogenic isotope (see isos_massfrac).

        See Also
        --------
        Radiogenics.isos_massfrac
        """
        return self._isos_concentration

    @isos_concentration.setter
    def isos_concentration(self, value):
        raise ConfigPropertyChangeError

    @property
    def radiogenic_layer_mass_frac(self) -> float:
        """ Fraction of the layer's mass where radiogenic isotopes are concentrated (defaults to 1) """
        return self._radiogenic_layer_mass_frac

    @radiogenic_layer_mass_frac.setter
    def radiogenic_layer_mass_frac(self, value):
        raise ConfigPropertyChangeError

    # # State properties
    @property
    def heating(self) -> 'FloatArray':
        """ Radiogenic heating rate [W] """
        return self._heating

    @heating.setter
    def heating(self, value):
        # TODO We could have the user set the radiogenic heating and then this setter could call the layer's
        #  thermal update.
        raise IncorrectMethodToSetStateProperty

    # # Outer-scope properties
    @property
    def time(self):
        """ Outer-scope Wrapper for world.time """
        return self.world.time

    @time.setter
    def time(self, value):
        raise OuterscopePropertySetError

    @property
    def mass(self):
        """ Outer-scope Wrapper for layer.mass including a scale for the radiogenic mass frac """
        return self.layer.mass * self.radiogenic_layer_mass_frac

    @mass.setter
    def mass(self, value):
        raise OuterscopePropertySetError
