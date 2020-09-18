""" Layered Tides Module
"""

from typing import TYPE_CHECKING, Dict

import numpy as np

from .base import TidesBase
from .defaults import tide_defaults
from ..mode_manipulation import DissipTermsArray
from ...exceptions import (NotYetImplementedError, ParameterValueError,
                           FailedForcedStateUpdate,
                           BadAttributeValueError, IncorrectMethodToSetStateProperty)
from ...utilities.types import FloatArray

if TYPE_CHECKING:
    from ...structures.worlds import LayeredWorld
    from ...structures.layers import PhysicsLayer


class LayeredTides(TidesBase):
    """ LayeredTides
    Class used for layered planets (icy or rocky worlds)

    Tides class stores model parameters and methods for heating and torque which are functions of
        (T, P, melt_frac, w, e, theata)

    See Also
    --------
    TidalPy.tides.classes.TidesBase
    """

    default_config = tide_defaults['layered']

    def __init__(self, world: 'LayeredWorld', store_config_in_world: bool = True, initialize: bool = True):
        """ Constructor for TidesBase class

        Parameters
        ----------
        world : TidalWorldType
            The world where tides are being calculated.
        store_config_in_world : bool = True
            Flag that determines if the final model's configuration dictionary should be copied into the
            `world.config` dictionary.
        initialize : bool = True
            If `True`, then an initial call to the tide's reinit method will be made at the end of construction.
        """

        super().__init__(world, store_config_in_world=store_config_in_world, initialize=initialize)

        # State properties
        self._tidal_heating_by_layer = {layer: None for layer in self.world}
        self._negative_imk_by_layer = {layer: None for layer in self.world}

        # Ensure the tidal order and orbital truncation levels make sense
        if self.max_tidal_order_lvl > 7:
            raise NotYetImplementedError(f'Tidal order {self.max_tidal_order_lvl} has not been implemented yet.')
        if self.eccentricity_truncation_lvl % 2 != 0:
            raise ParameterValueError('Orbital truncation level must be an even integer.')
        if self.eccentricity_truncation_lvl <= 2:
            raise ParameterValueError('Orbital truncation level must be greater than or equal to 2.')
        if self.eccentricity_truncation_lvl not in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20):
            if self.max_tidal_order_lvl == 2 and self.eccentricity_truncation_lvl == 22:
                # This was implemented in v0.2.1
                pass
            else:
                raise NotYetImplementedError(f'Orbital truncation level of {self.eccentricity_truncation_lvl} has not '
                                              f'been implemented yet.')

        # Pull out information from the planet's layers
        self._tidal_inputs_by_layer = dict()
        for layer in self.world:
            if layer.is_tidal:
                tidal_scale = layer.tidal_scale
                # This system assumes that density, radius, and gravity will not change after initialization
                radius = layer.radius
                bulk_density = layer.density_bulk
                gravity_surf = layer.gravity_surface

                for param in [tidal_scale, radius, bulk_density, gravity_surf]:
                    if param is None:
                        # How did that happen?
                        raise BadAttributeValueError

                self._tidal_inputs_by_layer[layer] = (tidal_scale, radius, bulk_density, gravity_surf)
            else:
                # Layer does not contribute to tides. This will be marked by a None in this list
                self._tidal_inputs_by_layer[layer] = None

        # Pull out planet properties that may be used based on the configuration
        self._planet_tidal_inputs = None
        if self.config['use_planet_params_for_love_calc']:
            # TODO: These are used to calculate the effective rigidity. Should these be for the layer or for the planet
            #    as a whole?
            planet_radius = self.world.radius
            planet_gravity = self.world.gravity_surface
            planet_density = self.world.density_bulk
            self._planet_tidal_inputs = (planet_radius, planet_density, planet_gravity)

    def clear_state(self):
        """ Clear the state properties for the layered tides model """

        super().clear_state()

        # Clear tidal results stored for each layer
        self._tidal_heating_by_layer = {layer: None for layer in self.world}
        self._negative_imk_by_layer = {layer: None for layer in self.world}

    def collapse_modes(self, force_update: bool = True) -> DissipTermsArray:
        """ Calculate Global Love number based on current thermal state.

        Requires a prior update_orbit_spin() call as unique frequencies are used to calculate the complex compliances
            used to calculate the Love numbers.

        Parameters
        ----------
        force_update : bool = True
            If True, will raise an error if the update does not successful complete.
            Failed completions are usually a result of a missing state property(ies).

        Returns
        -------
        tidal_heating : np.ndarray
            Tidal heating [W]
            This could potentially restricted to a layer or for an entire planet.
        dUdM : np.ndarray
            Tidal potential derivative with respect to the mean anomaly [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.
        dUdw : np.ndarray
            Tidal potential derivative with respect to the pericentre [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.
        dUdO : np.ndarray
            Tidal potential derivative with respect to the planet's node [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.

        See Also
        --------
        TidalPy.tides.Tides.update_orbit_spin
        """

        # Check to see if all the needed state properties are present
        all_values_present = True

        # Check to see if tidal terms have been loaded.
        if self.tidal_terms_by_frequency is None:
            # Attempt to update the orbit to see if we can load them
            self.update_orbit_spin(force_update=False)

            if self.tidal_terms_by_frequency is None:
                # Update did not help.
                all_values_present = False

        if all_values_present:
            nonNone_love_number = list()
            nonNone_neg_imk = list()
            nonNone_tidal_heating = list()
            nonNone_dUdM = list()
            nonNone_dUdw = list()
            nonNone_dUdO = list()

            love_number_by_layer = dict()
            neg_imk_by_layer = dict()
            tidal_heating_by_layer = dict()
            dUdM_by_layer = dict()
            dUdw_by_layer = dict()
            dUdO_by_layer = dict()

            broke_out = False
            for layer, other_inputs in self._tidal_inputs_by_layer.items():
                if other_inputs is None:
                    # Not a tidal layer
                    love_number_by_layer[layer] = None
                    neg_imk_by_layer[layer] = None
                    tidal_heating_by_layer[layer] = None
                    dUdM_by_layer[layer] = None
                    dUdw_by_layer[layer] = None
                    dUdO_by_layer[layer] = None

                else:
                    tidal_scale, radius, bulk_density, gravity_surf = other_inputs
                    if self._planet_tidal_inputs is not None:
                        radius, bulk_density, gravity_surf = self._planet_tidal_inputs

                    # Pull out variables that change often
                    shear_modulus = layer.shear_modulus
                    complex_compliances_by_frequency_list = layer.complex_compliance_by_frequency_list

                    if shear_modulus is None or complex_compliances_by_frequency_list is None:
                        # uh oh
                        broke_out = True
                        break

                    # Mode collapse will parse through tidal order-l and all unique frequencies and calculate global and
                    #    localized dissipation values
                    tidal_heating, dUdM, dUdw, dUdO, love_number, negative_imk = \
                        self.collapse_modes_func(gravity_surf, radius, bulk_density, shear_modulus,
                                                 complex_compliances_by_frequency_list,
                                                 self.tidal_terms_by_frequency, self.tidal_susceptibility,
                                                 self.tidal_host.mass,
                                                 tidal_scale, cpl_ctl_method=False)

                    # These will be summed for global values
                    nonNone_love_number.append(love_number)
                    nonNone_neg_imk.append(negative_imk)
                    nonNone_tidal_heating.append(tidal_heating)
                    nonNone_dUdM.append(dUdM)
                    nonNone_dUdw.append(dUdw)
                    nonNone_dUdO.append(dUdO)

                    # Accessed by layer classes for thermal evolution
                    tidal_heating_by_layer[layer] = tidal_heating
                    neg_imk_by_layer[layer] = negative_imk

                    # TODO: These are not accessible at the moment. I suppose it would be useful info?
                    love_number_by_layer[layer] = love_number
                    dUdM_by_layer[layer] = dUdM
                    dUdw_by_layer[layer] = dUdw
                    dUdO_by_layer[layer] = dUdO

            if not broke_out:
                # Loop finished successfully. Store info in accessible containers
                self._tidal_heating_by_layer = tidal_heating_by_layer
                self._negative_imk_by_layer = neg_imk_by_layer

                self._tidal_heating_global = sum(nonNone_tidal_heating)
                self._dUdM = sum(nonNone_dUdM)
                self._dUdw = sum(nonNone_dUdw)
                self._dUdO = sum(nonNone_dUdO)
                self._negative_imk_global = sum(nonNone_neg_imk)

                # Now tell other methods to update now that derivatives and heating has been altered
                # TODO: orbit derivatives
                # TODO: layer thermal evolution
                self.set_spin_derivative()

                # Return tidal heating and derivatives
                return self.tidal_heating_global, self.dUdM, self.dUdw, self.dUdO

            else:
                if force_update:
                    raise FailedForcedStateUpdate

        else:
            if force_update:
                raise FailedForcedStateUpdate

    # # State properties
    @property
    def tidal_heating_by_layer(self) -> Dict['PhysicsLayer', FloatArray]:
        """ Tidal heating stored by layer instance [W] """
        return self._tidal_heating_by_layer

    @tidal_heating_by_layer.setter
    def tidal_heating_by_layer(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def negative_imk_by_layer(self) -> Dict['PhysicsLayer', FloatArray]:
        """ -Im[k2] stored by layer instance """
        return self._negative_imk_by_layer

    @negative_imk_by_layer.setter
    def negative_imk_by_layer(self, value):
        raise IncorrectMethodToSetStateProperty