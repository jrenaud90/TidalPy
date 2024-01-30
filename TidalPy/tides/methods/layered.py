""" Layered Tides Module
"""

from typing import Dict, List, TYPE_CHECKING

import numpy as np

from TidalPy.exceptions import (BadAttributeValueError, IncorrectMethodToSetStateProperty, MissingAttributeError)

from .base import TidesBase

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray
    from TidalPy.structures.world_types import LayeredWorld
    from TidalPy.structures.layers import PhysicalLayerType

    from ..modes.mode_manipulation import DissipTermsArray


class LayeredTides(TidesBase):
    """ LayeredTides
    Class used for layered planets (icy or rocky world_types)

    Tides class stores model parameters and methods for heating and torque which are functions of
        (T, P, melt_frac, w, e, theata)

    Attributes
    ----------
    tidal_heating_by_layer
    negative_imk_by_layer_by_orderl

    See Also
    --------
    TidalPy.tides.methods.TidesBase
    """

    model = 'layered'

    def __init__(self, world: 'LayeredWorld', store_config_in_world: bool = True, initialize: bool = True):
        """ Constructor for LayeredTides class

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

        super().__init__(world, store_config_in_world=store_config_in_world, initialize=False)

        # State properties
        self._tidal_heating_by_layer = {layer: None for layer in self.world}
        self._negative_imk_by_layer = {layer: None for layer in self.world}

        # Configuration properties
        self._tidal_input_getters_by_layer = dict()
        self._world_tidal_input_getters = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        """ Load configurations into the Tides class and import any config-dependent functions.

        This reinit process is separate from the __init__ method because the Orbit class may need to overload some
            configurations after class initialization.

        Parameters
        ----------
        initial_init : bool = False
            This should be set to True the first time reinit is called.

        Raises
        ------
        BadAttributeValueError

        """

        super().reinit(initial_init)

        # Reset configuration properties
        if not initial_init:
            self._tidal_input_getters_by_layer = dict()
            self._world_tidal_input_getters = None

        # Pull out tidal inputs
        for layer in self.world:
            if layer.is_tidal:
                get_tidal_scale = lambda: layer.tidal_scale
                # This system assumes that density, radius, and gravity will not change after initialization
                get_radius = lambda: layer.radius
                get_bulk_density = lambda: layer.density_bulk
                get_surf_gravity = lambda: layer.gravity_surface

                for param in [get_tidal_scale, get_radius, get_bulk_density, get_surf_gravity]:
                    if param is None:
                        # How did that happen?
                        raise BadAttributeValueError

                self._tidal_input_getters_by_layer[layer] = \
                    (get_tidal_scale, get_radius, get_bulk_density, get_surf_gravity)
            else:
                # Layer does not contribute to tides. This will be marked by a None in this list
                self._tidal_input_getters_by_layer[layer] = None

        # Pull out planet properties that may be used based on the configuration
        if self.config['use_planet_params_for_love_calc']:
            # TODO: These are used to calculate the effective rigidity. Should these be for the layer or for the planet
            #    as a whole?
            # TODO: Tidal scale for world? -> planet_tidal_scale = lambda: self.world.tidal
            get_world_radius = lambda: self.world.radius
            get_world_density = lambda: self.world.density_bulk
            get_world_gravity = lambda: self.world.gravity_surface
            self._world_tidal_input_getters = (get_world_radius, get_world_density, get_world_gravity)

    def complex_compliances_changed(self, collapse_tidal_modes: bool = True):
        """ The complex compliances have changed. Make any necessary updates.

        Parameters
        ----------
        collapse_tidal_modes : bool = True
            If `True`, then the world will tell its tides model to collapse tidal modes.
        """

        # This is called from bottom-to-top starting in the ComplexCompliances class inside Rheology.

        if collapse_tidal_modes:
            self.collapse_modes()

    def clear_state(self):
        """ Clear the state properties for the layered tides model """

        super().clear_state()

        # Clear tidal results stored for each layer
        self._tidal_heating_by_layer = {layer: None for layer in self.world}
        self._negative_imk_by_layer = {layer: None for layer in self.world}

    def collapse_modes(self) -> 'DissipTermsArray':
        """ Calculate Global Love number based on current thermal state.

        Requires a prior orbit_spin_changed() call as unique frequencies are used to calculate the complex compliances
            used to calculate the Love numbers.

        Returns
        -------
        tidal_heating : np.ndarray
            Tidal heating [W]
            This may be restricted to a specific layer or for an entire planet.
        dUdM : np.ndarray
            Tidal potential derivative with respect to the mean anomaly [J kg-1 radians-1]
            This may be restricted to a specific layer or for an entire planet.
        dUdw : np.ndarray
            Tidal potential derivative with respect to the pericentre [J kg-1 radians-1]
            This may be restricted to a specific layer or for an entire planet.
        dUdO : np.ndarray
            Tidal potential derivative with respect to the planet's node [J kg-1 radians-1]
            This may be restricted to a specific layer or for an entire planet.

        Raises
        ------
        MissingAttributeError

        See Also
        --------
        TidalPy.tides.Tides.orbit_spin_changed
        """

        super().collapse_modes()

        # Check to see if all the needed state properties are present and then begin calculations
        if self.tidal_terms_by_frequency is not None:
            nonNone_love_number = list()
            nonNone_neg_imk = list()
            nonNone_effective_q = list()
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

            # Clear global parameters
            self._effective_q_by_orderl = dict()
            self._global_negative_imk_by_orderl = dict()
            self._global_love_by_orderl = dict()

            broke_out = False
            for layer, other_inputs in self._tidal_input_getters_by_layer.items():
                if other_inputs is None:
                    # Not a tidal layer
                    love_number_by_layer[layer] = None
                    neg_imk_by_layer[layer] = None
                    tidal_heating_by_layer[layer] = None
                    dUdM_by_layer[layer] = None
                    dUdw_by_layer[layer] = None
                    dUdO_by_layer[layer] = None

                else:
                    # Find getter functions
                    tidal_scale_getter, radius_getter, bulk_density_getter, gravity_surf_getter = other_inputs

                    if self._world_tidal_input_getters is None:
                        # Use layer properties - set above do nothing here
                        pass
                    else:
                        # Use world properties
                        radius_getter, bulk_density_getter, gravity_surf_getter = self._world_tidal_input_getters

                    # Use getters to find key parameters
                    tidal_scale  = tidal_scale_getter()
                    radius       = radius_getter()
                    bulk_density = bulk_density_getter()
                    gravity_surf = gravity_surf_getter()

                    # Pull out variables that change often
                    shear_modulus = layer.shear_modulus
                    complex_compliances_by_frequency_list = layer.rheology.complex_compliances

                    if shear_modulus is None or complex_compliances_by_frequency_list is None:
                        # uh oh
                        broke_out = True
                        break

                    # Mode collapse will parse through tidal order-l and all unique frequencies and calculate global and
                    #    localized dissipation values
                    tidal_heating, dUdM, dUdw, dUdO, love_number_by_orderl, negative_imk_by_orderl, \
                    effective_q_by_orderl = \
                        self.collapse_modes_func(
                            gravity_surf, radius, bulk_density, shear_modulus,
                            tidal_scale, self.tidal_host.mass, self.tidal_susceptibility,
                            complex_compliances_by_frequency_list,
                            self.tidal_terms_by_frequency, self.max_tidal_order_lvl,
                            cpl_ctl_method=False
                            )

                    # These will be summed for global values
                    nonNone_love_number.append(love_number_by_orderl)
                    nonNone_neg_imk.append(negative_imk_by_orderl)
                    nonNone_effective_q.append(effective_q_by_orderl)

                    nonNone_tidal_heating.append(tidal_heating)
                    nonNone_dUdM.append(dUdM)
                    nonNone_dUdw.append(dUdw)
                    nonNone_dUdO.append(dUdO)

                    # Accessed by layer methods for thermal evolution
                    tidal_heating_by_layer[layer] = tidal_heating
                    neg_imk_by_layer[layer] = negative_imk_by_orderl
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

                # Sum the parameters that are stored by order-l
                for order_l in range(2, self.max_tidal_order_lvl + 1):

                    # TODO: Should all of these simply be sums of each layer's contribution?
                    self._effective_q_by_orderl[order_l] = \
                        sum(
                            [nonNone_effective_q_for_layer[order_l] for nonNone_effective_q_for_layer
                             in nonNone_effective_q]
                            )
                    self._global_negative_imk_by_orderl[order_l] = \
                        sum(
                            [nonNone_neg_imk_for_layer[order_l] for nonNone_neg_imk_for_layer
                             in nonNone_neg_imk]
                            )
                    self._global_love_by_orderl[order_l] = \
                        sum(
                            [nonNone_love_number_for_layer[order_l] for nonNone_love_number_for_layer
                             in nonNone_love_number]
                            )

                # Now tell other methods to update now that derivatives and heating has been altered
                self.world.dissipation_changed()

            else:
                # Broke out from the layer loop for some reason. Likely because this function was called when not
                #  everything was set.
                # Let it go.
                pass
                # raise MissingAttributeError

        # Return tidal heating and derivatives
        return self.tidal_heating_global, self.dUdM, self.dUdw, self.dUdO

    # # State properties
    @property
    def tidal_heating_by_layer(self) -> Dict['PhysicalLayerType', 'FloatArray']:
        """ Tidal heating stored by layer instance [W] """
        return self._tidal_heating_by_layer

    @tidal_heating_by_layer.setter
    def tidal_heating_by_layer(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def negative_imk_by_layer_by_orderl(self) -> Dict['PhysicalLayerType', List['FloatArray']]:
        """ -Im[k2] stored by layer instance and by order l"""
        return self._negative_imk_by_layer

    @negative_imk_by_layer_by_orderl.setter
    def negative_imk_by_layer_by_orderl(self, value):
        raise IncorrectMethodToSetStateProperty
