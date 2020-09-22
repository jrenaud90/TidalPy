from typing import Tuple

import burnman

from .layered import LayeredWorld
from ...exceptions import InitiatedPropertyChangeError


class BurnManWorld(LayeredWorld):

    """ BurnManWorld
    Tidal world_types that have layers whose state properties are initialized by user-provided equations of state. These
        calculations are done by the Burnman package.


    See Also
    --------
    Parent Class:
        TidalPy.structures.world_types.LayeredWorld
    """

    world_class = 'burnman'

    def __init__(self, world_config: dict, burnman_world: burnman.Planet, burnman_layers: Tuple[burnman.Layer, ...],
                 name: str = None, initialize: bool = True):
        """ BurnManWorld constructor

        Parameters
        ----------
        world_config : dict
            Configuration file used to build the world. User provided configs override default configurations that
                TidalPy assumes.
            Please see files stored in <TidalPy directory>/structures/world_configs for example configuration dict.
        burnman_world : burnman.Planet
            An initialized burnman planet object.
        burnman_layers : Tuple[burnman.Layer, ...]
            List of initialized burnman layer objects.
        name : str = None
            Name of the world. If None, will use name provided in world_config.
        initialize : bool = True
            Determines if initial reinit should be performed on the world (loading in data from world_config).
        """

        # Store BurnMan information in state properties
        self._bm_world = burnman_world
        self._bm_layers = burnman_layers

        # Setup TidalPy layers (part of parent class)
        super().__init__(world_config, name, initialize=False)

        # Initialize world
        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False, reinit_geometry: bool = True, setup_simple_tides: bool = False,
               set_by_burnman: bool = True, reinit_layers: bool = True):
        """ Initialize or Reinitialize the world based on changes to its configurations.

        This must be called at least once before an instance can be used. The constructor will automatically make an
            initial call to reinit unless told to not to.

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this function has been called.
        reinit_geometry : bool = True
            If `True`, the initializer will automatically call the `set_geometry()` method.
        set_by_burnman : bool = True
            Set to `True` if called from a burnman world.
        setup_simple_tides : bool = True
            Set to `True` if a global CPL/CTL tidal calculation is desired.
        reinit_layers : bool = True
            If `True`, calls to the world's layers' reinit() method.
        """

        # Pull out burnman mass and radius - use these for the TidalPy planet as well
        self._mass = self.bm_world.mass
        self._radius = self.bm_world.radius_planet
        self._moi = self.bm_world.moment_of_inertia

        # Setup Geometry
        if reinit_geometry:
            self.set_geometry(self.radius, self.mass)

        # Make call to parent reinit
        super().reinit(initial_init, reinit_geometry, setup_simple_tides, set_by_burnman, reinit_layers)


    # # Initiated properties
    @property
    def bm_world(self):
        return self._bm_world

    @bm_world.setter
    def bm_world(self, value):
        raise InitiatedPropertyChangeError

    @property
    def bm_layers(self):
        return self._bm_layers

    @bm_layers.setter
    def bm_layers(self, value):
        raise InitiatedPropertyChangeError
