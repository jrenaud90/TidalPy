from .tidal import TidalWorld
from ... import log
from ...stellar.stellar import efftemp_from_luminosity, luminosity_from_mass, luminosity_from_efftemp


# TODO: Implement a fixed-q tides class/method for stellar and gas planets. Wait it is a tidal world...

class StarWorld(TidalWorld):

    world_class = 'star'

    def __init__(self, star_config: dict, name: str = None, initialize: bool = True):

        super().__init__(star_config, name=name, initialize=False)

        # StarWorld-specific attributes
        self.luminosity = None
        self.effective_temperature = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False, reinit_geometry: bool = True, set_by_burnman: bool = False,
               setup_simple_tides: bool = True):
        """ Initialize or Reinitialize the star based on changes to its configurations.

        This must be called at least once before an instance can be used. The constructor will automatically make an
            initial call to reinit unless told to not to.

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this function has been called.
        reinit_geometry : bool = True
            If `True`, the initializer will automatically call the `set_geometry()` method.
        set_by_burnman : bool = False
            Set to `True` if called from a burnman world.
        setup_simple_tides : bool = True
            Set to `True` if a global CPL/CTL tidal calculation is desired.
        """

        super().reinit(initial_init=initial_init, reinit_geometry=reinit_geometry,
                       set_by_burnman=set_by_burnman, setup_simple_tides=setup_simple_tides)

        self.luminosity = self.config.get('luminosity', None)
        self.effective_temperature = self.config.get('effective_temperature', None)

        if self.luminosity is None:
            # If no luminosity provided: Try to convert effective surface temperature
            if self.effective_temperature is None:
                # if that fails, try to estimate from mass
                log.info('Luminosity and effective temperature of {self.name} was not provided. '
                         'Estimating from stellar mass.')
                self.luminosity = luminosity_from_mass(self.mass)
                self.effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)
            else:
                self.luminosity = luminosity_from_efftemp(self.effective_temperature, self.radius)
        else:
            if self.effective_temperature is None:
                self.effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)
