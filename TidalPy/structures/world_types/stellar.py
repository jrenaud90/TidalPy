from .tidal import TidalWorld
from ... import log
from ...exceptions import ConfigPropertyChangeError
from ...stellar.stellar import efftemp_from_luminosity, luminosity_from_mass, luminosity_from_efftemp


# TODO: Implement a fixed-q tides class/method for stellar and gas planets. Wait it is a tidal world...

class StarWorld(TidalWorld):

    """ StarWorld
    Stars can dissipate tidal energy through the CPL/CTL method (or not at all) and have methods to calculate
        insolation heating.


    See Also
    --------
    Parent Class:
        TidalPy.structures.world_types.TidalWorld
    """

    world_class = 'star'

    def __init__(self, star_config: dict, name: str = None, initialize: bool = True):

        super().__init__(star_config, name=name, initialize=False)

        # Configuration properties
        self._luminosity = None
        self._effective_temperature = None

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

        self._luminosity = self.config.get('luminosity', None)
        self._effective_temperature = self.config.get('effective_temperature', None)

        if self.luminosity is None:
            # If no luminosity provided: Try to convert effective surface temperature
            if self.effective_temperature is None:
                # if that fails, try to estimate from mass
                log.info(f'Luminosity and effective temperature of {self} was not provided. '
                         'Estimating these values from the stellar mass.')
                self._luminosity = luminosity_from_mass(self.mass)
                self._effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)
            else:
                self._luminosity = luminosity_from_efftemp(self.effective_temperature, self.radius)
        else:
            if self._effective_temperature is None:
                self._effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)


    # # Configuration properties
    @property
    def luminosity(self) -> float:
        """ Stars luminosity [W] """
        return self._luminosity

    @luminosity.setter
    def luminosity(self, value):
        raise ConfigPropertyChangeError

    @property
    def effective_temperature(self) -> float:
        """ Star's effective surface temperature [K] """
        return self._effective_temperature

    @effective_temperature.setter
    def effective_temperature(self, value):
        raise ConfigPropertyChangeError
