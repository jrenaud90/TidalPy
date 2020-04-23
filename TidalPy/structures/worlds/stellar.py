from .tidal import TidalWorld


# TODO: Implement a fixed-q tides class/method for stellar and gas planets. Wait it is a tidal world...

class StarWorld(TidalWorld):

    world_class = 'star'

    def __init__(self, star_config: dict, name: str = None, initialize: bool = True):

        super().__init__(star_config, name=name, initialize=False)

        # StarWorld-specific attributes
        self.luminosity = None
        self.effective_temperature = None

        if initialize:
            self.reinit()

    def reinit(self):

        super().reinit()

        self.luminosity = self.config.get('luminosity', None)
        self.effective_temperature = self.config.get('effective_temperature', None)

        if self.luminosity is None:
            # If no luminosity provided: Try to convert effective surface temperature
            if self.effective_temperature is None:
                # if that fails, try to estimate from mass
                log(
                        'Luminosity and effective temperature of {self.name} was not provided. Estimating from stellar mass.')
                self.luminosity = luminosity_from_mass(self.mass)
                self.effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)
            else:
                self.luminosity = luminosity_from_efftemp(self.effective_temperature, self.radius)
        else:
            if self.effective_temperature is None:
                self.effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)

    def __str__(self):

        name = self.name
        if name is None:
            name = 'Unknown'
        else:
            name = name.title()
        return name

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name} {self.__class__} object at {hex(id(self))}'
        return f'{self.__class__} object at {hex(id(self))}'