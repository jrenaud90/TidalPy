from __future__ import annotations

import os
from typing import TYPE_CHECKING, Tuple

import numpy as np
from ...rheology.complexCompliance.compliance_models import fixed_q, fixed_q_array
from TidalPy.utilities.types import FloatArray

from ...utilities.numpyHelper import reshape_help
from .defaults import world_defaults
from .. import PhysicalObjSpherical
from ... import debug_mode, use_disk, tidalpy_dir
from ...configurations import (auto_save_planet_config_to_rundir, auto_save_planet_config_to_tidalpydir,
                               overwrite_configs)
from ...exceptions import (ImproperPropertyHandling, UnusualRealValueError,
                           IncorrectAttributeType, AttributeNotSetError, AttributeException,
                           IOException)
from ...initialize import log
from ...stellar.stellar import (equilibrium_insolation_functions, equilibrium_temperature)

planet_config_loc = os.path.join(tidalpy_dir, 'planets', 'planet_configs')

if TYPE_CHECKING:
    from ...orbit import Orbit


FREQ_TO_DAY = (2. * np.pi ) / 86400.

class WorldBase(PhysicalObjSpherical):

    """ WorldBase Class - Base class used to build other world classes.
    """

    default_config = world_defaults
    world_class = 'base'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):

        # Load in defaults
        self.default_config = self.default_config[self.world_class]

        super().__init__(config=world_config)

        # Attributes
        self.name = name

        # Independent State variables
        self._spin_freq = None
        self._obliquity = None
        self._time = None
        self._albedo = None
        self._emissivity = None
        self._surf_temperature = None
        self._int2surface_heating = None
        # Global shape is overridden if an orbit is present (see property definition)
        self._global_shape = None

        # Orbit reference
        self._orbit = None

        # Other flags
        self.is_host = False
        self.is_spin_sync = False

        # Additional orbit information
        # orbit_location is an index for where a planet is at in the orbit's geometry. This is set by the initialization
        # of an orbit class.
        #   0 == StarWorld location
        #   1 == Host location
        #   2+ == Various target bodies
        self.orbit_location = None

        # Models to be initialized later
        self.equilibrium_insolation_func = None

        self.pyname = f'{self.name}_{self.world_class}'

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        """ Re-initialize the basic world based on changed to its configuration."""

        super().reinit()

        if not initial_init:
            self.clear_state()

        # Pull out some basic information from the config
        if self.name is None and 'name' in self.config:
            self.name = self.config['name']
        self._albedo = self.config['albedo']
        self._emissivity = self.config['emissivity']
        self.is_spin_sync = self.config['force_spin_sync']

        if self.orbit is not None:
            self.update_orbit()

    def clear_state(self, preserve_orbit: bool = False):

        super().clear_state()

        # Clear world-specific state properties
        self._spin_freq = None
        self._obliquity = None
        self._time = None
        self._surf_temperature = None
        self._int2surface_heating = None

        if not preserve_orbit and self.orbit is not None:
            # Purge orbital properties for this world from the Orbit class.
            self.orbit.set_orbit(self, purge_data=True)

        # Clear the global shape now that all state data has been cleared
        #     (probably the reason this function was called in the first place)
        self._global_shape = None

    def config_update(self):
        """ Config update changes various state parameters of an object based on the self.config
        """
        super().config_update()

        # Call reinit as some configurations may have changed
        self.reinit()

        # Setup functions and models
        self.equilibrium_insolation_func = equilibrium_insolation_functions[self.config['equilibrium_insolation_model']]

    def set_geometry(self, radius: float, mass: float, thickness: float = None):

        # Thickness of a world will always be equal to its radius.
        del thickness
        super().set_geometry(radius, mass, thickness=radius)

    def change_shape(self, new_shape: Tuple[int, ...], force_override: bool = False):
        """ Update the global shape parameter

        Global shape parameter is used to set the shape of arrays throughout the planet.

        It is overridden by the orbit's global shape, if an Orbit class has been applied.

        Parameters
        ----------
        new_shape : Tuple[int, ...]
            New array shape provided by numpy.ndarray.shape
        force_override : bool = False
            Force an override if global_shape is already set.
        """

        if self.orbit is not None:
            self.orbit.change_shape(new_shape, set_by_planet=self)

        if self._global_shape is None:
            self._global_shape = new_shape
        else:
            if force_override:
                self._global_shape = new_shape
            raise AttributeException('global_shape is already set.')

    def update_surface_temperature(self):
        """ Update the planet's surface temperature based on its distance from its host star.

        Orbit must be applied before surface temperature can be updated.
        """

        if self.orbit is not None:
            surface_heating = self.insolation_heating
            if self.int2surface_heating is not None:
                surface_heating += self.int2surface_heating

            self._surf_temperature = equilibrium_temperature(surface_heating, self.radius, self.emissivity)

    def update_orbit(self):
        """ Performs state updates whenever the planet's orbital parameters are changed
        """

        # Tell the orbit class to update this planet's insolation heating
        if self.orbit is not None:
            self.orbit.calculate_insolation(self.orbit_location)

        # A change to the orbit will also change tides
        self.update_tides()

        # A change to semi-major axis for a planet orbiting a star will cause a change to its surface temperature
        self.update_surface_temperature()

    def update_tides(self):
        """ Update the tides module

        Only applicable for some child class types
        """

        # Not applicable for this class. Do nothing
        pass

    def set_state(self, orbital_freq: FloatArray = None, semi_major_axis: FloatArray = None,
                  eccentricity: FloatArray = None, obliquity: FloatArray = None,
                  spin_freq: FloatArray = None, time: FloatArray = None):
        """ Set multiple orbital parameters at once, this reduces the number of calls to self.orbit_change

        This contains a wrapper to the orbit class method set_state. It extends that function by including the spin
        frequency.

        This function has better performance than the individual setters for these parameters if (and only if) you are
        changing two or more of them at the same time.

        Parameters
        ----------
        orbital_freq : FloatArray (Optional)
            Mean orbital frequency of an object in [rads s-1]
            Exclusive w/ semi_major_axis
        semi_major_axis : FloatArray (Optional)
            Semi-major axis of an object in [m]
            Exclusive w/ orbital_freq
        eccentricity : FloatArray (Optional)
            Orbital eccentricity
        obliquity : FloatArray (Optional)
            Planet's axial tilt (or inclination) in [rads]
        spin_freq : FloatArray (Optional)
            Object's rotation frequency in [rads s-1]
        time : FloatArray (Optional)
            Object's calculation time (used in things like radiogenics) in [s]
        """

        update_tides_flag = False

        for planet_param_name, planet_param in {'time': time, '_spin_freq': spin_freq, '_obliquity': obliquity}.items():

            if planet_param is None:
                continue

            if planet_param_name in ['_spin_freq', '_obliquity']:
                update_tides_flag = True

            new_shape, planet_param = \
                reshape_help(planet_param, self.global_shape, call_locale=f'{self}.set_state.{planet_param_name}')
            if new_shape:
                self.change_shape(new_shape)

            setattr(self, planet_param_name, planet_param)

        if orbital_freq is not None or semi_major_axis is not None or eccentricity is not None:
            if self.orbit is not None:
                raise AttributeNotSetError('Trying to set orbital parameters at the planet '
                                           'level before an Orbit class has been applied.')

            self.orbit.set_orbit(self.orbit_location, orbital_freq, semi_major_axis, eccentricity, set_by_planet=True)

            # A call to update_orbit will automatically call update tides.
            #     So we do not need to worry about another call to it.
            self.update_orbit()
        else:
            # Update tides has not been called at this point. Call it if needed
            if update_tides_flag:
                self.update_tides()

    def save_world(self, save_dir: str = None, no_cwd: bool = False, save_to_tidalpy_dir: bool = False):
        """ Save the world's configuration file to a specified directory.

        Parameters
        ----------
        save_dir : str = None
            The current working directory will be preappended unless no_cwd is set to true.
            If no directory is provided it will be saved os.getcwd()
        no_cwd : bool = False
            If True, the current working directory will *not* be prepended to the provided directory.
        save_to_tidalpy_dir : bool = False
            If True, the config will be saved to the TidalPy directory as well as the CWD.
        """

        if not use_disk:
            raise IOException("User attempted to save a world's configurations when TidalPy has been set to "
                              "not use_disk")
        log(f'Saving world {self.name}...', level='debug')
        save_locales = list()
        if save_to_tidalpy_dir:
            save_locales.append(planet_config_loc)

        if save_dir is not None:
            if no_cwd:
                save_locales.append(save_dir)
            else:
                save_locales.append(os.path.join(os.getcwd(), save_dir))
        else:
            save_locales.append(os.getcwd())

        # No need to save to run dir as that will automatically happen when the planet is killed.
        self.save_config(save_to_run_dir=False, additional_save_dirs=save_locales, overwrite=overwrite_configs)

    def kill_world(self):
        """ Performs saving tasks when the world is about to be deleted due to end of run

        The exit_planets variable in the main TidalPy configurations controls rather or not this method ever gets
        called automatically.
        """

        log(f'Killing world {self.name}...', level='debug')
        # Save configuration file
        if use_disk:
            if auto_save_planet_config_to_tidalpydir:
                tidalpy_planet_cfg_dir = [planet_config_loc]
            else:
                tidalpy_planet_cfg_dir = list()
            self.save_config(save_to_run_dir=auto_save_planet_config_to_rundir,
                             additional_save_dirs=tidalpy_planet_cfg_dir, overwrite=overwrite_configs)


    # State properties
    @property
    def global_shape(self) -> Tuple[int, ...]:
        if self.orbit is None:
            return self._global_shape
        else:
            return self.orbit.global_shape

    @global_shape.setter
    def global_shape(self, value):
        raise ImproperPropertyHandling

    @property
    def time(self) -> np.ndarray:
        if self.orbit is None:
            return self._time
        else:
            return self.orbit.time

    @time.setter
    def time(self, new_time: FloatArray):
        """ Update time for a world

        Parameters
        ----------
        new_time : FloatArray
            Time [Myr]
        """
        if self.orbit is None:
            new_shape, new_time = reshape_help(new_time, self.global_shape, call_locale=f'{self}.time.setter')
            if new_shape:
                self.change_shape(new_shape)

            if debug_mode:
                if np.any(abs(new_time)) > 1.e6:
                    raise UnusualRealValueError(f'Time should be entered in units of [Myr]. '
                                                f'|{new_time}| seems very large.')

            self._time = new_time
        else:
            raise ImproperPropertyHandling('Time must be set at the Orbit-level once an orbit is applied.')

    @property
    def spin_freq(self) -> np.ndarray:
        return self._spin_freq

    @spin_freq.setter
    def spin_freq(self, new_spin_frequency: np.ndarray):
        """ Update spin frequency of a world

        Parameters
        ----------
        new_spin_frequency: np.ndarray
            New spin frequency (inverse of period) [rads s-1]
        """
        new_shape, new_spin_frequency = reshape_help(new_spin_frequency, self.global_shape,
                                                     call_locale=f'{self}.spin_freq.setter')
        if new_shape:
            self.change_shape(new_shape)

        if debug_mode:
            if np.any(np.abs(new_spin_frequency)) > 1.e-3:
                raise UnusualRealValueError(f'Spin-frequency should be entered in units of [rads s-1]. '
                                            f'|{new_spin_frequency}| seems very large.')

        self._spin_freq = new_spin_frequency
        if not self.is_spin_sync:
            # A change to the spin-rate during NSR will change tidal modes
            self.update_tides()

    @property
    def albedo(self) -> float:
        return self._albedo

    @albedo.setter
    def albedo(self, value: float):

        if debug_mode:
            if type(value) != float:
                raise IncorrectAttributeType

        self._albedo = value

        if self.orbit is not None:
            # This will change the planet's surface temperature
            self.update_surface_temperature()

    @property
    def emissivity(self) -> float:
        return self._emissivity

    @emissivity.setter
    def emissivity(self, value: float):

        if debug_mode:
            if type(value) != float:
                raise IncorrectAttributeType

        self._emissivity = value

        if self.orbit is not None:
            # This will change the planet's surface temperature
            self.update_surface_temperature()

    @property
    def surface_temperature(self) -> np.ndarray:

        return self._surf_temperature

    @surface_temperature.setter
    def surface_temperature(self, value):
        raise ImproperPropertyHandling

    @property
    def obliquity(self) -> np.ndarray:
        return self._obliquity

    @obliquity.setter
    def obliquity(self, new_obliquity: FloatArray):
        new_shape, new_obliquity = reshape_help(new_obliquity, self.global_shape,
                                                call_locale=f'{self}.obliquity.setter')
        if new_shape:
            self.change_shape(new_shape)

        if debug_mode:
            if np.any(new_obliquity > 7.):
                raise UnusualRealValueError('Obliquity should be entered in radians. '
                                            f'A value of {np.max(new_obliquity)} seems unusually large.')
        self._obliquity = new_obliquity
        self.update_tides()

    @property
    def int2surface_heating(self) -> np.ndarray:
        return self._int2surface_heating

    @int2surface_heating.setter
    def int2surface_heating(self, value):
        raise ImproperPropertyHandling

    @property
    def orbit(self) -> Orbit:
        return self._orbit

    @orbit.setter
    def orbit(self, value: Orbit):

        self._orbit = value

        # Now that orbit is set, various other updates can be performed.
        self.update_orbit()


    # Outer-scope properties
    # # Orbit Class
    @property
    def semi_major_axis(self):
        return self.orbit.get_semi_major_axis(self.orbit_location)

    @semi_major_axis.setter
    def semi_major_axis(self, new_semi_major_axis: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError('Can not set semi-major axis until an Orbit class has been applied to the planet.')

        new_shape, new_semi_major_axis = reshape_help(new_semi_major_axis, self.global_shape,
                                                      call_locale=f'{self}.semi_major_axis.setter' )
        if new_shape:
            self.change_shape(new_shape)

        self.orbit.set_semi_major_axis(self.orbit_location, new_semi_major_axis, set_by_planet=True)
        self.update_orbit()

    @property
    def orbital_freq(self):
        return self.orbit.get_orbital_freq(self.orbit_location)

    @orbital_freq.setter
    def orbital_freq(self, new_orbital_frequency: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError('Can not set orbital frequency (or period) until an Orbit class has been applied to the planet.')

        new_shape, new_orbital_frequency = reshape_help(new_orbital_frequency, self.global_shape,
                                                        call_locale=f'{self}.orbital_frequency.setter')
        if new_shape:
            self.change_shape(new_shape)

        self.orbit.set_orbital_freq(self.orbit_location, new_orbital_frequency, set_by_planet=True)
        self.update_orbit()

    @property
    def eccentricity(self):
        return self.orbit.get_eccentricity(self.orbit_location)

    @eccentricity.setter
    def eccentricity(self, new_eccentricity: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError('Can not set eccentricity until an Orbit class has been applied to the planet.')

        new_shape, new_eccentricity = reshape_help(new_eccentricity, self.global_shape,
                                                   call_locale=f'{self}.eccentricity.setter')
        if new_shape:
            self.change_shape(new_shape)

        self.orbit.set_eccentricity(self.orbit_location, new_eccentricity, set_by_planet=True)
        self.update_orbit()

    @property
    def insolation_heating(self):

        if self.orbit is not None:
            return self.orbit.get_insolation(self)
        return None

    @insolation_heating.setter
    def insolation_heating(self, value):
        raise ImproperPropertyHandling


    # Aliased properties
    @property
    def orbital_motion(self):
        return self.orbital_freq

    @orbital_motion.setter
    def orbital_motion(self, value):
        self.orbital_freq = value

    @property
    def orbital_period(self):
        return FREQ_TO_DAY / self.orbital_freq

    @orbital_period.setter
    def orbital_period(self, new_orbital_period: FloatArray):
        """ Set orbital period in days

        Parameters
        ----------
        new_orbital_period : FloatArray
            Orbital period [days]
        """

        if debug_mode:
            unusual_value = False
            if type(new_orbital_period) is np.ndarray:
                if np.any(new_orbital_period < 0.01):
                    unusual_value = True
            else:
                if new_orbital_period < 0.01:
                    unusual_value = True
            if unusual_value:
                raise UnusualRealValueError(f'Unusually small orbital period encountered (should be entered in days): {new_orbital_period}')

        new_orbital_freq =  FREQ_TO_DAY / new_orbital_period
        self.orbital_freq = new_orbital_freq


    # Dunder properties
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


class GeometricWorld(WorldBase):

    """ GeometricWorld Class - Most basic type of world that has a geometry (mass, radius, etc.).
    """

    world_class = 'geometric'

    def __init__(self, planet_config: dict, name: str = None, initialize: bool = True):

        super().__init__(planet_config, name=name, initialize=False)

        self.set_geometry(planet_config['radius'], planet_config['mass'])

        if initialize:
            self.reinit()










class TidalWorld(GeometricWorld):

    """ TidalWorld - Provides a simple base to build tidally dissipative worlds off of.
    """

    world_class = 'simple_tide'

    def __init__(self, planet_config: dict, name: str = None, initialize: bool = True):

        super().__init__(planet_config, name=name, initialize=False)

        # State Properties
        self._fixed_q = None

        self.fixed_q_func = fixed_q
        self.fixed_q_func_array = fixed_q_array()

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):

        super().reinit(initial_init)

        self._fixed_q = self.config['quality_factor']


    # State properties
    @property
    def fixed_q(self) -> float:
        return self._fixed_q

    @fixed_q.setter
    def fixed_q(self, new_fixed_q: float):

        if type(new_fixed_q) is not float:
            raise IncorrectAttributeType

        self._fixed_q = new_fixed_q
        self.update_tides()

