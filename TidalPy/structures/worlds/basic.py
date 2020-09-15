from __future__ import annotations

import os
from typing import TYPE_CHECKING, Tuple

import numpy as np

from .defaults import world_defaults
from .. import PhysicalObjSpherical
from ... import debug_mode, use_disk, tidalpy_loc, configurations, log
from ...exceptions import (ImproperPropertyHandling, UnusualRealValueError,
                           IncorrectAttributeType, AttributeNotSetError, AttributeException,
                           IOException)
from ...stellar.stellar import (equilibrium_insolation_functions, equilibrium_temperature)
from ...utilities.graphics import geotherm_plot
from ...utilities.numpyHelper import reshape_help
from ...utilities.types import FloatArray

planet_config_loc = os.path.join(tidalpy_loc, 'planets', 'planet_configs')

if TYPE_CHECKING:
    from ...orbit import Orbit


FREQ_TO_DAY = (2. * np.pi ) / 86400.

class BaseWorld(PhysicalObjSpherical):

    """ WorldBase Class - Base class used to build other world classes.


    See Also
    --------
    Parent Class:
        TidalPy.structures.PhysicalObjSpherical
    Child Classes:
        TidalPy.structures.world.TidalWorld
        TidalPy.structures.world.GasGiantWorld
        TidalPy.structures.world.StarWorld
        TidalPy.structures.world.LayeredWorld
        TidalPy.structures.world.GasGiantLayeredWorld
        TidalPy.structures.world.BurnManWorld
    """

    default_config = world_defaults
    world_class = 'base'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):
        """ BaseWorld constructor (child of PhysicalObjSpherical)

        Parameters
        ----------
        world_config : dict
            Configuration file used to build the world. User provided configs override default configurations that
                TidalPy assumes.
            Please see files stored in <TidalPy directory>/structures/worldConfigs for example configuration dict.
        name : str = None
            Name of the world. If None, will use name provided in world_config.
        initialize : bool = True
            Determines if initial setup should be performed on the world (loading in data from world_config).
        """

        # Load in defaults
        self.default_config = self.default_config[self.world_class]

        # Key Attributes
        self.name = name
        log.debug(f'Setting up new world: {self.name}; class type = {self.world_class}.')

        super().__init__(config=world_config)

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
        self._is_spin_sync = False

        # Additional orbit information
        # orbit_location is an index for where a planet is at in the orbit's geometry. This is set by the initialization
        # of an orbit class.
        #   0 == StarWorld location
        #   1 == Host location
        #   2+ == Various target bodies
        self.orbit_location = None

        # Models to be initialized later
        self.equilibrium_insolation_func = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False, reinit_geometry: bool = True, set_by_burnman: bool = False):
        """ Initialize or Reinitialize the world based on changes to its configurations.

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
        """

        super().reinit(initial_init=initial_init, set_by_burnman=set_by_burnman)

        # Pull out some basic information from the config
        if self.name is None and 'name' in self.config:
            self.name = self.config['name']
        self._albedo = self.config['albedo']
        self._emissivity = self.config['emissivity']
        self._is_spin_sync = self.config['force_spin_sync']
        self._pressure_above = self.config['surface_pressure']

        # Setup geometry
        if reinit_geometry:
            self.set_geometry(self.config['radius'], self.config['mass'])
            self.set_static_pressure(pressure_above=self.pressure_above, build_slices=True)

        if self.orbit is not None:
            self.update_orbit()

    def clear_state(self, preserve_orbit: bool = False):

        log.debug(f'Clear state called for {self}. Orbit preserved = {preserve_orbit}.')

        super().clear_state()

        # Clear world-specific state properties
        self._spin_freq = None
        self._obliquity = None
        self._time = None
        self._surf_temperature = None
        self._int2surface_heating = None

        # Clear the global shape now that all state data has been cleared
        #     (probably the reason this function was called in the first place)
        self._global_shape = None

        # Setup functions and models
        self.equilibrium_insolation_func = equilibrium_insolation_functions[self.config['equilibrium_insolation_model']]

        # TODO:
        # Reset any orbits this world is connected to
        if not preserve_orbit and self.orbit is not None:
            # Purge orbital properties for this world from the Orbit class.
            self.orbit.set_orbit(self, purge_data=True)


    def set_geometry(self, radius: float, mass: float, thickness: float = None, mass_below: float = 0.,
                     update_state_geometry: bool = True, build_slices: bool = True):
        """ Calculates and sets the world's physical parameters based on user provided input.

        Assumptions
        -----------
        Spherical Geometry

        Parameters
        ----------
        radius : float
            Outer radius of object [m]
        mass : float
            Mass of object [kg]
        thickness : float = None
            Thickness of the object [m]
        mass_below : float = 0.
            Mass below this object (only applicable for shell-like structures)
            Used in gravity and pressure calculations
        update_state_geometry : bool = True
            Update the class' state geometry
        build_slices : bool = True
            If True, method will attempt to calculate gravities, densities, etc. for each slice.

        """

        # Thickness of a world will always be equal to its radius.
        del thickness
        super().set_geometry(radius, mass, thickness=radius, mass_below=0.,
                             update_state_geometry=update_state_geometry, build_slices=build_slices)

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

    def paint(self, depth_plot: bool = False, auto_show: bool = True):
        """ Create a geotherm or depth plot of the planet's gravity, pressure, and density
        Parameters
        ----------
        depth_plot : bool = False
            If true the plot will be versus depth rather than radius.
        auto_show : bool = False
            Calls plt.show() if true.
        Returns
        -------
        figure: matplotlib.pyplot.figure
        """

        figure = geotherm_plot(self.radii, self.gravity_slices, self.pressure_slices, self.density_slices,
                               bulk_density=self.density_bulk, planet_name=self.name,
                               planet_radius=self.radius, depth_plot=depth_plot, auto_show=auto_show)

        return figure

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

        log.debug(f'Saving world configurations for {self}.')

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
        self.save_config(save_to_run_dir=False, additional_save_dirs=save_locales,
                         overwrite=configurations['overwrite_configs'])

    def kill_world(self):
        """ Performs saving tasks when the world is about to be deleted due to end of run

        The exit_planets variable in the main TidalPy configurations controls rather or not this method ever gets
        called automatically.
        """

        log.debug(f'Killing world {self}.')

        # Save configuration file
        if use_disk:
            if configurations['auto_save_planet_config_to_tidalpydir']:
                tidalpy_planet_cfg_dir = [planet_config_loc]
            else:
                tidalpy_planet_cfg_dir = list()
            self.save_config(save_to_run_dir=configurations['auto_save_planet_config_to_rundir'],
                             additional_save_dirs=tidalpy_planet_cfg_dir,
                             overwrite=configurations['overwrite_configs'])


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
        """ The time of either the Orbit or the object (used for radiogenic calculations) [Myr] """
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
        """ Spin (Sidereal Rotation) Frequency of the World [rad s-1] """
        return self._spin_freq

    @spin_freq.setter
    def spin_freq(self, new_spin_frequency: np.ndarray):
        """ Update Spin (Sidereal Rotation) Frequency of the World

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
        """ World's Albedo """
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
        """ World's Emissivity """
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
        """ World's Surface Temperature [K] """
        # TODO: Alias this with PhysicalObject's self.temperature_outer?
        return self._surf_temperature

    @surface_temperature.setter
    def surface_temperature(self, value):
        raise ImproperPropertyHandling

    @property
    def obliquity(self) -> np.ndarray:
        """ World's Obliquity [rad]

        This obliquity must be relative to the orbital plane defined by the tidal target and tidal host [1]_.
            If the star is neither the host nor target, then it should not be used as a reference object.

        References
        ----------
        .. [1] J. P. Renaud et al, "Tidal Dissipation in Dual-Body, Highly Eccentric, and
           Non-synchronously Rotating Systems: Applications to Pluto-Charon and the Exoplanet TRAPPIST-1e"
           The Planetary Science Journal, vol. 22, pp. TBA, 2020.
        """

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
        """ The Orbit Class Associated with this World """
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
        """ World's Semi-Major Axis (stored in the world's Orbit class) [m] """
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
    def orbital_frequency(self):
        """ World's Orbital Frequency (inverse of orbital period; stored in the world's Orbit class) [rad s-1]"""
        return self.orbit.get_orbital_freq(self.orbit_location)

    @orbital_frequency.setter
    def orbital_frequency(self, new_orbital_frequency: FloatArray):

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
        """ World's Orbital Eccentricity (inverse of orbital period; stored in the world's Orbit class) """
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
        """ World's Surface Heating due to Stellar Insolation [W] """
        if self.orbit is not None:
            return self.orbit.get_insolation(self)
        return None

    @insolation_heating.setter
    def insolation_heating(self, value):
        raise ImproperPropertyHandling


    # # Aliased properties
    @property
    def orbital_freq(self):
        """ Alias of BaseWorld.orbital_frequency """
        return self.orbital_frequency

    @orbital_freq.setter
    def orbital_freq(self, new_orbital_frequency):
        self.orbital_frequency = new_orbital_frequency

    @property
    def n(self):
        """ Alias of BaseWorld.orbital_frequency """
        return self.orbital_frequency

    @n.setter
    def n(self, new_orbital_frequency):
        self.orbital_frequency = new_orbital_frequency

    @property
    def orbital_motion(self):
        """ Alias of BaseWorld.orbital_frequency """
        return self.orbital_frequency

    @orbital_motion.setter
    def orbital_motion(self, value):
        self.orbital_frequency = value

    @property
    def orbital_period(self):
        """ Calculated from BaseWorld.orbital_frequency [days]

        .. math:: ( 2 \pi  / 86400 ) / n
        """

        return FREQ_TO_DAY / self.orbital_frequency

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
        self.orbital_frequency = new_orbital_freq


    # Dunder properties
    def __str__(self):

        name = self.name
        if name is None:
            name = f'Unknown World'

        name = f'{name} ({self.world_class} world)'

        return name

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name}, {self.__class__} object at {hex(id(self))}'
        return f'Unknown World, {self.__class__} object at {hex(id(self))}'
