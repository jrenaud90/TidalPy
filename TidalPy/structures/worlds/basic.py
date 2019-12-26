from __future__ import annotations

import copy
import os
from typing import List, TYPE_CHECKING, Tuple, Union

import burnman
import dill
import numpy as np

from ...types import NoneType
from ...dynamics.modes_l2 import nsr_modes, spin_sync_modes
from ...utilities.numpy_help import value_np_cleanup
from .defaults import world_defaults
from .. import PhysicalObjSpherical
from ... import debug_mode, use_disk, tidalpy_dir, __version__
from ...configurations import (auto_save_planet_config_to_rundir, auto_save_planet_config_to_tidalpydir,
                              auto_save_planet_dill_to_rundir, auto_save_planet_dill_to_tidalpydir, overwrite_configs,
                              overwrite_dills)
from ...dynamics import spin_rate_derivative
from ...exceptions import (ImproperAttributeHandling, ParameterMissingError, ReinitError, UnusualRealValueError,
                           IncompatibleModelError, IncorrectAttributeType)
from ...graphics.planet_plot import geotherm_plot
from ...initialize import log
from ...io import inner_save_dir
from ...physics.stellar import (efftemp_from_luminosity, equilibrium_insolation_functions, equilibrium_temperature,
                               luminosity_from_efftemp, luminosity_from_mass)
from ...utilities.dill_helper import dill_file_path
from ...structures.layers import ThermalLayer

planet_config_loc = os.path.join(tidalpy_dir, 'planets', 'planet_configs')

if TYPE_CHECKING:
    from ...orbit import Orbit


class WorldBase(PhysicalObjSpherical):

    """ WorldBase Class - Base class used to build other world classes.
    """

    default_config = world_defaults
    world_class = 'base'

    def __init__(self, world_config: dict, name: str = None):

        # Load in defaults
        self.default_config = self.default_config[self.world_class]

        super().__init__(config=world_config)

        # Independent State variables
        self.name = name
        self._spin_freq = None
        self._time = None
        self._albedo = None
        self._emissivity = None

        # Orbit reference
        self._orbit = None  # type: Union[Orbit, NoneType]

        # Other flags
        self.is_host = False
        self.is_spin_sync = False

        # Additional orbit information
        # orbit_location is an index for where a planet is at in the orbit's geometry. This is set by the initialization
        # of an orbit class.
        #   0 == Star location
        #   1 == Host location
        #   2+ == Various target bodies
        self.orbit_location = None

        # Models to be initialized later
        self.equilibrium_insolation_func = None

        self.pyname = f'{self.name}_{self.world_class}'

    def reinit(self):
        """ Re-initialize the basic world based on changed to its configuration."""

        # Pull out some basic information from the config
        if self.name is None and 'name' in self.config:
            self.name = self.config['name']
        self._albedo = self.config['albedo']
        self._emissivity = self.config['emissivity']
        self._is_spin_sync = self.config['force_spin_sync']
        self.update_orbit()

    def config_update(self):
        """ Config update changes various state parameters of an object based on the self.config
        """
        super().config_update()

        # Call reinit as some configurations may have changed
        self.reinit()

        # Setup functions and models
        self.equilibrium_insolation_func = equilibrium_insolation_functions[self.config['equilibrium_insolation_model']]

    def set_geometry(self, radius: float, mass: float):

        # Thickness of a world will always be equal to its radius.
        super().set_geometry(radius, mass, thickness=radius)

    def update_orbit(self):
        """ Performs state updates whenever the planet's orbital parameters are changed
        """

        # Tell the orbit class to update this planet's insolation heating
        if self.orbit is not None:
            self.orbit.calculate_insolation(self.orbit_location)

    def set_state(self, orbital_freq: np.ndarray = None, semi_major_axis: np.ndarray = None,
                  eccentricity: np.ndarray = None, inclination: np.ndarray = None,
                  spin_freq: np.ndarray = None, time: np.ndarray = None):
        """ Set multiple orbital parameters at once, this reduces the number of calls to self.orbit_change

        This contains a wrapper to the orbit class method set_state. It extends that function by including the spin
        frequency.

        This function has better performance than the individual setters for these parameters if (and only if) you are
        changing two or more of them at the same time.

        Parameters
        ----------
        orbital_freq : np.ndarray
            Mean orbital frequency of an object in [rads s-1]
            Optional, exclusive w/ semi_major_axis
        semi_major_axis : np.ndarray
            Semi-major axis of an object in [m]
            Optional, exclusive w/ semi_major_axis
        eccentricity : np.ndarray
            Orbital eccentricity
            Optional
        inclination : np.ndarray
            Orbital inclination in [rads]
            Optional
        spin_freq : np.ndarray
            Object's rotation frequency in [rads s-1]
            Optional
        time : np.ndarray
            Object's calculation time (used in things like radiogenics) in [s]
            Optional
        """

        if time is not None:
            self.time = time
        if spin_freq is not None:
            self.spin_freq = spin_freq

        self.orbit.set_orbit(self.orbit_location, orbital_freq, semi_major_axis, eccentricity, inclination,
                             set_by_planet=True)
        self.update_orbit()

    def save_world(self):

        #TODO
        pass

    #TODO: new way to kill world that doesn't throw so many errors?

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

        # Dill this object and save it
        if use_disk:
            dill_name, tidalpy_dill_path = dill_file_path(self.name)
            if auto_save_planet_dill_to_tidalpydir:
                if os.path.isfile(tidalpy_dill_path) and not overwrite_dills:
                    pass
                else:
                    with open(tidalpy_dill_path, 'wb') as dill_file:
                        dill.dump(self, dill_file)
            if auto_save_planet_dill_to_rundir:
                rundir_dill_path = os.path.join(inner_save_dir, dill_name)
                if os.path.isfile(rundir_dill_path) and not overwrite_dills:
                    pass
                else:
                    with open(rundir_dill_path, 'wb') as dill_file:
                        dill.dump(self, dill_file)


    @property
    def time(self) -> np.ndarray:
        return self._time

    @time.setter
    def time(self, new_time: np.ndarray):
        """ Update spin frequency of a world

        Parameters
        ----------
        new_time: np.ndarray
            New time [Myr]
        """

        new_time = value_np_cleanup(new_time)
        if debug_mode:
            if np.any(abs(new_time)) > 1.e6:
                raise UnusualRealValueError(f'Time should be entered in units of [Myr]. '
                                            f'|{new_time}| seems very large.')

        self._time = new_time

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
        new_spin_frequency = value_np_cleanup(new_spin_frequency)
        if debug_mode:
            if np.any(np.abs(new_spin_frequency)) > 1.e-3:
                raise UnusualRealValueError(f'Spin-frequency should be entered in units of [rads s-1]. '
                                            f'|{new_spin_frequency}| seems very large.')

        self._spin_freq = new_spin_frequency
        if not self.is_spin_sync and self.orbit is not None:
            # A change to the spin-rate during NSR will change tidal modes
            self.update_orbit()

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
            self.orbit.calculate_insolation(self)

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
            self.orbit.calculate_insolation(self)

    # Class properties
    @property
    def orbit(self) -> Orbit:
        return self._orbit

    @orbit.setter
    def orbit(self, value):
        raise ImproperAttributeHandling('The orbit is set once the planet is passed to an Orbit class constructor.')

    # Outerscope property references
    @property
    def semi_major_axis(self) -> np.ndarray:
        return self.orbit.get_semi_major_axis(self.orbit_location)

    @semi_major_axis.setter
    def semi_major_axis(self, value: np.ndarray):
        self.orbit.set_semi_major_axis(self.orbit_location, value, set_by_planet=True)
        self.update_orbit()

    @property
    def orbital_freq(self) -> np.ndarray:
        return self.orbit.get_orbital_freq(self.orbit_location)

    @orbital_freq.setter
    def orbital_freq(self, value: np.ndarray):

        self.orbit.set_orbital_freq(self.orbit_location, value, set_by_planet=True)
        self.update_orbit()

    @property
    def eccentricity(self) -> np.ndarray:
        return self.orbit.get_eccentricity(self.orbit_location)

    @eccentricity.setter
    def eccentricity(self, value: np.ndarray):
        self.orbit.set_eccentricity(self.orbit_location, value, set_by_planet=True)
        self.update_orbit()

    @property
    def inclination(self) -> np.ndarray:
        return self.orbit.get_inclination(self.orbit_location)

    @inclination.setter
    def inclination(self, value: np.ndarray):
        self.orbit.set_inclination(self.orbit_location, value, set_by_planet=True)
        self.update_orbit()

    @property
    def insolation_heating(self) -> np.ndarray:

        inso_heating = None
        if self.orbit is not None:
            inso_heating = self.orbit.get_insolation(self)

        return inso_heating

    @insolation_heating.setter
    def insolation_heating(self, value):
        raise ImproperAttributeHandling

    @property
    def surface_temperature(self) -> np.ndarray:

        surf_temp = None
        if self.orbit is not None:
            surf_temp = self.orbit.get_surface_temperature(self)

        return surf_temp

    @surface_temperature.setter
    def surface_temperature(self, value):
        raise ImproperAttributeHandling

    # Aliased Parameter Names
    @property
    def orbital_motion(self):
        return self.orbital_freq

    @orbital_motion.setter
    def orbital_motion(self, value):
        self.orbital_freq = value

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

    def __init__(self, planet_config: dict, name: str = None):

        super().__init__(planet_config, name=name)

        self.set_geometry(planet_config['radius'], planet_config['mass'])

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