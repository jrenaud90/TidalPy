import atexit
import copy
from typing import TextIO, Union

import dill
import json5
import numpy as np

from .dilled_planets import dilled_planets_loc
from .planet_configs import planet_config_loc
from .. import other_data_locs
from ..burnman_interface.build import build_planet as build_bm_planet
from ..configurations import exit_planets
from ..exceptions import MissingArgumentError, UnknownModelError, ImplementationException
from ..initialize import log
from ..structures.worlds import TidalWorld, world_types
from ..utilities.pathing import get_all_files_of_type


def check_for_duplicates(dict_to_check: dict):
    potential_dups = list()
    for planet_name, planet_filepath in dict_to_check.items():
        possible_dup_index = planet_name.split('_')[-1]
        try:
            int(possible_dup_index)
        except ValueError:
            continue
        else:
            if planet_name not in potential_dups:
                potential_dups.append(planet_name)

    for planet_name in potential_dups:
        log.warn(f'Possible duplicate saved planet found: {planet_name}. Ensure you use the correct subscript.')


# Locate all dilled planets
known_planets_dill = dict()
for potential_data_location in other_data_locs + [dilled_planets_loc]:
    known_planets_dill = {**known_planets_dill, **get_all_files_of_type(potential_data_location, ['dill', 'pickle'])}
check_for_duplicates(known_planets_dill)

# Locate all planet configurations
known_planets_cfg = get_all_files_of_type(planet_config_loc, ['cfg', 'json', 'json5'])
check_for_duplicates(known_planets_cfg)
_configs_are_dict = False


def _cfgpath_to_json():
    global known_planets_cfg
    global _configs_are_dict
    for planet_name, config_path in known_planets_cfg.items():
        with open(config_path, 'r') as config_file:
            known_planets_cfg[planet_name] = json5.load(config_file)

    _configs_are_dict = True


# Check for conflict nam

def build_planet(planet_name: str, planet_config: Union[dict, TextIO] = None, force_build: bool = True,
                 planet_dill_path: str = None):
    log(f'Preparing to find and/or build world: {planet_name}')

    if not force_build:
        # TODO: Once this is implemented and tested change force_build to default=False
        raise NotImplementedError('TidalPy 0.1.0 has not fully implemented dill/pickle loading. '
                                  'You must force_build planets for now.')

    # If planet_config is a file then load it through json and get a dict
    if planet_config is not None:
        if type(planet_config) != dict:
            planet_config = json5.load(planet_config)

        # Make a copy of the dict so any subsequent changes do not affect the original
        planet_config = copy.deepcopy(planet_config)

    # See if dilled planet exists
    need_to_build = force_build
    if not force_build:

        if planet_dill_path is None:
            if planet_name in known_planets_dill:
                planet_dill_path = known_planets_dill[planet_name]
            elif planet_name.lower() in known_planets_dill:
                planet_dill_path = known_planets_dill[planet_name.lower()]

        if planet_dill_path is not None:
            log(f'Dilled planet was found! This will save a lot of time...', level='debug')
            # Planet was found! This will save a lot of time.
            with open(planet_dill_path, 'r') as planet_file:
                planet = dill.load(planet_file)

            # If no new planet_config is provided then we are done
            if planet_config is None:
                log(f'No new configurations were provided. Returning dilled planet', level='debug')
                return planet
            else:
                log(f'New configurations were provided. Attempting to load them into dilled planet.', level='debug')
                planet.replacement_config = planet_config
                planet.reinit()
                log(f'New configurations were successful loaded with no obvious issues.', level='debug')
        else:
            log(f'Dilled version of {planet_name} was not found. Attempting to build.', level='debug')
            need_to_build = True

    planet = None
    if need_to_build:
        if planet_config is None:
            log(f'No manual configuration dictionary was provided for {planet_name}, '
                f'attempting to locate saved configuration file.', level='debug')
            if not _configs_are_dict:
                _cfgpath_to_json()
            if planet_name in known_planets_cfg:
                planet_config = known_planets_cfg[planet_name]
            elif planet_name.lower() in known_planets_cfg:
                planet_config = known_planets_cfg[planet_name.lower()]
            else:
                raise MissingArgumentError('No Planet configuration found or provided. '
                                           'Not enough information to build planet.')
            log(f'Configuration file found', level='debug')

        # Determine world type
        planet_type = planet_config['type']
        if planet_type not in world_types:
            raise UnknownModelError('Unknown world type encountered.')
        planet_class = world_types[planet_type]

        if planet_class == TidalWorld:
            log('Burnman planet type detected. Attempting to build BurnMan planet. This may take a while.',
                level='debug')
            # Build BurnMan Planet first
            burnman_layers, burnman_planet = build_bm_planet(planet_config)
            log('Burnman planet build complete!')
            planet = planet_class(planet_config, burnman_planet, burnman_layers)
        else:
            log(f'{planet_class.class_type} planet type detected.', level='debug')
            planet = planet_class(planet_config)

    if planet is None:
        raise Exception

    log('World Load Successful!')
    log('Note that the first calculations with a new world will be slow as numba compiles functions. '
        'Subsequent calls should speed up.', level='info')

    if exit_planets:
        atexit.register(planet.kill_world)

    return planet


def clean_planet_config(planet_config: dict, make_copy: bool = True):
    """ Provides a clean copy of a planet's configuration, deleting items initialized by TidalPy

    Parameters
    ----------
    planet_config : dict
        Planet's configuration dictionary
    make_copy : bool = True
        Determines if a copy of the dictionary is made. Otherwise changes in this function will affect any other
        pointers to this dict object.
    Returns
    -------
    cleaned_dict : dict
        Cleaned dictionary with no TidalPy entries
    """

    if make_copy:
        cleaned_dict = copy.deepcopy(planet_config)
    else:
        cleaned_dict = planet_config

    if 'TidalPy_version' in planet_config:
        del cleaned_dict['TidalPy_version']
    for layer_name, layer_dict in planet_config['layers'].items():
        if 'radii' in layer_dict:
            del cleaned_dict['layers'][layer_name]['radii']
        for thermal_param in ['thermal_conductivity', 'thermal_diffusivity', 'thermal_expansion', 'stefan',
                              'shear_modulus']:
            if thermal_param in layer_dict:
                del cleaned_dict['layers'][layer_name][thermal_param]

    return cleaned_dict


def build_from_planet(old_planet, new_config: dict, new_name: str = None):
    """ Constructs a new planet based on a previously built one.

    Parameters
    ----------
    old_planet :
        Already initialized planet object (this should have characteristics close to the one you desire)
    new_config : dict
        Planet config for the changes to the new planet. This will override configs found in old_planet.
    new_name : str
        Optional name provided to the new planet.

    Returns
    -------
    new_planet :
        The newly configured and initialized planet object.
    """

    # Make a deep copy of the original planet's config
    old_config = old_planet.config

    # Delete items that are most likely going to change
    old_config_copy = clean_planet_config(old_config, make_copy=True)

    # Combine new and old dictionaries allowing the new dict to over write the old
    combo_dict = {**old_config_copy, **new_config}

    # Make any additional changes to the configs before planet build
    change_name = True
    if 'name' in new_config:
        # User provided a new name, do nothing.
        new_name = new_config['name']
        change_name = False
    if new_name is None:
        # Come up with a new name so that it is clear that this new planet was built from a different one
        new_name = combo_dict['name']
        new_name += '_variant'
    if change_name:
        combo_dict['name'] = new_name

    # Build the planet
    new_planet = build_planet(planet_name=new_name, planet_config=combo_dict, force_build=True)
    return new_planet


def scale_from_planet(old_planet, mass_scale: float = None, radius_scale: float = None, new_name: str = None):
    """ Creates a new planet that is a scaled up/down version of an old planet

    The old_planet should not be too different from your goal planet in terms of mass, radius, layer structure, and
    composition. For example, if you are trying to make a slightly larger/smaller Earth, don't use Io or Europa as your
    old_planet. Likewise, if you want to make a 2x mass Pluto, don't use the Earth as your base.

    It is recommended to call the new planet's .paint() method before using it in any science to ensure its structure
    makes sense.

    Parameters
    ----------
    old_planet :
        Planet to base the scale off of. This should be an already initialized TidalPy planet object
    mass_scale : FloatNone
        Providing this will treat mass fractions of the planet's layer's as constants. Radii will be estimated from density
    radius_scale : FloatNone
        Providing this will treat volume fractions of the planet's layer's as constants. Masses will be calculated from EOS
    new_name : str
        New name to call the planet

    Returns
    -------
    new_planet :
        The newly scaled planet.
    """

    if mass_scale is not None:
        # TODO: This may require some sort of iterative method w/ Burnman
        raise ImplementationException('Right now planets can only be scaled by their radius')

    if radius_scale is None:
        raise MissingArgumentError('radius_scale must be provided to scale_planet function (until mass_scale is implemented).')

    log(f'Scaling planet from {old_planet.name}, with radius scale of {radius_scale}.', level='debug')
    print(f'Scaling planet from {old_planet.name}, with radius scale of {radius_scale}.')
    scaled_config = clean_planet_config(old_planet.config, make_copy=True)

    # Change layer radius
    scaled_config['radius'] = radius_scale * old_planet.config['radius']
    prev_layer_radius = 0.
    for layer_name, layer_dict in old_planet.config['layers'].items():
        old_radius = layer_dict['radius_upper']
        scaled_config['layers'][layer_name]['radius_upper'] = radius_scale * old_radius
        scaled_config['layers'][layer_name]['radius_lower'] = prev_layer_radius

        # Use this layer's upper radius as the next layer's lower radius
        prev_layer_radius = scaled_config['layers'][layer_name]['radius_upper']

        # Update other items
        scaled_config['layers'][layer_name]['radius'] = scaled_config['layers'][layer_name]['radius_upper']
        scaled_config['layers'][layer_name]['thickness'] = scaled_config['layers'][layer_name]['radius_upper'] - \
            scaled_config['layers'][layer_name]['radius_lower']

        # Delete incorrect information
        if 'cooling' in layer_dict:
            if 'density' in layer_dict['cooling']:
                del scaled_config['layers'][layer_name]['cooling']['density']
            if 'gravity' in layer_dict['cooling']:
                del scaled_config['layers'][layer_name]['cooling']['gravity']

    # Make new planet based on the scaled dictionary
    new_planet = build_from_planet(old_planet, new_config=scaled_config, new_name=new_name)

    return new_planet