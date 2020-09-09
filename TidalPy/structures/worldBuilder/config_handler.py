import copy

import json5

from ... import log
from ...utilities.io.pathing import get_all_files_of_type
from ..worldConfigs import WORLD_CONFIG_LOC

def clean_world_config(world_config: dict, make_copy: bool = True):
    """ Provides a clean copy of a world's configuration, deleting any items initialized by TidalPy

    Parameters
    ----------
    world_config : dict
        World's configuration dictionary
    make_copy : bool = True
        Determines if a copy of the dictionary is made. Otherwise changes in this function will affect any other
        pointers to this dict object.
    Returns
    -------
    cleaned_dict : dict
        Cleaned dictionary with no TidalPy entries
    """

    if make_copy:
        cleaned_dict = copy.deepcopy(world_config)
    else:
        cleaned_dict = world_config

    # Determine if this is a BurnMan world type
    is_burnman = world_config['type'] == 'burnman'

    if 'TidalPy_version' in world_config:
        del cleaned_dict['TidalPy_version']
    for layer_name, layer_dict in world_config['layers'].items():

        if 'radii' in layer_dict:
            del cleaned_dict['layers'][layer_name]['radii']

        if is_burnman:
            # Burnman world's will have some material properties set by the EOS. Delete any that are in the config.
            for thermal_param in ['thermal_conductivity', 'thermal_diffusivity', 'thermal_expansion', 'stefan',
                                  'shear_modulus']:
                if thermal_param in layer_dict:
                    del cleaned_dict['layers'][layer_name][thermal_param]

    return cleaned_dict

def check_for_duplicate_worlds(world_configs: dict):
    """ Check for duplicate worlds in the world config listing.



    Parameters
    ----------
    world_configs : dict
        Dictionary of worlds (in the format [<name>: <config filepath>])

    """

    potential_dups = list()
    for world_name, world_filepath in world_configs.items():
        possible_dup_index = world_name.split('_')[-1]
        try:
            int(possible_dup_index)
        except ValueError:
            continue
        else:
            if world_name not in potential_dups:
                potential_dups.append(world_name)

    for planet_name in potential_dups:
        log.warn(f'Possible duplicate saved planet found: {planet_name}. Ensure you use the correct subscript.')

# Find all planet configurations and import their config files
# Locate all planet configurations
known_planets_cfg = get_all_files_of_type(WORLD_CONFIG_LOC, ['cfg', 'json', 'json5'])
check_for_duplicate_worlds(known_planets_cfg)
_configs_are_dict = False

def _cfgpath_to_json():
    global known_planets_cfg
    global _configs_are_dict
    for planet_name, config_path in known_planets_cfg.items():
        with open(config_path, 'r') as config_file:
            known_planets_cfg[planet_name] = json5.load(config_file)

    _configs_are_dict = True