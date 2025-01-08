import copy

import toml

from TidalPy.paths import get_worlds_dir
from TidalPy.utilities.io.pathing import get_all_files_of_type

from TidalPy.logger import get_logger
log = get_logger("TidalPy")

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

    # Remove any old tidalpy version in case this was loaded from a pickled save.
    if 'TidalPy_version' in world_config:
        del cleaned_dict['TidalPy_version']

    # If there are layers, then clean the layer dicts
    if 'layers' in world_config:
        for layer_name, layer_dict in world_config['layers'].items():
            if 'radii' in layer_dict:
                del cleaned_dict['layers'][layer_name]['radii']

    return cleaned_dict


def check_for_duplicate_worlds(world_configs: dict):
    """ Check for duplicate world_types in the world config listing.

    Parameters
    ----------
    world_configs : dict
        Dictionary of world_types (in the format [<name>: <config filepath>])

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
        log.warning(f'Possible duplicate saved planet found: {planet_name}. Ensure you use the correct subscript.')


# Find all planet configurations and import their config files
# Locate all planet configurations
known_worlds_files = get_all_files_of_type(get_worlds_dir(), ['toml'])
known_worlds_cfg = dict()
check_for_duplicate_worlds(known_worlds_cfg)
_configs_loaded = False


def _cfgpath_to_toml():
    global known_worlds_cfg
    global _configs_loaded

    if not _configs_loaded:
        log.debug('Loading pre-built world configurations into the known world_types config dictionary.')
        for world_name, config_path in known_worlds_files.items():
            with open(config_path, 'r') as config_file:
                # Store with filename
                known_worlds_cfg[world_name] = toml.load(config_file)
            # Also store with name used in config
            if 'name' in known_worlds_cfg[world_name]:
                config_name = known_worlds_cfg[world_name]['name']
                if config_name != world_name:
                    known_worlds_cfg[config_name] = known_worlds_cfg[world_name]

    _configs_loaded = True


def get_world_configs():
    global known_worlds_cfg
    global _configs_loaded

    if not _configs_loaded:
        _cfgpath_to_toml()

    return known_worlds_cfg
