"""Helper functions for loading TidalPy's configurations.

It is recommended that you only change these once you have some experience with the package
You can check their default values by examining the same file at https://github.com/jrenaud90/TidalPy/TidalPy/defaultc.py
"""

import os
import warnings
from typing import Union
from itertools import islice

import toml

import TidalPy
from TidalPy import version
from TidalPy.exceptions import ConfigurationException, InitializationError
from TidalPy.paths import get_config_dir, get_worlds_dir, unique_path
from TidalPy.defaultc import default_config_str


def dict_replace_value(d_old: dict, d_new: dict) -> dict:
    merged_dict = {}
    for k, v in d_old.items():
        if isinstance(v, dict):
            if k in d_new:
                v = dict_replace_value(v, d_new[k])
        else:
            if k in d_new:
                v = d_new[k]
        merged_dict[k] = v
    return merged_dict

def save_dict_to_toml(dict_to_save: dict,
              file_path: str,
              overwrite: bool = True):
    """Saves a python dictionary to a toml file at the specified file path.

    Parameters
    ----------
    dict_to_save : dict
        Python dictionary.
    file_path : str
        Filepath to save to.
    overwrite : bool, default = True
        If True, then the file will be overwritten if already present. by default True
    """

    if '.toml' not in file_path:
        raise AttributeError('Please provide a toml file path (include ".toml" extension).')
    
    if type(dict_to_save) is not dict:
        raise AttributeError(f'Can only save python dictionaries to toml files, not {type(dict_to_save)}.')

    toml_output = None
    if os.path.isfile(file_path):
        if overwrite:
            os.remove(file_path)
        else:
            # Append a number to the config name until one is found that is not already in use.
            file_path = unique_path(file_path, is_dir=False, make_dir=False)
    
    with open(file_path, 'w') as toml_file:
        toml_output = toml.dump(dict_to_save, toml_file)
    return toml_output

def check_config_version(
        config_path: str,
        allow_bugfix_difference: bool = True,
        warn_on_false: bool = True,
        raise_on_false: bool = False) -> bool:
    """ Checks a TidalPy configuration file to ensure that it is compatible with this version of TidalPy.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file to test.
    allow_bugfix_difference : bool, default = True
        If true, then a config with version A.B.C will still be allowed for TidalPy version A.B.D
    warn_on_false : bool, default = True
        If true, then a warning message will be shown if the config version check fails.
    raise_on_false : bool, default = False
        If true, then an error message will be raised if the config version check fails.
    
    Returns
    -------
    compatible : bool
        Flag for if this configuration file is compatible.
    """
    compatible = False
    with open(config_path, 'r') as config_file:
        config_version_found = False
        for line in islice(config_file, 0, 10):  # Assume the version number is in the first 10 lines
            if 'version:' in line.lower():
                config_version = line.split(': ')[1].split('\n')[0].strip()
                config_version_found = True
                break
            
        if not config_version_found:
            if raise_on_false:
                raise ConfigurationException(f'Can not find configuration version in {config_file}.')
            else:
                if warn_on_false:
                    message = f'Could not determine version for TidalPy configuration file, {config_path}. ' + \
                              f'It may not be compatible with this version of TidalPy.\n'
                    warnings.warn(message)
                return False
        
        if config_version == version:
            compatible = True
        elif allow_bugfix_difference:
            config_sub_vers = config_version.split('.')
            tpy_sub_vers = version.split('.')
            if (config_sub_vers[0] == tpy_sub_vers[0]) and (config_sub_vers[1] == tpy_sub_vers[1]):
                compatible = True

    if not compatible:
        message = f'TidalPy configuration file, {config_path}, was built for a different version of TidalPy ' + \
                  f'({config_version} vs. {version}). Unexpected behavior may arise.\n'
        if raise_on_false:
            raise ConfigurationException(message)
        elif warn_on_false:
            warnings.warn(message)

    return compatible

def get_default_config() -> dict:
    """ Loads TidalPy configurations that are found on the local disk.
    If no configuration file is found (likely when TidalPy is used for the first time) then default configurations
    will be saved to disk first.
    """

    config_dir = get_config_dir()
    config_path = os.path.join(config_dir, 'TidalPy_Configs.toml')
    # Check if TidalPy's config file is not present.
    if not os.path.isfile(config_path):
        # Create toml file with default configurations.
        with open(config_path, 'w') as config_file:
            config_file.write(f'#===========================================================#\n')
            config_file.write(f'#    TidalPy Default Configurations for Version: {version}\n')
            config_file.write(f'#===========================================================#\n\n')
            config_file.write(default_config_str)
    else:
        # Check if configuration file is for the correct version of TidalPy.
        check_config_version(config_path)
            
    # Load configurations (these may have been changed by the user) to dict
    config_dict = toml.load(config_path)

    # Update path
    TidalPy._config_path = config_path

    return config_dict

def set_config(new_config_path: Union[str, dict]) -> dict:
    """Sets TidalPy's configuration based on a provided configuration file path.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file the user wishes to use. 
        if set to "default" then the default config will be used.
    """
    
    new_config_name = 'unknown'
    if isinstance(new_config_path, dict):
        new_config = new_config_path
        new_config_name = 'User-provided dict'
    elif isinstance(new_config_path, str):
        new_config_name = f'{new_config_path}'
        if new_config_path.lower() == 'default':
            # Use default path.
            new_config = get_default_config()
        else:
            # Check if file exists
            if not os.path.isfile(new_config_path):
                raise InitializationError(f'Provided configuration path is not a file: {new_config_path}.')
        
            # Check if the provided configuration file is for the correct version of TidalPy.
            check_config_version(new_config_path, warn_on_false=True, raise_on_false=False)

            # Update path
            TidalPy._config_path = new_config_path
            
            # Load configurations (these may have been changed by the user) to dict
            new_config = toml.load(new_config_path)
    else:
        raise TypeError("Unexpected type found for TidalPy config replacement. Expected configuration file filepath (str) or config (dict).")

        TidalPy.config = get_default_config()

    # Set or override configurations with this new config file.
    if TidalPy.config is None:
        # No config has been loaded. Use this as the base config.
        TidalPy.config = new_config
    else:
        # A base config has already been loaded, override the base with any items from this new config.
        TidalPy.config = dict_replace_value(TidalPy.config, new_config)
        if TidalPy._tidalpy_init:
            from TidalPy.logger import get_logger
            log = get_logger('TidalPy')
            log.debug(f"TidalPy Configs overridden by {new_config_name}.")

def get_default_world_dir() -> str:
    """ Find the directory containing TidalPy's world configuration files.
    If no directory is found (likely when TidalPy is used for the first time) then default configurations
    will be saved to disk first.
    """

    worlds_dir = get_worlds_dir()

    install_worlds = True
    # Use a test world file to check that the default worlds are installed.
    # TODO: Update extension if/when converting world configs to toml.
    io_config = os.path.join(worlds_dir, 'io.toml')
    if os.path.isfile(io_config):
        # TODO: Have a check here to see if world config version matches tidalpy and rebuild if it doesn't?
        install_worlds = False
    
    if install_worlds:
        # Install worlds to world config.
        tpy_path = os.path.dirname(os.path.realpath(__file__))
        world_config_zip = os.path.join(tpy_path, 'WorldPack', 'WorldPack.zip')
        if not os.path.isfile(world_config_zip):
            raise InitializationError("Can not find TidalPy's WorldPack. " + \
                                      "There may have been an issue during TidalPy's installation.")
        import zipfile
        with zipfile.ZipFile(world_config_zip, 'r') as zip_ref:
            zip_ref.extractall(worlds_dir)
        
        # Re-perform Io test.
        if not os.path.isfile(io_config):
            raise InitializationError("Can not find Io configuration after WorldPack installation.")
    
    return worlds_dir

def set_world_dir(world_dir_path: str):
    """Sets TidalPy's worlds config file directory based on a provided directory path.
    
    Parameters
    ----------
    world_dir_path : str
        Path to the worlds directory the user wishes to use. 
        if set to "default" then the default directory will be used.
    """

    if world_dir_path.lower() == 'default':
        # Use default path.
        TidalPy.world_config_dir = get_default_world_dir()
    else:
        # TODO: Check if the provided directory has files compatible with the correct version of TidalPy.
        TidalPy.world_config_dir = world_dir_path
