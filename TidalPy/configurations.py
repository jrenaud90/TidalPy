"""Helper functions for loading TidalPy's configurations.

It is recommended that you only change these once you have some experience with the package
You can check their default values by examining the same file at https://github.com/jrenaud90/TidalPy/TidalPy/defaultc.py
"""

import os
import warnings
import toml

import TidalPy
from TidalPy import version
from TidalPy.exceptions import ConfigurationException, InitializationError
from TidalPy.paths import get_config_dir, get_worlds_dir
from TidalPy.defaultc import default_config_str


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
        config_version = config_file.readline().split(': ')[1].split('\n')[0].strip()
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
            config_file.write(f'# TidalPy Configurations for version: {version}\n\n')
            config_file.write(default_config_str)
    else:
        # Check if configuration file is for the correct version of TidalPy.
        check_config_version(config_path)
            
    # Load configurations (these may have been changed by the user) to dict
    config_dict = toml.load(config_path)

    return config_dict

def set_config(config_path: str) -> dict:
    """Sets TidalPy's configuration based on a provided configuration file path.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file the user wishes to use. 
        if set to "default" then the default config will be used.
    """
    
    if config_path.lower() == 'default':
        # Use default path.
        TidalPy.config = get_default_config()
    else:
        # Check if file exists
        if not os.path.isfile(config_path):
            raise InitializationError(f'Provided configuration path is not a file: {config_path}.')

        # Check if the provided configuration file is for the correct version of TidalPy.
        check_config_version(config_path)
    
        # Load configurations (these may have been changed by the user) to dict
        TidalPy.config = toml.load(config_path)

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
