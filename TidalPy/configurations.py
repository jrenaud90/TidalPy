import os
import toml
import warnings

import TidalPy
from TidalPy.paths import get_config_dir
from TidalPy.defaultc import default_config_str
# General configurations that are used by TidalPy
# Only change these once you have some experience with the package
# You can check their default values by examining the same file at https://github.com/jrenaud90/TidalPy

def get_config() -> dict:
    """
    Loads TidalPy configurations that are found on the local disk.
    If no configuration file is found (likely when TidalPy is used for the first time) then default configurations
    will be saved to disk first.
    """

    config_dir = get_config_dir()
    config_path = os.path.join(config_dir, 'TidalPy_Configs.toml')
    # Check if TidalPy's config file is not present.
    if not os.path.isfile(config_path):
        # Create toml file with default configurations.
        with open(config_path, 'w') as config_file:
            config_file.write(f'# TidalPy Configurations for version: {TidalPy.__version__}\n\n')
            config_file.write(default_config_str)
    else:
        # Check if configuration file is for the correct version.
        with open(config_path, 'r') as config_file:
            config_version = config_file.readline().split(': ')[1].split('\n')[0].strip()
            if config_version != TidalPy.__version__:
                warnings.warn('TidalPy configuration file was built for a different version of TidalPy '
                              f'({config_version} vs. {TidalPy.__version__}). Unexpected behavior may result.\n'
                              f'It is suggested that you delete the configuration file at {config_path} so it can be '
                              'reinitialized with the correct version of TidalPy.')
            
    # Load configurations (these may have been changed by the user) to dict
    config_dict = toml.load(config_path)

    return config_dict
