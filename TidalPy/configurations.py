import os
import toml

import TidalPy
from TidalPy.paths import get_config_dir
from TidalPy.defaultc import default_config_str
# General configurations that are used by TidalPy
# Only change these once you have some experience with the package
# You can check their default values by examining the same file at https://github.com/jrenaud90/TidalPy

# Check if configuration file can be found

# Important! Do not change `default_configurations`! These are overridden by user stored configs which are stored
# on your local disk. Location depends on OS:
#   - Windows: C:\Documents and Settings\<User>\Application Data\Local Settings\<AppAuthor>\<AppName>
#   - Linux: ~/.local/share/<AppName>
#   - MaxOS: ~/Library/Application Support/<AppName>
default_configurations = {

    # Log File
    # 'save_log_cwd'
    'save_log_cwd' : True, 

    # OLD CONFIGS -- REVIEW
    # Debug mode is an optional runtime mode which will call upon many more checks. The checks will help minimize bugs,
    #   but can slow down TidalPy. It is recommended that you always run a case (at least partially) with
    #   `debug_mode = True` first, and then once you have confirmed that no obvious bugs are present,
    #   you can re-run/finish the run with it off.
    'debug_mode'                           : True,

    # Save to Disk Preferences
    # # Can TidalPy save data to disk
    'save_to_disk'                         : False,
    # # Which directory should TidalPy use? If `None`, then TidalPy will make a new directory in the CWD
    'save_dir'                             : None,

    # Object/Planet save preferences
    'auto_save_object_config_to_rundir'    : True,
    'auto_save_object_config_to_tidalpydir': False,
    'overwrite_configs'                    : False,
    'give_configs_subscript'               : False,

    # numba.njit can speed up many functions, but it also makes debugging and error tracing more difficult.
    #  If you are having problems try setting this to False.
    'use_numba'                            : True,
    'cache_numba'                          : True,

    # Formatting preferences
    'format_numpy_floats'                  : True,

    # atexit invocations
    'exit_planets'                         : True,

    # Logger Config
    # # Determines if log will be printed to console or saved to drive when TidalPy is being used in a Jupyter Notebook
    'print_log_in_jupyter'                 : False,
    'write_log_in_jupyter'                 : True,
    # # Logging levels (set to 'DEBUG' during development)
    'stream_level'                         : 'DEBUG',
    'stream_err_level'                     : 'WARNING',
    'error_logfile_level'                  : 'WARNING',
    'regular_logfile_level'                : 'INFO',

    ## Physics Configurations
    'MIN_FREQ'                             : 1.0e-12,  # Minimum frequency (in rad s-1). Frequencies below this value will be treated as zero.
    }

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
            
    # Load configurations (these may have been changed by the user) to dict
    config_dict = toml.load(config_path)

    return config_dict
