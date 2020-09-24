# General configurations that are used by TidalPy
# Only change these once you have some experience with the package
# You can check their default values by examining the same file at https://github.com/jrenaud90/TidalPy

# Important Notes Regarding Changing Configurations
#     - You should always re-build planets when configurations have changed!
#           Previously built planets may not perform as anticipated under new configurations.

configurations = {
    # Debug mode is an optional runtime mode which will call upon many more checks. The checks will help minimize bugs,
    #   but can slow down TidalPy. It is recommended that you always run a case (at least partially) with
    #   `debug_mode = True` first, and then once you have confirmed that no obvious bugs are present,
    #   you can re-run/finish the run with it off.
    'debug_mode'                           : True,

    # Save to Disk Preferences
    'save_to_disk'                         : False,  # Can TidalPy save data to disk
    'save_dir'                             : None,
    # Which directory should TidalPy use? If `None`, then TidalPy will make a new directory in the CWD
    'save_to_disk_in_jupyter'              : False,
    # Overrides the above save_to_disk if TidalPy is being used in a Jupyter Notebook

    #    Object/Planet save preferences
    'auto_save_object_config_to_rundir'    : True,
    'auto_save_object_config_to_tidalpydir': False,
    'overwrite_configs'                    : False,
    'give_configs_subscript'               : False,

    # Configurations related to BurnMan Planet Building, Saving, and Loading
    'save_burnman_planets_to_tidalpy_dir'  : False,
    'save_burnman_planets_to_run_folder'   : True,
    'raise_on_changed_config'              : True,
    'force_burnman_quiet'                  : True,

    # Configurations related to how BurnMan results are used in TidalPy
    'burnman_interpolation_method'         : 'mid',  # Options: mid, avg, median
    'burnman_interpolation_N'              : 100,

    # numba.njit can speed up many functions, but it also makes debugging and error tracing more difficult.
    #  If you are having problems try setting this to False.
    'use_numba'                            : True,
    'cache_numba'                          : True,

    # Formatting preferences
    'format_numpy_floats'                  : True,

    # atexit invocations
    'exit_planets'                         : True,

    # Logger Config
    #    Determines if log will be printed to console or saved to drive when TidalPy is being used in a Jupyter Notebook
    'print_log_in_jupyter'                 : False,
    'write_log_in_jupyter'                 : True,
    #    Logging levels (set to 'DEBUG' during development)
    'stream_level'                         : 'DEBUG',
    'stream_err_level'                     : 'WARNING',
    'error_logfile_level'                  : 'WARNING',
    'regular_logfile_level'                : 'INFO'
}
