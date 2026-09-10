"""Functions used to initialize/reinitialize TidalPy"""
import os
from pathlib import Path

import toml

def is_notebook() -> bool:
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

def build_logging_x_config() -> dict:
    """Map the classic ``[logging]`` settings onto the new backend's C++ (spdlog) logger configuration.

    Both loggers read the same section: the console and file levels carry over, the console is silenced in a
    notebook unless ``print_log_notebook`` is set, and the file sink is enabled only when ``write_log_to_disk``
    is set (and ``write_log_notebook`` in a notebook) outside test mode. The file is a separate, timestamped
    ``TidalPy_x`` log in the classic log directory, so the two loggers never share a file.

    Returns
    -------
    dict
        Keys accepted by ``TidalPy.Utilities_x.logging_x.init_logger`` (empty when no configuration is loaded).
    """
    import TidalPy
    from TidalPy.paths import get_log_dir, timestamped_str

    if not TidalPy.config or TidalPy.config.get('logging') is None:
        return {}
    logging_config = TidalPy.config['logging']
    in_notebook = is_notebook()

    console_level = logging_config['console_level']
    if in_notebook and not logging_config['print_log_notebook']:
        console_level = 'off'

    log_to_file = bool(logging_config['write_log_to_disk']) and not TidalPy._test_mode
    if in_notebook and not logging_config['write_log_notebook']:
        log_to_file = False

    log_file_path = ''
    if log_to_file:
        if logging_config['use_cwd']:
            log_dir = os.path.join(TidalPy._output_dir, 'Logs')
        else:
            log_dir = get_log_dir()
        Path(log_dir).mkdir(parents=True, exist_ok=True)
        log_name = timestamped_str('TidalPy_x', date=True, time=True, second=True, millisecond=False,
                                   preappend=False) + '.log'
        log_file_path = os.path.join(log_dir, log_name)

    return {
        'console_level': console_level,
        'file_level': logging_config['file_level'],
        'log_to_file': log_to_file,
        'log_file_path': log_file_path,
    }


def initialize(provided_config_file = None):
    """ Initialize (or reinitialize) TidalPy based on information stored in TidalPy.config

    Items in TidalPy.config are identical to those in the TidalPy_Config.toml unless the user changed them and called
        TidalPy.reinit()
    
    See more information about TidalPy_Config.toml in TidalPy.configurations.py
    """
    import TidalPy
    from TidalPy.constants import update_constants, update_constants_x

    # Are we in a Jupyter Notebook?
    running_in_jupyter = is_notebook()
    TidalPy._in_jupyter = running_in_jupyter

    # Set TidalPy configurations if they are not already set.
    from TidalPy.configurations import set_config
    if TidalPy.config is None:
        # No configuration dictionary has been set.
        set_config('default')

    # Load (or create) the configuration for the new `_x` class system. Stored on
    # TidalPy.config_x and written to TidalPy_Configs_x.toml on first use.
    from TidalPy.configurations import get_default_config_x
    if TidalPy.config_x is None:
        get_default_config_x()

    # Update default configs with any in the CWD
    if TidalPy.config['configs']['use_cwd_for_config']:
        set_config(os.path.join(os.getcwd(), 'TidalPy_Configs.toml'))
    
    # Update default configs with any provided directly to initialize.
    if provided_config_file is not None:
        set_config(provided_config_file)

    # Setup pathing
    from TidalPy.paths import timestamped_str
    if TidalPy.config:
        output_dir = os.path.join(os.getcwd(), TidalPy.config['pathing']['save_directory'])
        if TidalPy.config['pathing']['append_datetime']:
            output_dir = timestamped_str(output_dir, date=True, time=True, second=False, millisecond=False, preappend=False)
        TidalPy._output_dir = output_dir
    else:
        TidalPy._output_dir = os.getcwd()

    # Setup world configuration directory path
    from TidalPy.configurations import set_world_dir
    if TidalPy.config:
        if TidalPy.config['configs']['use_cwd_for_world_dir']:
            set_world_dir(os.getcwd())
        else:
            set_world_dir('default')
    else:
        set_world_dir('default')

    # Setup logging
    from TidalPy.logger import initialize_handlers, get_logger
    initialize_handlers()
    log = get_logger("TidalPy")

    from TidalPy.Utilities_x.logging_x.logger import init_logger
    init_logger(build_logging_x_config())
    # Reset initialization status
    if TidalPy._tidalpy_init:
        TidalPy._tidalpy_init = False
        log.debug('TidalPy reinitializing...')
    else:
        log.debug('TidalPy initializing...')

    if TidalPy.config:
        if TidalPy.config['configs']['save_configs_locally'] or TidalPy.config['logging']['write_log_to_disk']:
            log.debug(f'Output directory: {TidalPy._output_dir}')

    # Save a copy of TidalPy's current configurations to the save_dir
    if TidalPy.config:
        if TidalPy.config['configs']['save_configs_locally']:
            # Create output directory if it does not exist
            Path(TidalPy._output_dir).mkdir(parents=True, exist_ok=True)

            # Save TidalPy configurations to that directory.
            config_file_name = 'TidalPy_Configs.toml'
            config_file_path = os.path.join(TidalPy._output_dir, config_file_name)
            with open(config_file_path, 'w') as config_file:
                toml.dump(TidalPy.config, config_file)

    # Load any other parameters from config to top-level program
    if TidalPy.config:
        TidalPy.extensive_logging = TidalPy.config['debug']['extensive_logging']
        TidalPy.extensive_checks  = TidalPy.config['debug']['extensive_checks']
    else:
        TidalPy.extensive_logging = False
        TidalPy.extensive_checks  = False

    # Update constant values. update_constants_x() runs after the legacy update so
    # the shared C++ config singleton carries the `_x` config's numerical settings
    # for all `_x` modules.
    update_constants()
    if TidalPy.config_x:
        update_constants_x()

    # Finish initialization
    TidalPy._tidalpy_init = True
    log.debug('TidalPy initialization complete.')
