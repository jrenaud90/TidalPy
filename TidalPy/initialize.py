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

def initialize():
    """ Initialize (or reinitialize) TidalPy based on information stored in TidalPy.config

    Items in TidalPy.config are identical to those in the TidalPy_Config.toml unless the user changed them and called
        TidalPy.reinit()
    
    See more information about TidalPy_Config.toml in TidalPy.configurations.py
    """
    import TidalPy

    # Are we in a Jupyter Notebook?
    running_in_jupyter = is_notebook()
    TidalPy._in_jupyter = running_in_jupyter

    # Set TidalPy configurations if they are not already set.
    if TidalPy.config is None:
        # No configuration dictionary has been set.
        from TidalPy.configurations import set_config
        set_config('default')
    elif TidalPy.config['configs']['use_cwd_for_config']:
        set_config(os.path.join(os.getcwd(), 'TidalPy_Configs.toml'))

    # Setup pathing
    from TidalPy.paths import timestamped_str
    output_dir = os.path.join(os.getcwd(), TidalPy.config['pathing']['save_directory'])
    if TidalPy.config['pathing']['append_datetime']:
        output_dir = timestamped_str(output_dir, date=True, time=True, second=False, millisecond=False, preappend=False)
    TidalPy._output_dir = output_dir

    # Setup world configuration directory path
    from TidalPy.configurations import set_world_dir
    if TidalPy.config['configs']['use_cwd_for_world_dir']:
        set_world_dir(os.getcwd())
    else:
        set_world_dir('default')

    # Setup logging
    from TidalPy.logger import get_logger
    log = get_logger(__name__)
    # Reset initialization status
    if TidalPy._tidalpy_init:
        TidalPy._tidalpy_init = False
        log.info('TidalPy reinitializing...')
    else:
        log.info('TidalPy initializing...')

    if TidalPy.config['configs']['save_configs_locally'] or TidalPy.config['logging']['write_log_to_disk']:
        log.info(f'Output directory: {TidalPy._output_dir}')

    # Save a copy of TidalPy's current configurations to the save_dir
    if TidalPy.config['configs']['save_configs_locally']:
        # Create output directory if it does not exist
        Path(TidalPy._output_dir).mkdir(parents=True, exist_ok=True)

        # Save TidalPy configurations to that directory.
        config_file_name = f'TidalPy_Configs.toml'
        config_file_path = os.path.join(TidalPy._output_dir, config_file_name)
        with open(config_file_path, 'w') as config_file:
            toml.dump(TidalPy.config, config_file)

    # Load any other parameters from config to top-level program
    TidalPy.extensive_logging = TidalPy.config['debug']['extensive_logging']
    TidalPy.extensive_checks  = TidalPy.config['debug']['extensive_checks']

    # Finish initialization
    TidalPy._tidalpy_init = True
    log.info('TidalPy initialization complete.')
