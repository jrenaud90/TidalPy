"""Functions used to initialize/reinitialize TidalPy"""
import os

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
    """ Initialize (or reinitialize) TidalPy based on information stored in TidalPy._config

    Items in TidalPy._config are identical to those in the TidalPy_Config.toml unless the user changed them and called
        TidalPy.reinit()
    
    See more information about TidalPy_Config.toml in TidalPy.configurations.py
    """
    import TidalPy

    # Are we in a Jupyter Notebook?
    running_in_jupyter = is_notebook()
    TidalPy._in_jupyter = running_in_jupyter

    # Get TidalPy configurations
    from TidalPy.configurations import get_config
    TidalPy._config = get_config()

    # Setup pathing
    from TidalPy.paths import timestamped_str
    output_dir = os.path.join(os.getcwd(), TidalPy._config['pathing']['outer_dir_name'])
    if TidalPy._config['append_datetime']:
        output_dir = timestamped_str(output_dir, date=True, time=True, second=False, millisecond=False, preappend=False)
    TidalPy._output_dir = output_dir

    # Setup logging
    from TidalPy.logger import get_logger
    log = get_logger(__name__)
    # Reset initialization status
    if TidalPy._tidalpy_init:
        TidalPy._tidalpy_init = False
        log.info('TidalPy reinitializing...')
    else:
        log.info('TidalPy initializing...')
    log.info(f'Output directory: {TidalPy._output_dir}')

    # Save a copy of TidalPy's current configurations to the save_dir
    if TidalPy._config['configs']['save_local_configs']:
        config_file_name = f'TidalPy_Configs.toml'
        config_file_path = os.path.join(TidalPy._output_dir, config_file_name)
        with open(config_file_path, 'w') as config_file:
            toml.dump(TidalPy._config, config_file)

    # Load any other parameters from config to top-level program
    TidalPy.extensive_logging = TidalPy._config['debug']['extensive_logging']

    # Finish initialization
    TidalPy._tidalpy_init = True
    log.info('TidalPy initialization complete.')
