"""Functions used to initialize and reinitialize TidalPy."""
import os
from pathlib import Path


def is_notebook() -> bool:
    """Whether TidalPy is running inside a Jupyter notebook (or qtconsole) kernel."""
    try:
        shell = get_ipython().__class__.__name__
        # ZMQInteractiveShell is a notebook or qtconsole; TerminalInteractiveShell is a terminal IPython.
        return shell == 'ZMQInteractiveShell'
    except NameError:
        # A standard Python interpreter.
        return False


def build_logging_config() -> dict:
    """Map the ``[logging]`` configuration section onto the C++ (spdlog) logger configuration.

    The console and file levels carry over. In a notebook the console keeps its messages for
    :func:`register_notebook_log_printer` to print below each cell, and prints only those at
    ``notebook_console_level`` and above (warnings, by default) unless ``print_log_notebook`` is set. The file sink
    is enabled only when ``write_log_to_disk`` is set (and ``write_log_notebook`` in a notebook) outside test mode.
    The log file is timestamped and lives in the run output directory (``use_cwd``) or the TidalPy data directory's
    ``Logs`` folder; without a usable data directory that second choice writes none.

    Returns
    -------
    dict
        Keys accepted by ``TidalPy.Utilities.logging.init_logger`` (empty when no configuration is loaded).
    """
    import TidalPy
    from TidalPy.paths import get_log_dir, timestamped_str
    from TidalPy.Utilities.logging.logger import resolve_log_level

    if not TidalPy.config or TidalPy.config.get('logging') is None:
        return {}
    logging_config = TidalPy.config['logging']
    in_notebook = is_notebook()

    console_level = logging_config['console_level']
    if in_notebook and not logging_config['print_log_notebook']:
        # The stricter of the two, so a notebook prints the warnings but never more than the console would.
        console_level = max(
            resolve_log_level(console_level), resolve_log_level(logging_config['notebook_console_level']))

    log_to_file = bool(logging_config['write_log_to_disk']) and not TidalPy._test_mode
    if in_notebook and not logging_config['write_log_notebook']:
        log_to_file = False

    log_file_path = ''
    if log_to_file:
        if logging_config['use_cwd']:
            log_dir = os.path.join(TidalPy._output_dir, 'Logs')
            Path(log_dir).mkdir(parents=True, exist_ok=True)
        else:
            # None when the data directory cannot be created (already warned about).
            log_dir = get_log_dir()
        if log_dir is None:
            log_to_file = False
        else:
            log_name = timestamped_str('TidalPy', date=True, time=True, second=True, millisecond=False,
                                       preappend=False) + '.log'
            log_file_path = os.path.join(log_dir, log_name)

    return {
        'console_level': console_level,
        'file_level': logging_config['file_level'],
        'log_to_file': log_to_file,
        'log_file_path': log_file_path,
        'console_pending': in_notebook,
    }


def register_notebook_log_printer():
    """Print the log messages a notebook's console keeps below each cell once the cell has run.

    Registers ``TidalPy.Utilities.logging.logger.print_pending_messages`` as an IPython ``post_run_cell`` hook, once
    per kernel. A notebook kernel does not show what C++ writes to the process's stdout on every platform, so the
    console sink keeps its messages (``console_pending``) and this hook writes them to the cell through
    ``sys.stderr``.
    """
    from TidalPy.Utilities.logging.logger import print_pending_messages

    events = get_ipython().events
    if print_pending_messages not in events.callbacks['post_run_cell']:
        events.register('post_run_cell', print_pending_messages)


def initialize(provided_config=None):
    """ Initialize (or reinitialize) TidalPy from its configuration, ``TidalPy.config``.

    Loads the configuration when none is loaded yet (the packaged defaults with the user's ``TidalPy_Configs.toml``
    merged over them), merges a ``TidalPy_Configs.toml`` found in the working directory when ``[configs]
    use_cwd_for_config`` is set, merges any override, sets up the run output directory and the logger, and pushes
    the numerical settings into the C++ config singleton. ``TidalPy.reinit`` is this function.

    Parameters
    ----------
    provided_config : str or dict, optional
        A configuration file path or dict merged over ``TidalPy.config``. ``"default"`` reloads the packaged defaults
        merged with the user's ``TidalPy_Configs.toml``. A file written by :func:`TidalPy.save_config` restores the
        settings of the run that saved it.
    """
    import TidalPy
    from TidalPy.configurations import get_default_config, save_config, set_config
    from TidalPy.constants import update_constants
    from TidalPy.paths import timestamped_str
    from TidalPy.Utilities.logging.logger import init_logger, log_debug

    TidalPy._in_jupyter = is_notebook()

    # Load the configuration if it is not already loaded.
    if TidalPy.config is None:
        get_default_config()

    # Merge a configuration found in the working directory, then one provided directly. A working directory without
    # one leaves the configuration as loaded.
    cwd_config_path = os.path.join(os.getcwd(), 'TidalPy_Configs.toml')
    if TidalPy.config['configs']['use_cwd_for_config'] and os.path.isfile(cwd_config_path):
        set_config(cwd_config_path)
    if provided_config is not None:
        set_config(provided_config)

    # Run output directory.
    output_dir = os.path.join(os.getcwd(), TidalPy.config['pathing']['save_directory'])
    if TidalPy.config['pathing']['append_datetime']:
        output_dir = timestamped_str(output_dir, date=True, time=True, second=False, millisecond=False,
                                     preappend=False)
    TidalPy._output_dir = output_dir

    # Logging.
    init_logger(build_logging_config())
    if TidalPy._in_jupyter:
        register_notebook_log_printer()
    if TidalPy._tidalpy_init:
        TidalPy._tidalpy_init = False
        log_debug('TidalPy reinitializing...')
    else:
        log_debug('TidalPy initializing...')

    # Save a copy of the loaded configuration to the run output directory.
    if TidalPy.config['configs']['save_configs_locally']:
        Path(TidalPy._output_dir).mkdir(parents=True, exist_ok=True)
        save_config(os.path.join(TidalPy._output_dir, 'TidalPy_Configs.toml'))
        log_debug(f'Output directory: {TidalPy._output_dir}')

    # Push the numerical settings into the C++ config singleton.
    update_constants()

    TidalPy._tidalpy_init = True
    log_debug('TidalPy initialization complete.')
