import logging
import os
import sys
import warnings
from datetime import datetime

LOGGING_LEVELS = {
    # Critical: A serious error, indicating that the program itself may be unable to continue running.
    'CRITICAL': logging.CRITICAL,
    # Error: Due to a more serious problem, the software has not been able to perform some function.
    'ERROR': logging.ERROR,
    # Warning: An indication that something unexpected happened, or indicative of some problem in the near future (e.g.,
    #    `dis space low`). The software is still working as expected.
    'WARNING': logging.WARNING,
    # Info: Confirmation that things are working as expected.
    'INFO': logging.INFO,
    # Debug: Detailed information, typically of interest only when diagnosing problems.
    'DEBUG': logging.DEBUG
}


class LevelFilter(logging.Filter):
    def __init__(self, low, high):
        self._low = low
        self._high = high
        super().__init__()
    def filter(self, record):
        normal_filter = super().filter(record)
        if normal_filter:
            # Normal filter is okay, how about logging levels...
            if self._low <= record.levelno < self._high:
                return True
        return False


def log_setup(write_to_disk: bool = False, write_locale: str = None, running_in_jupyter: bool = False):
    """ Setup Python's logging module based on user provided information as well as built-in TidalPy settings

    Look at TidalPy.config or the /configurations.py file for switches that control the logging level and if the log
        automatically saves to disk.

    Parameters
    ----------
    write_to_disk : bool = False
        If True, the logger will save to the write_locale.
    write_locale : str = None
        Location that the logger will attempt to save to.
        If set to None, the logger will save to the current working directory.
    running_in_jupyter : bool = False
        If `True`, then the logger will not get a console handler (no logs printed to console)

    Returns
    -------
    log : logging.logger
        Global logger to be used throughout the TidalPy package.
    """
    from TidalPy import config
    from . import __version__

    # Build header text
    now = datetime.now()
    now_str = now.strftime('%x at %X')
    HEADER_TEXT = (
        f'----------------------------------------------------------------------------------',
        f'TidalPy - Tidal Heating Calculator and Orbital Evolver',
        f'Version: {__version__}',
        f'Primary Development by Joe Renaud, ca. 2016--2020',
        f'Found a bug or have a suggestion? Open a new issue at github.com/jrenaud90/TidalPy',
        f'----------------------------------------------------------------------------------',
        f'Run made on {now_str}.',
        f'Using Python {sys.version} on {sys.platform}.\n##\n\n'
    )
    HEADER_TEXT = '\n'.join(HEADER_TEXT)

    # Setup a global logger
    tidalpy_log = logging.getLogger('tidalpy')
    tidalpy_log.setLevel(LOGGING_LEVELS['DEBUG'])
    tidalpy_log.handlers = list()

    # Setup the log's format
    #    How the saved file looks...
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')
    #    How the console output looks...
    stream_formatter = logging.Formatter('TidalPy Log - %(levelname)s: %(message)s')

    # Setup logger filenames
    regular_log_filepath = None
    error_log_filepath = None
    if write_to_disk:
        outer_dir = write_locale
        if outer_dir is None:
            outer_dir = os.getcwd()
        regular_log_filepath = os.path.join(outer_dir, 'info_log.txt')
        error_log_filepath = os.path.join(outer_dir, 'error_log.txt')

        # Add TidalPy info text to the top of the regular log file.
        with open(regular_log_filepath, 'w') as regular_log_file:
            regular_log_file.write(HEADER_TEXT)

    # Setup handlers
    #    Console printer
    if not running_in_jupyter:
        # We do not want to print to console when we are running in a jupyter notebook
        reg_stream_handler = logging.StreamHandler(stream=sys.stdout)
        reg_stream_handler.setFormatter(stream_formatter)
        reg_stream_handler.addFilter(LevelFilter(LOGGING_LEVELS[config['stream_level']],
                                                 LOGGING_LEVELS[config['stream_err_level']]))
        tidalpy_log.addHandler(reg_stream_handler)

        err_stream_handler = logging.StreamHandler(stream=sys.stderr)
        err_stream_handler.setFormatter(stream_formatter)
        err_stream_handler.setLevel(config['stream_err_level'])
        tidalpy_log.addHandler(err_stream_handler)
    else:
        # Suppress warnings to console
        warnings.filterwarnings("ignore")
        tidalpy_log.propagate = False
        tidalpy_log.disabled = True

    #    File printer
    if write_to_disk:
        regular_file_handler = logging.FileHandler(regular_log_filepath)
        regular_file_handler.setFormatter(file_formatter)
        regular_file_handler.setLevel(config['regular_logfile_level'])
        error_file_handler = logging.FileHandler(error_log_filepath)
        error_file_handler.setFormatter(file_formatter)
        error_file_handler.setLevel(config['error_logfile_level'])

        # Add handlers to log
        tidalpy_log.addHandler(regular_file_handler)
        tidalpy_log.addHandler(error_file_handler)

    return tidalpy_log
