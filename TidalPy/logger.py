import os
import sys
import logging
from datetime import datetime
from pathlib import Path

import TidalPy
from TidalPy.paths import get_log_dir, timestamped_str

FILE_HANDLER = None
STREAM_HANDLER = None
STREAM_ERR_HANDLER = None
LOG_FILE_INIT = False

LOGGING_LEVELS = {
    # Critical: A serious error, indicating that the program itself may be unable to continue running.
    'CRITICAL': logging.CRITICAL,
    # Error: Due to a more serious problem, the software has not been able to perform some function.
    'ERROR'   : logging.ERROR,
    # Warning: An indication that something unexpected happened, or indicative of some problem in the near future (e.g.,
    #    `dis space low`). The software is still working as expected.
    'WARNING' : logging.WARNING,
    # Info: Confirmation that things are working as expected.
    'INFO'    : logging.INFO,
    # Debug: Detailed information, typically of interest only when diagnosing problems.
    'DEBUG'   : logging.DEBUG
    }

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

def get_header_text() -> str:
    # Build header text
    now = datetime.now()
    now_str = now.strftime('%x at %X')
    header = (
        f'----------------------------------------------------------------------------------',
        f'TidalPy - Tidal Heating Calculator and Orbital Evolver',
        f'Version: {TidalPy.__version__}',
        f'Primary Development by Joe Renaud, ca. 2016--2025',
        f'Found a bug or have a suggestion? Open a new issue at github.com/jrenaud90/TidalPy',
        f'----------------------------------------------------------------------------------',
        f'Run made on {now_str}.',
        f'Using Python {sys.version} on {sys.platform}.\n##\n'
        )
    header = '\n'.join(header)

    return header

class DeltaTimeFormatter(logging.Formatter):
    def format(self, record):
        duration = datetime.utcfromtimestamp(record.relativeCreated / 1000)
        record.delta = duration.strftime("%H:%M:%S::%f")
        return super().format(record)
FORMATTER = DeltaTimeFormatter('%(asctime)s(+%(delta)s) - %(levelname)-9s: %(message)s', "%Y-%m-%d %H:%M:%S")

def get_console_handler(error_stream=False):
    
    if (not TidalPy.config['logging']['print_log_notebook']) and is_notebook():
        return None
    
    if error_stream:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(LOGGING_LEVELS[TidalPy.config['logging']['console_error_level']])
    else:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(LOGGING_LEVELS[TidalPy.config['logging']['console_level']])
    
    console_handler.setFormatter(FORMATTER)
    return console_handler

def get_file_handler() -> logging.FileHandler:
    """ Get file handler for TidalPy's logger. """
    if TidalPy.config is None:
        return None

    if TidalPy.config['logging'] is None:
        return None

    if not TidalPy.config['logging']['write_log_to_disk']:
        # User does not want log written to disk.
        return None
    
    if (not TidalPy.config['logging']['write_log_notebook']) and is_notebook():
        # User does not want log written while using Jupyter notebook; which we are in.
        return None
    if TidalPy._test_mode:
        # TidalPy tests are being run, don't write to disk.
        return None

    if TidalPy.config['logging']['use_cwd']:
        log_dir = os.path.join(TidalPy._output_dir, 'Logs')
    else:
        log_dir = get_log_dir()
    # Ensure directory exists
    Path(log_dir).mkdir(parents=True, exist_ok=True)

    # Find log filepath
    log_name = timestamped_str('TidalPy', date=True, time=True, second=True, millisecond=False, preappend=False)
    log_name += '.log'
    log_path = os.path.join(log_dir, log_name)

    # Create handler
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(FORMATTER)
    file_handler.setLevel(LOGGING_LEVELS[TidalPy.config['logging']['file_level']])
    return file_handler

# Initialize root logger. Its handlers will be used by other modules logger. 
_root_logger = logging.getLogger('TidalPy')

# Initialize handlers
def initialize_handlers():
    global LOG_FILE_INIT
    global FILE_HANDLER
    global STREAM_HANDLER
    global STREAM_ERR_HANDLER

    # Clear any handlers that might be present
    _root_logger.handlers = list()

    # Set base logging level to the lowest one (it will be overridden by the tidalpy config via handlers)
    _root_logger.setLevel(1)

    # Log file handler
    FILE_HANDLER = get_file_handler()

    if FILE_HANDLER is not None:
        _root_logger.addHandler(FILE_HANDLER)
        # Check if log file has been initialized
        if not LOG_FILE_INIT:
            # Add header text to log file
            with open(FILE_HANDLER.baseFilename, 'w') as log_file:
                log_file.write(get_header_text())
            LOG_FILE_INIT = True

    # Console handler
    STREAM_HANDLER = get_console_handler(error_stream=False)

    if STREAM_HANDLER is not None:
        _root_logger.addHandler(STREAM_HANDLER)

    # Console Error handler
    STREAM_ERR_HANDLER = get_console_handler(error_stream=True)

    if STREAM_ERR_HANDLER is not None:
        _root_logger.addHandler(STREAM_ERR_HANDLER)

def get_logger(logger_name: str) -> logging.Logger:
    # Get logger class
    logger = logging.getLogger(logger_name)

    # Perform any adjustments to the logger
    # None are currently required

    return logger

# Intercept exceptions and log them
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    
    log = get_logger("TidalPy")

    log.error("Uncaught Exception", exc_info=(exc_type, exc_value, exc_traceback))

    # raise exc_type(exc_value)

sys.excepthook = handle_exception
