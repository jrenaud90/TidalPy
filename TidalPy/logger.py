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

def get_header_text() -> str:
    # Build header text
    now = datetime.now()
    now_str = now.strftime('%x at %X')
    header = (
        f'----------------------------------------------------------------------------------',
        f'TidalPy - Tidal Heating Calculator and Orbital Evolver',
        f'Version: {TidalPy.__version__}',
        f'Primary Development by Joe Renaud, ca. 2016--2023',
        f'Found a bug or have a suggestion? Open a new issue at github.com/jrenaud90/TidalPy',
        f'----------------------------------------------------------------------------------',
        f'Run made on {now_str}.',
        f'Using Python {sys.version} on {sys.platform}.\n##\n\n'
        )
    header = '\n'.join(header)

    return header

class DeltaTimeFormatter(logging.Formatter):
    def format(self, record):
        duration = datetime.datetime.utcfromtimestamp(record.relativeCreated / 1000)
        record.delta = duration.strftime("%H:%M:%S::%f")
        return super().format(record)
FORMATTER = DeltaTimeFormatter('%(asctime)s(+%(delta)s) - %(levelname)-9s: %(message)s', "%Y-%m-%d %H:%M:%S")

def get_console_handler(error_stream=False):
    if error_stream:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(LOGGING_LEVELS[TidalPy._config['console_error_level']])
    else:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(LOGGING_LEVELS[TidalPy._config['console_level']])
    console_handler.setFormatter(FORMATTER)
    return console_handler

def get_file_handler() -> logging.FileHandler:
    """ Get file handler for TidalPy's logger. """
    assert TidalPy._config is not None

    if not TidalPy._config['logging']['write_log']:
        # User does not want log written to disk.
        return None
    if not TidalPy._config['logging']['write_log_notebook'] and TidalPy._in_jupyter:
        # User does not want log written while using Jupyter notebook; which we are in.
        return None

    if TidalPy._config['logging']['use_cwd']:
        log_dir = os.path.join(TidalPy._output_dir, 'Logs')
    else:
        log_dir = get_log_dir()
    # Ensure directory exists
    Path(log_dir).mkdir(parents=True, exist_ok=True)

    # Find log filepath
    log_name = timestamped_str('TidalPy', date=True, time=True, second=True, millisecond=False, preappend=False)
    log_name += '.log'
    log_path = os.path.join(log_dir, log_name)

    # Add header test to log file
    with open(log_path, 'w') as log_file:
        log_file.write(get_header_text())

    # Create handler
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(FORMATTER)
    file_handler.setLevel(LOGGING_LEVELS[TidalPy._config['file_level']])
    return file_handler

def get_logger(logger_name: str) -> logging.Logger:
    global FILE_HANDLER
    global STREAM_HANDLER
    global STREAM_ERR_HANDLER

    # Ensure handlers are set
    if FILE_HANDLER is None:
        FILE_HANDLER = get_file_handler()
    if STREAM_HANDLER is None:
        STREAM_HANDLER = get_console_handler(error_stream=False)
    if STREAM_ERR_HANDLER is None:
        STREAM_ERR_HANDLER = get_console_handler(error_stream=True)
    
    # Get logger class
    logger = logging.getLogger(logger_name)

    # Clear any handlers that might be present
    logger.handlers = list()

    # Add handlers
    if FILE_HANDLER is not None:
        logger.addHandler(FILE_HANDLER)
    if STREAM_HANDLER is not None:
        logger.addHandler(STREAM_HANDLER)
    if STREAM_ERR_HANDLER is not None:
        logger.addHandler(STREAM_ERR_HANDLER)
    
    # Turn off propagation
    logger.propagate = False

    return logger
