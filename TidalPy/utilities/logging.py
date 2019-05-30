import warnings
import time
import sys
import atexit

from ..io import unique_path
from .. import verbose_level as global_verbose_level
from .. import logging_level as global_logging_level
from .. import __version__
from typing import Union

HEADER_TEXT = (
    f'----------------------------------------------------------------------------------'
    f'TidalPy - Tidal Heating Calculator and Orbital Evolver'
    f'Version: {__version__}'
    f'Made by Joe Renaud, ca 2015--2020'
    f'Found a bug or have a suggestion? Open a new issue at github.com/jrenaud90/TidalPy'
    f'----------------------------------------------------------------------------------'
)
HEADER_TEXT = '\n'.join(HEADER_TEXT)

LEVEL_CONVERSION = {
    'debug': 5,
    'info': 4,
    'warning': 3,
    'error': 2,
    'critical': 1
}

class TidalLogger:


    def __init__(self, verbose_level=global_verbose_level, logging_level=global_logging_level,
                 use_timestamps: bool = True):

        self.init_time = time.time()
        self.verbose_level = verbose_level
        self.log_level = logging_level
        self.use_timestamps = use_timestamps

        self._last_offload = time.time()
        self._buffer = ''

        # Create log file
        self.filepath = unique_path('TidalPy.log', is_dir=False, preappend_run_dir=True)
        with open(self.filepath, 'w') as logfile:
            logfile.write(HEADER_TEXT)

        # Make sure log is cleaned up whenever python is closed or the process is finished
        atexit.register(self.cleanup)

    def record(self, text: str, level: Union[int, str] = 'info', warning: bool = False):
        """ Records text to buffer and prints to display if verbose_level exceeds message level"""

        # Find printing level
        if type(level) == str:
            level = LEVEL_CONVERSION[level]

        # Determine if text should be printed to console
        print_text = False
        if self.verbose_level > 0:
            level = max(level, 5)
            if warning and level > 3:
                level = 3
            if self.verbose_level >= level:
                print_text = True

        # Determine if text should be included in log file
        log_text = False
        if self.log_level > 0:
            level = max(level, 5)
            if warning and level > 3:
                level = 3
            if self.log_level >= level:
                log_text = True

        # Save to buffer
        if log_text:
            if self._buffer == '':
                self._buffer += '\n'
            if self.use_timestamps:
                buffer_text = ''
                for l_i, line in enumerate(text.split('\n')):
                    if l_i == 0:
                        buffer_text = f'{self.timestamp()},L{level} - {line}'
                    else:
                        # Put spaces that match the [timestamp=6] + [level info=3] + [spacer=3]
                        buffer_text += '\n' + ' ' * 6 + ' ' * 3 + ' ' * 3 + line
            else:
                buffer_text = text
            self._buffer += buffer_text
            self.check_buffer()

        if print_text:
            if warning:
                warnings.warn(text)
            else:
                print(text)

    def check_buffer(self):
        """ Checks if buffer should be off loaded to file """

        # If the buffer is too large then it must be offloaded or deleted to avoid memory overloads
        if self._buffer == '':
            return False

        if time.time() - self._last_offload > 30:
            # Only offload every 30 seconds or so.
            return self.offload()

    def offload(self, must_offload: bool = False):
        """ Offloads buffer contents to log file """

        if self._buffer == '':
            return True

        offloaded = False
        if self.filepath in [False, None]:
            if must_offload:
                warnings.warn('Can not offload buffer when no log filepath is provided. Buffer will be lost')
        else:
            with open(self.filepath, 'a') as logfile:
                logfile.write(self._buffer)
            offloaded = True

        # Clear Buffer
        self.clear()
        return offloaded

    def clear(self):
        """ Clears Text Buffer """

        self._buffer = ''

    def timestamp(self):
        """ Returns a time stamp in seconds since log started.
        :return: <str>
        """

        return f'{(time.time() - self.init_time):06}'

    def cleanup(self):
        """ Cleans up log when program is about to exit """

        self.record('TidalPy Closing...', level=3)
        self.offload(must_offload=True)

    def warn(self, text: str):

        self.record(text, warning=True)

    def __del__(self):

        self.cleanup()

    def __call__(self, text: str, level: Union[int, str] = 'info', warning: bool = False):

        self.record(text, level, warning)
