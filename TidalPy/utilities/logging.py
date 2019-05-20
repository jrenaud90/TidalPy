import warnings
import time
import sys
import atexit

from ..io import unique_path
from .. import verbose_level as global_verbose_level
from .. import __version__

HEADER_TEXT = (
    f'----------------------------------------------------------------------------------'
    f'TidalPy - Tidal Heating Calculator and Orbital Evolver'
    f'Version: {__version__}'
    f'Made by Joe Renaud, ca 2015--2020'
    f'Found a bug or have a suggestion? Open a new issue at github.com/jrenaud90/TidalPy'
    f'----------------------------------------------------------------------------------'
)
HEADER_TEXT = '\n'.join(HEADER_TEXT)


class TidalLogger:


    def __init__(self, verbose_level=global_verbose_level, use_timestamps: bool = True):

        self.init_time = time.time()
        self.verbose_level = verbose_level
        self.use_timestamps = use_timestamps

        self._last_offload = time.time()
        self._buffer = ''

        # Create log file
        self.filepath = unique_path('TidalPy.log', is_dir=False, preappend_run_dir=True)
        with open(self.filepath, 'w') as logfile:
            logfile.write(HEADER_TEXT)

        # Make sure log is cleaned up whenever python is closed or the process is finished
        atexit.register(self.cleanup)

    def record(self, text: str, level: int = 3, warning: bool = False):
        """ Records text to buffer and prints to display if verbose_level exceeds message level"""

        # Save to buffer
        if self._buffer == '':
            self._buffer += '\n'
        if self.use_timestamps:
            buffer_text = ''
            for l_i, line in enumerate(text.split('\n')):
                if l_i == 0:
                    buffer_text = f'{self.timestamp()} - {line}'
                else:
                    # Put spaces that match the timestamp space
                    buffer_text += '\n' + ' '*6 + ' '*3 + line
        else:
            buffer_text = text
        self._buffer += buffer_text
        self.check_buffer()

        # Print to console if applicable
        if self.verbose_level > 0:
            level = max(level, 5)
            if warning:
                level = 1
            if level == 0:
                return None
            if self.verbose_level >= level:
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
        if self.filepath:
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

    def __call__(self, text: str, level: int = 3, warning: bool = False):

        self.record(text, level, warning)
