# Find Version Number
import importlib.metadata
__version__ = importlib.metadata.version("TidalPy")
version = __version__

# Set test_mode to False (used to turn off logging during testing)
_test_mode = False

import os
if os.environ.get('TIDALPY_TEST_MODE', '').strip().lower() in ('1', 'true', 'yes', 'on'):
    _test_mode = True

import time

# Initial Runtime
init_time = time.time()

# Various properties to be set by the configuration file and initializer (these should not be changed by user)
_tidalpy_init = False
_in_jupyter = False
_output_dir = None
_config_path = None

# TidalPy configurations
config = None

# World configuration directory
world_config_dir = None

# Public properties that can be changed by user
extensive_logging = False
extensive_checks = False

# Load the TidalPy initializer and run it (user can run it later so load it with the handle `reinitialize`)
from TidalPy.initialize import initialize as reinit

# Call reinit for the first initialization
reinit()

# Import module functions
from .cache import clear_cache as clear_cache
from .cache import clear_data as clear_data

# Announce the backend change coming in TidalPy 0.8.0, once per session.
import warnings as _warnings
from TidalPy.exceptions import TidalPyDeprecationWarning

_warnings.warn(
    "TidalPy 0.8.0 will replace this version's backend with a new C++ backend. The current modules (structures, tides," \
    "RadialSolver, rheology, ...) will be removed and module, class, and function names and signatures will change," \
    "so code written for TidalPy 0.7.X will need to be updated. Pin 'TidalPy<0.8' to keep using the current API." \
    "Silence this message with" \
    "warnings.filterwarnings('ignore', category=TidalPy.exceptions.TidalPyDeprecationWarning).",
    TidalPyDeprecationWarning,
    stacklevel=2)

def test_mode():
    """ Turn on test mode and reinitialize TidalPy """
    global _test_mode

    if _test_mode:
        # Don't need to do anything.
        pass
    else:
        _test_mode = True
        reinit()

def log_to_file():
    """ Quick switch to turn on saving logs to file """
    if not config['logging']['write_log_to_disk']:
        config['logging']['write_log_to_disk'] = True
        reinit()

# Helper function that provides directories to TidalPy's (and CyRK's) C++ headers
def get_include():
    """ Return the include directories of TidalPy's C++ source files, plus CyRK's, for dependent builds. """
    import CyRK

    # Since we depend on CyRK to build TidalPy; we likely want to include its headers as well.
    tidalpy_dirs = CyRK.get_include()

    tidalpy_dir = os.path.dirname(__file__)
    tidalpy_dirs += [
        # Utilities
        os.path.join(tidalpy_dir, 'utilities', 'arrays'),
        os.path.join(tidalpy_dir, 'utilities', 'dimensions'),

        # RadialSolver
        os.path.join(tidalpy_dir, 'RadialSolver'),

        # Material
        os.path.join(tidalpy_dir, 'Material', 'eos')
    ]

    return tidalpy_dirs
