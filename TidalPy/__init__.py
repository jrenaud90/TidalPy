# Find Version Number
import importlib.metadata
__version__ = importlib.metadata.version("TidalPy")
version = __version__

# Set test_mode to False (used to turn off logging during testing)
_test_mode = False

import os
if 'TIDALPY_TEST_MODE' in os.environ:
    if os.environ['TIDALPY_TEST_MODE']:
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

# Helper function that provides directories to CyRK c++ headers
def get_include():
    import os
    # Since we depend on CyRK to build TidalPy; we likely want to include its headers as well.
    import CyRK
    tidalpy_dirs = CyRK.get_include()

    import TidalPy
    tidalpy_dir = os.path.dirname(TidalPy.__file__)

    # Utilities
    tidalpy_dirs += [
        # Utilities
        os.path.join(tidalpy_dir, 'utilities', 'lookups'),
        os.path.join(tidalpy_dir, 'utilities', 'arrays'),
        os.path.join(tidalpy_dir, 'utilities', 'dimensions'),

        # RadialSolver
        os.path.join(tidalpy_dir, 'RadialSolver'),

        # Material
        os.path.join(tidalpy_dir, 'Material', 'eos')
    ]

    return tidalpy_dirs
