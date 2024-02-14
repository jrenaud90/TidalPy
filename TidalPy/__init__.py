# Find Version Number
import importlib.metadata
__version__ = importlib.metadata.version("TidalPy")
version = __version__

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
from .cache import clear_cache, clear_data
