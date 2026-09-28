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

# TidalPy configurations (loaded from TidalPy_Configs.toml)
config = None

# Load the TidalPy initializer and run it (user can run it later so load it with the handle `reinit`)
from TidalPy.initialize import initialize as reinit

# Call reinit for the first initialization
reinit()

# Import module functions
from .cache import clear_cache as clear_cache
from .cache import clear_data as clear_data

# Save the effective configuration, headed by the package versions that produced it.
from .configurations import save_config as save_config


def test_mode():
    """ Turn on test mode and reinitialize TidalPy """
    global _test_mode

    if not _test_mode:
        _test_mode = True
        reinit()


def log_to_file():
    """ Quick switch to turn on saving logs to file """
    if not config['logging']['write_log_to_disk']:
        config['logging']['write_log_to_disk'] = True
        reinit()


def get_include() -> list:
    """ Directories holding TidalPy's C++ headers, and CyRK's, for packages that compile against them.

    Similar to ``numpy.get_include``. The headers include one another both by relative path and by bare file name,
    so every TidalPy directory that holds a header is listed, starting with the package root (``constants_.hpp``).
    The headers also need the header-only libraries TidalPy builds with (Eigen, xsf, and spdlog), which are not
    installed with TidalPy.

    Returns
    -------
    list of str
        CyRK's include directories followed by TidalPy's.
    """
    import CyRK

    include_dirs = list(CyRK.get_include())
    tidalpy_dir = os.path.dirname(os.path.abspath(__file__))
    for directory, sub_directories, file_names in os.walk(tidalpy_dir):
        sub_directories[:] = sorted(name for name in sub_directories if name != '__pycache__')
        if any(file_name.endswith('.hpp') for file_name in file_names):
            include_dirs.append(directory)
    return include_dirs
