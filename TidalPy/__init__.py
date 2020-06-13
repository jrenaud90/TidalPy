import os
import sys
import time

from .version import version


use_disk = True
_tidalpy_init = False
__version__ = version
_vers_major, _vers_minor, _vers_hotfix = version.split('.')
compatibility_signature = _vers_major + _vers_minor

# Location of TidalPy Install
tidalpy_loc = os.path.dirname(os.path.abspath(__file__))


def is_compatible(test_version: str):
    """ Tests rather or not a feature made in test_version is likely to be compatible in current __version__

    """
    test_vers_major, test_vers_minor, test_vers_hotfix = test_version.split('.')
    if int(test_vers_major) > int(_vers_major):
        return False
    if int(test_vers_major) < int(_vers_major):
        # Major releases may not be backwards compatible
        return False
    if int(test_vers_minor) > int(_vers_minor):
        return False
    if int(test_vers_minor) <= int(_vers_minor):
        # Minor releases should be backwards compatible
        pass

    return True


# Initial Runtime
init_time = time.time()

# debug_mode is an optional runtime mode which will call upon many more checks. The checks will help minimize bugs,
#   but can slow down TidalPy. It is reccomended that you always run a case (at least partially) with debug_mode=True
#   first, and then once you have confirmed that no obvious bugs are present, you can re-run/finish the run with it off.
debug_mode = True

# auto_write determines if TidalPy will automatically create output directories and save logs, data, plots, etc. to it.
auto_write = False

# verbose_level determines what logging information will be displayed to standard output.
# Set to 0 or False for no printing including warnings
# Set to 5 for all printing
# verbose_level = 1 only permits warnings
# Does not affect print()
if debug_mode:
    verbose_level = 5
    logging_level = 5
else:
    logging_level = 4
    verbose_level = 3

# Check to see if this is running in a jupyter notebook. If it is then it is usually not ideal to print log information
#   to the console.
from . import configurations
# Alias name for configurations file
config = configurations

if config.format_numpy_floats:
    # Format numpy array's so their floats don't print a ton of digits
    float_formatter = lambda x: f'{x:0.3e}'
    import numpy as np
    np.set_printoptions(formatter={'float_kind': float_formatter})


running_in_jupyter = False
for path in sys.path:
    if 'ipython' in path.lower():
        running_in_jupyter = True
    if 'jupyter' in path.lower():
        running_in_jupyter = True
    if running_in_jupyter:
        break
if running_in_jupyter and not config.print_log_in_jupyter:
    verbose_level = 0
if running_in_jupyter and not config.write_log_in_jupyter:
    logging_level = 0

# Some data files are stored in the TidalPy directory chain. Need to know where TidalPy is located to find those files.
tidalpy_dir = os.path.dirname(os.path.abspath(__file__))

# Data could be saved in other locations the user can tell TidalPy where to search for things like dilled planets, etc.
other_data_locs = [os.getcwd()]

# Bring up various functionality to top-level
from .utilities.performance import clear_cache

# Initialize logger
from .initialize import log