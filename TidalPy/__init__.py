import time
import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

# Initial Runtime
init_time = time.time()

# Version information
from .version import version
__version__ = version

# Load configuration dictionary
from .configurations import configurations
#    Alias the configurations dictionary
config = configurations

# Various properties to be set by the initializer
_tidalpy_init = False
use_disk = False
running_in_jupyter = False
tidalpy_loc = None  # type: str
disk_loc = None # type: str
log = None
#    debug_mode is an optional runtime mode which will call upon many more checks. The checks will help minimize bugs,
#       but can slow down TidalPy. It is reccomended that you always run a case (at least partially) with debug_mode=True
#       first, and then once you have confirmed that no obvious bugs are present, you can re-run/finish the run with it off.
debug_mode = False

# Load the TidalPy initializer and run it (user can run it later so load it with the handle `reinitialize`)
from .initialize import initialize_tidalpy as reinitialize
#    Alias the function
reinit = reinitialize
reinit()

# Try and find the world configurations data folder (if it has not been removed)
world_config_folder_found = False
world_config_loc = os.path.join(tidalpy_loc, 'WorldConfigs')
if os.path.isdir(world_config_loc):
    world_config_folder_found = True
else:
    # Make a new directory
    world_config_loc = os.path.join(tidalpy_loc, 'WorldConfigs-in_place')
    if not os.path.exists(world_config_loc):
        os.makedirs(world_config_loc)

# Bring up various functionality to top-level
#    Performance / Basic functionality
from .utilities.performance import clear_cache
#    Graphics
from .utilities.graphics import cmaps
#    Physics

#    OOP
from .structures import build_world, build_from_world, scale_from_world, Orbit