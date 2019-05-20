from ..setup import version

import time
__version__ = version

# Initial Runtime
init_time = time.time()

# debug_mode is an optional runtime mode which will call upon many more checks. The checks will help minimize bugs,
#   but can slow down TidalPy. It is reccomended that you always run a case (at least partially) with debug_mode=True
#   first, and then once you have confirmed that no obvious bugs are present, you can re-run/finish the run with it off.
debug_mode = True

# auto_write determines if TidalPy will automatically create output directories and save logs, data, plots, etc. to it.
auto_write = True

# verbose_level determines what logging information will be displayed to standard output.
# Set to 0 or False for no printing including warnings
# Set to 5 for all printing
# verbose_level = 1 only permits warnings
# Does not affect print()
verbose_level = 3

# Make logger that will be used by other
from .utilities.logging import TidalLogger
log = TidalLogger()

