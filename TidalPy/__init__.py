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

