from .utilities.io.logging import TidalLogger


# Make logger that will be used throughout TidalPy
log = TidalLogger()

# Let everyone know that the logger and any other initializations are done.
_tidalpy_init = True
