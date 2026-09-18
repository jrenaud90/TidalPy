# Constants and Runtime Parameters
TidalPy uses mathematical, scientific, and computational constants, along with parameters set from the configuration file the first time TidalPy is imported or reinitialized. Unlike the rest of the Utilities section, these are defined at the project level in the "constants_.hpp" header and its Cython wrapper, "constants.pyx / .pxd". In C++ the data sits in a struct with static members, so every part of TidalPy's C++ code reads them efficiently while Cython and Python can still update them from user-provided configuration files.

Some parameters are loaded from third-party packages instead. Newton's gravitational constant, for example, is read from [SciPy](https://docs.scipy.org/doc/scipy/reference/constants.html) during TidalPy's (re)initialization.

For the current values, see the latest [constants_.hpp](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/constants_.hpp) and the default configuration file builder [defaultc.py](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py).
