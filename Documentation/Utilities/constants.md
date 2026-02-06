# Constants and Runtime Parameters
TidalPy uses several mathematical, scientific, and computational constants. As well as a number of parameters that are
set by the configuration file the first time TidalPy is imported (or it is reinitialized). Unlike other functions and
data in the Utilities section, these properties are defined at the project level in the "constants_.hpp" header and its
Cython wrapper: "constants.pyx / .pxd". In C++, data is stored in a struct with static members so they can be accessed
efficiently by all C++ code in TidalPy but still be updated via Cython/Python from user-provided configuration files.

In addition to constants set at compile time and those stored in `TidalPy.config`, there are some static parameters
that are loaded from 3rd party packages. An example would be Newton's gravitational constant which is loaded during
TidalPy's (re)initialization from [SciPy](https://docs.scipy.org/doc/scipy/reference/constants.html).

To avoid having to update this document every time a parameter is changed or added, please review the values found in
the latest version of [constants_.hpp](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/constants_.hpp) as well
as the default configuration file builder
[defaultc.py](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py).
