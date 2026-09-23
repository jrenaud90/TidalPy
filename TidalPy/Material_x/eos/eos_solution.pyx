# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# The headers this module compiles read the shared TidalPy configuration, so its pointer is wired here as well.
from TidalPy.constants cimport get_shared_config_address, set_tidalpy_config_ptr
set_tidalpy_config_ptr(get_shared_config_address())
