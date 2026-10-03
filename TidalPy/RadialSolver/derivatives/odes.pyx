# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# The headers this module compiles read the shared TidalPy configuration, so its pointer is wired here as well.
from TidalPy.constants cimport get_shared_config_address, set_tidalpy_config_ptr
set_tidalpy_config_ptr(get_shared_config_address())


def find_num_shooting_solutions(int layer_type, int is_static, int is_incompressible):
    """
    Return the number of independent shooting solutions for a layer.

    Parameters
    ----------
    layer_type : int
        0 = solid, 1 = liquid.
    is_static : int
        1 = static, 0 = dynamic.
    is_incompressible : int
        1 = incompressible, 0 = compressible.

    Returns
    -------
    num_solutions : int
    """
    return c_find_num_shooting_solutions(layer_type, is_static, is_incompressible)

