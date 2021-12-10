import os
import shutil

from ... import tidalpy_loc


def clear_cache(verbose: bool = True):
    """ Clears TidalPy's cached functions (python cache and cached numba functions).

    Parameters
    ----------
    verbose : bool = True
        Prints the name of pycache directories as they are cleared.

    Returns
    -------
    success: bool
    """

    if verbose:
        print('TidalPy Directory:', tidalpy_loc)
        print('Clearing TidalPy Cache...')

    for subdir, dirs, files in os.walk(tidalpy_loc):

        # Python and Numba caches save to the __pycache__ dir
        if '__pycache__' in dirs:
            cache_dir = os.path.join(subdir, '__pycache__')
            if verbose:
                print('Deleting:', cache_dir)
            shutil.rmtree(cache_dir)

    return True
