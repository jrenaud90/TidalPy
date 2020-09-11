import os
import shutil

from ... import tidalpy_loc


def clear_cache(verbose: bool = False):
    """ Clears TidalPy Cache including cached numba functions
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
