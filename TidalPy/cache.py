import os
import shutil

from TidalPy.paths import get_config_dir, get_log_dir, get_worlds_dir

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

    # Get install directory for TidalPy
    tidalpy_loc = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)

    if verbose:
        print('TidalPy Directory:', tidalpy_loc)
        print('Clearing TidalPy Cache...')

    for subdir, dirs, files in os.walk(tidalpy_loc):

        # Python and Numba caches save to the __pycache__ dir
        if '__pycache__' in dirs:
            cache_dir = os.path.join(subdir, '__pycache__')
            if verbose:
                print('Deleting: ', cache_dir)
            shutil.rmtree(cache_dir)

    return True

def clear_data(verbose: bool = True):
    """ Clears TidalPy's data files.
    
    Parameters
    ----------
    verbose : bool = True
        Prints the name of directories as they are cleared.

    Returns
    -------
    success: bool
    """

    dirs_to_del = list()
    dirs_to_del.append(get_config_dir())
    dirs_to_del.append(get_log_dir())
    dirs_to_del.append(get_worlds_dir())
    dir_str = '\n\t'.join(dirs_to_del)

    confirmation = input("Confirm that you would like to delete TidalPy's data directories? " + \
                         "Any edits made to TidalPy configurations and world configs will be lost. " + \
                         "Advise making backups first. " + \
                         f"The following directories will be removed: {dir_str}" + \
                         "\nProceed? (Y/N): ")
    
    if confirmation.lower() == 'y':
        # Delete files.
        for dir in dirs_to_del:
            if verbose:
                print('Deleting: ', dir)
            shutil.rmtree(dir)

        return True
    else:
        return False
