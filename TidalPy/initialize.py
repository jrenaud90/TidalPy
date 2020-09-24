""" Functions used to reinit TidalPy for a new run.
"""
import os
import sys

import json5
import numpy as np


def initialize_tidalpy():
    """ Initialize TidalPy based on information stored in TidalPy.config

    Items in TidalPy.config are identical to those in configurations.py unless the user changed them and called
        TidalPy.reinit()
    """
    import TidalPy

    from .logger import log_setup
    from .version import version
    from .io_helper import timestamped_str, unique_path

    # Are we in a Jupyter Notebook?
    running_in_jupyter = False
    for path in sys.path:
        if 'ipython' in path.lower():
            running_in_jupyter = True
            break
        elif 'jupyter' in path.lower():
            running_in_jupyter = True
            break
    TidalPy.running_in_jupyter = running_in_jupyter

    # Load configurations into the TidalPy.__init__
    TidalPy.debug_mode = TidalPy.config['debug_mode']

    # Other configurations
    format_numpy_floats = TidalPy.config['format_numpy_floats']
    if format_numpy_floats:
        float_formatter = lambda x: f'{x:0.3e}'
        np.set_printoptions(formatter={'float_kind': float_formatter})

    # Setup pathing
    tidalpy_loc = os.path.dirname(os.path.abspath(__file__))
    TidalPy.tidalpy_loc = tidalpy_loc

    save_dir = TidalPy.config['save_dir']
    if save_dir is None:
        save_dir = os.path.join(os.getcwd(), 'TidalPy_Output')
    save_to_disk = False
    if TidalPy.config['save_to_disk']:
        save_to_disk = True
        if running_in_jupyter and not TidalPy.config['save_to_disk_in_jupyter:']:
            save_to_disk = False

    # Create directories
    disk_loc = None
    if save_to_disk:
        # Make Outer Directory
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        inner_dir_str = timestamped_str(string_to_stamp='TidalPyRun',
                                        date=True, time=True, second=False, millisecond=False,
                                        preappend=False, separation='_')

        # Make inner `run` directory. Make sure it is unique.
        inner_dir_path = os.path.join(save_dir, inner_dir_str)
        inner_dir = unique_path(inner_dir_path, is_dir=True, make_dir=True)
        disk_loc = inner_dir
    TidalPy.disk_loc = disk_loc

    # Save a copy of TidalPy's current configurations to the save_dir
    if save_to_disk:
        config_file_name = f'TidalPyVers{version}_configurations.json'
        config_file_path = os.path.join(disk_loc, config_file_name)
        with open(config_file_path, 'w') as config_file:
            json5.dump(TidalPy.config, config_file, indent=4)

    # Initialize the loggers
    tidalpy_log = log_setup(save_to_disk, write_locale=disk_loc, running_in_jupyter=running_in_jupyter)
    TidalPy.log = tidalpy_log


    # Finish initialization
    TidalPy._tidalpy_init = True
