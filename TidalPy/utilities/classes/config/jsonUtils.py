import copy
import os
from collections import OrderedDict
from typing import Tuple
from warnings import warn

import json5
import numpy as np

from ....io_helper import unique_path
from .... import use_disk, log

JSON5_KWARGS = {'indent': 4}

def clean_config_for_json(dict_to_clean: dict, keys_to_skip: Tuple[str, ...] = ('radii',)) -> OrderedDict:
    """ JSON does not like some data types. This is called right before a JSON dump, converting things as needed.

    Parameters
    ----------
    dict_to_clean : dict
        Dictionary to clean up.
        Ideally this will have strings as keys, but can have many different kinds of data types for values.
    keys_to_skip : Tuple[str, ...]
        Iterable of key strings which, if encountered, will be removed from the final cleaned dictionary.

    Returns
    -------
    json_config : OrderedDict
        Dictionary that should be suitable for JSON5 saving and loading.
        It is sorted alphabetically by key names.
    """

    json_config = dict()

    for key, item in dict_to_clean.items():

        dont_store = False
        if key in keys_to_skip:
            dont_store = True

        # Clean up the key
        input_key = None
        if type(key) != str:
            try:
                input_key = str(key)
            except TypeError:
                # Can't do much here. It won't be stored
                warn(f'Unusual key type encountered while trying to clean a dictionary for json. '
                     f'Key type = {type(key)}.')
                dont_store = True
        else:
            input_key = key

        # Clean up the value
        if dont_store:
            # Move on to next item
            continue

        if type(item) == dict or isinstance(item, OrderedDict):
            input_item = clean_config_for_json(item)
        elif type(item) == np.ndarray:
            if item.shape == tuple() or item.shape == (1,):
                input_item = float(item)
            else:
                input_item = item.tolist()

        elif type(item) in [np.float, np.float32, np.float64,
                            np.complex, np.complex64, np.complex128]:
            input_item = float(item)
        elif type(item) in [np.int8, np.int16, np.int32, np.int64,
                            np.uint8, np.uint16, np.uint32, np.uint64,
                            np.intp, np.uintp]:
            input_item = int(item)
        elif type(item) in [list, tuple, set, int, float, str, bool, complex, type(None)]:
            # Try to go forward with value
            input_item = item
        else:
            raise TypeError(f'Unrecognized type encountered while cleaning dictionary for json. Type = {type(item)}.')

        json_config[input_key] = input_item

    return OrderedDict(sorted(json_config.items()))


def save_dict_to_json(dict_to_save: dict, full_save_path: str, overwrite: bool = True):
    """ Save a python dictionary to disk

    Parameters
    ----------
    dict_to_save : dict
        Dictionary that should be saved to disk
    full_save_path : str
        Full, os-corrected including filename and extension, path to where file should be saved.
    overwrite : bool = True
        Flag if a file should be overwritten if already present.

    """

    if not use_disk:
        log.warning('Tried to write config to JSON file but use_disk set to False.')
    else:

        if not overwrite and os.path.isfile(full_save_path):
            warn('File found at config save path. Making a new unique save path.')
            full_save_path = unique_path(full_save_path, is_dir=False)

        # Make a copy of the dictionary to ensure that any cleaning or saving steps do not mutate the original dict.
        copy_of_dict = copy.deepcopy(dict_to_save)

        # Clean the dictionary for saving.
        clean_dict = clean_config_for_json(copy_of_dict)

        # Save to json
        with open(full_save_path, 'w') as config_file:
            json5.dump(clean_dict, config_file, **JSON5_KWARGS)
