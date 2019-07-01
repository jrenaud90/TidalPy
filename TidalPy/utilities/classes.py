import copy
import os
import warnings
from collections import OrderedDict
from typing import List

import json5
import numpy as np

from .. import __version__, auto_write, debug_mode
from ..configurations import give_configs_subscript
from ..exceptions import ImproperAttributeHandling, ParameterMissingError
from ..io import inner_save_dir, unique_path


class TidalPyClass:
    """ All functional classes used in TidalPy inherit from this class
    """

    def __init__(self):
        self.pyname = f'{self.__class__.__name__}'


class ConfigHolder(TidalPyClass):
    """ A helper class to take in and hold onto configurations that are both user provided and a default (hardcoded)
    """

    default_config = None

    def __init__(self, replacement_config: dict = None, call_reinit: bool = True):

        super().__init__()

        if debug_mode:
            assert type(replacement_config) in [dict, type(None)]
            assert type(self.default_config) in [dict, type(None)]

        # State Variables
        self._user_config = replacement_config
        self._config = None
        self._old_config = None

        # Flags
        self.config_constructed = False

        # Merge the configurations
        self._old_config = copy.deepcopy(self.config)
        self.update_config()

        # Install user provided config
        if call_reinit:
            self.reinit()

    def reinit(self):
        """ Performs any tasks that need to be done in order to reinitialize a class that might have been loaded from
            a dill/pickle file
        """

        pass

    def replace_config(self, new_config: dict, delete_and_replace: bool = False):
        """ Replace the current user_config and update self.config with changes

        Parameters
        ----------
        new_config : dict
            New configuration dictionary
        delete_and_replace : bool = False
            If true then the old config will not be used as a baseline

        Returns
        -------
        config : dict
            self.config that has the old config overwritten by the new one (depending on if delete_and_replace is true).

        """

        self._old_config = copy.deepcopy(self.config)

        if delete_and_replace:
            self.default_config = None
        else:
            self.default_config = self.config
        self._user_config = new_config
        self.update_config()

        return self.config

    def get_param(self, param_name: str, raise_missing: bool = True, fallback=None):
        """ Uses the outer and inner keys (if set) to parse into the configuration dictionary

        :param param_name:      <str> parameter or configuration name (capitalization and spaces matter!)
        :param raise_missing:   <bool> If set to true, an exception will be thrown if a parameter can not be found
        :param fallback:        Fall back if parameter can't be found and raise_missing is False
        :return:                Parameter or fallback
        """

        try:
            try:
                output = self.config[param_name]
            except KeyError:
                raise ParameterMissingError
        except ParameterMissingError as e:
            if raise_missing:
                raise ParameterMissingError(
                    f'Parameter {param_name} missing from user provided and default configurations.')
            else:
                output = fallback

        return output

    def update_config(self) -> dict:
        """ Combines the default and user provided configurations into one dictionary.

        User parameters override defaults"""

        if self.default_config is None:
            self.default_config = dict()

        if self.user_config is None:
            self.user_config = dict()

        self._config = {**self.default_config, **self.user_config}

        # Even though it will be copied all over the place, I think it is best to include as many refs to the version
        #   a config was made under. So we add it here to the defaults.
        self._config['TidalPy_version'] = __version__

        self.config_constructed = True
        return self.config

    def save_config(self, save_default: bool = False, save_to_run_dir: bool = True,
                    additional_save_dirs: list = None, overwrite: bool = False) -> List[str]:
        """ Saves final configuration to a JSON file """

        json5_kwargs = {'indent': 4}

        if not auto_write:
            warnings.warn('Tried to write config to JSON file but auto_write set to False.')

        save_dirs = list()
        if save_to_run_dir:
            save_dirs.append(inner_save_dir)

        if additional_save_dirs is not None:
            save_dirs += additional_save_dirs

        config_filepaths = list()
        for save_dir in save_dirs:
            if save_default:
                config_filepath = os.path.join(save_dir, f'{self.pyname}_default.cfg')
                config_to_save = self.clean_config_for_json(self.default_config)
                if os.path.isfile(config_filepath) and not overwrite:
                    if give_configs_subscript:
                        config_filepath = unique_path(config_filepath, is_dir=False)
                    else:
                        config_filepath = None
                if config_filepath is not None:
                    with open(config_filepath, 'w') as config_file:
                        json5.dump(config_to_save, config_file, **json5_kwargs)

            config_filepath = os.path.join(save_dir, f'{self.pyname}.cfg')
            config_to_save = self.clean_config_for_json(self.config)
            if os.path.isfile(config_filepath) and not overwrite:
                if give_configs_subscript:
                    config_filepath = unique_path(config_filepath, is_dir=False)
                else:
                    config_filepath = None
            if config_filepath is not None:
                with open(config_filepath, 'w') as config_file:
                    json5.dump(config_to_save, config_file, **json5_kwargs)
            config_filepaths.append(config_filepath)

        return config_filepaths

    @property
    def user_config(self) -> dict:
        return self._user_config

    @user_config.setter
    def user_config(self, new_user_config: dict):

        if debug_mode:
            assert type(new_user_config) == dict

        self._user_config = copy.deepcopy(new_user_config)
        self.update_config()

    @property
    def config(self) -> dict:
        return self._config

    @config.setter
    def config(self, value):
        raise ImproperAttributeHandling('To change configurations set the "config_user" attribute '
                                        'or run "update_config"')

    @staticmethod
    def clean_config_for_json(config: dict) -> OrderedDict:
        """ JSON does not like some data types. This is called right before a JSON dump, converting things as needed.

        """

        json_config = dict()

        for key, item in config.items():

            dont_store = False
            if key in ['radii']:
                dont_store = True

            # Clean up the key
            input_key = None
            if type(key) != str:
                print(type(key))
                try:
                    input_key = str(key)
                except TypeError:
                    # Can't do much here. It won't be stored
                    dont_store = True
            else:
                input_key = key

            # Clean up the value
            if type(item) == dict or isinstance(item, OrderedDict):
                input_item = ConfigHolder.clean_config_for_json(item)
            elif type(item) == np.ndarray:
                if item.shape == tuple() or item.shape == (1,):
                    input_item = float(item)
                else:
                    input_item = [float(i) for i in item]
            elif type(item) in [list, tuple, set, int, float, str, bool, complex, type(None)]:
                # Try to go forward with value
                input_item = item
            else:
                print(type(item))
                raise TypeError

            if not dont_store:
                json_config[input_key] = input_item

        return OrderedDict(sorted(json_config.items()))
