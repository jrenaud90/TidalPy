from TidalPy.exceptions import ImproperAttributeHandling, ParameterMissingError
from .. import debug_mode, auto_write, __version__
from ..io import run_dir
import copy
import warnings
import os
import json5

class TidalPyClass:
    """ All functional classes used in TidalPy inherit from this class
    """
    pass


class ConfigHolder(TidalPyClass):
    """ A helper class to take in and hold onto configurations that are both user provided and a default (hardcoded)
    """

    name = 'ConfigHolder'
    config_outer_key = None
    config_inner_key = None

    def __init__(self, user_config: dict = None):

        super().__init__()

        if debug_mode:
            assert type(user_config) == dict

        # State Variables
        self._config_default = None
        self._config_user = user_config
        self._config = None

        # Flags
        self.config_constructed = False

    def get_param(self, param_name: str, raise_missing: bool = True, fallback = None):
        """ Uses the outer and inner keys (if set) to parse into the configuration dictionary

        :param param_name:      <str> parameter or configuration name (capitalization and spaces matter!)
        :param raise_missing:   <bool> If set to true, an exception will be thrown if a parameter can not be found
        :param fallback:        Fall back if parameter can't be found and raise_missing is False
        :return:                Parameter or fallback
        """

        try:
            if self.config_outer_key is None:
                if self.config_inner_key is None:
                    try:
                        output = self.config[param_name]
                    except KeyError:
                        raise ParameterMissingError
                else:
                    try:
                        output = self.config[self.config_inner_key][param_name]
                    except KeyError:
                        raise ParameterMissingError
            else:
                if self.config_inner_key is None:
                    try:
                        output = self.config[self.config_outer_key][param_name]
                    except KeyError:
                        raise ParameterMissingError
                else:
                    try:
                        output = self.config[self.config_outer_key][self.config_inner_key][param_name]
                    except KeyError:
                        raise ParameterMissingError
        except ParameterMissingError as e:
            if raise_missing:
                raise e
            else:
                output = fallback

        return output

    def update_config(self) -> dict:
        """ Combines the default and user provided configurations into one dictionary.

        User parameters override defaults"""

        # Even though it will be copied all over the place, I think it is best to include as many refs to the version
        #   a config was made under. So we add it here to the defaults.
        self._config_default['TidalPy_version'] = __version__

        if self._config_user is None:
            self._config = copy.deepcopy(self.config_default)
        else:
            self._config = {**self.config_default, **self.config_user}

        self.config_constructed = True
        return self.config

    def save_config(self, save_default: bool = False) -> str:
        """ Saves final configuration to a JSON file """

        if not auto_write:
            warnings.warn('Tried to write config to JSON file but auto_write set to False.')

        if save_default:
            config_filepath = os.path.join(run_dir, f'{self.name}_default.cfg')
            config_to_save = self.config_default
        else:
            config_filepath = os.path.join(run_dir, f'{self.name}.cfg')
            config_to_save = self.config


        with open(config_filepath, 'w') as config_file:
            json5.dump(config_to_save, config_file)

        return config_filepath

    @property
    def config_default(self) -> dict:
        return self._config_default

    @config_default.setter
    def config_default(self, value: dict):

        if debug_mode:
            assert type(value) == dict

        if self._config_default is None:
            # We make a deep copy of the dictionary to ensure that no changes after class initialization propagate
            self._config_default = copy.deepcopy(value)
            self.update_config()
        else:
            raise ImproperAttributeHandling('Default config can only be set once!')

    @property
    def config_user(self) -> dict:
        return self._config_user

    @config_user.setter
    def config_user(self, value: dict):

        if debug_mode:
            assert type(value) == dict

        self._config_user = copy.deepcopy(value)
        self.update_config()

    @property
    def config(self) -> dict:
        return self._config

    @config.setter
    def config(self, value):
        raise ImproperAttributeHandling('To change configurations set the "config_user" attribute '
                                        'or run "update_config"')



