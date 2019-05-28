from TidalPy.exceptions import ImproperAttributeHandling, ParameterMissingError, ImplementedBySubclassError, ReinitNotAllowedError
from TidalPy.structures.layers import LayerType
from .. import debug_mode, auto_write, __version__
from ..io import inner_save_dir
from .search import ModelSearcher
import copy
import warnings
import os
import json5

class TidalPyClass:
    """ All functional classes used in TidalPy inherit from this class
    """

    def __init__(self):

        self.pyname = f'{self.__class__}'


class ConfigHolder(TidalPyClass):
    """ A helper class to take in and hold onto configurations that are both user provided and a default (hardcoded)
    """

    default_config = None

    def __init__(self, replacement_config: dict = None, automate: bool = False):

        super().__init__()

        if debug_mode:
            assert type(replacement_config) in [dict, None]
            assert type(self.default_config) in [dict, None]

        # State Variables
        self._user_config = replacement_config
        self._config = None
        self._old_config = None

        # Flags
        self.config_constructed = False

        # Install user provided config
        if automate:
            self.update_config()

    def get_param(self, param_name: str, raise_missing: bool = True, fallback = None):
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
                raise e
            else:
                output = fallback

        return output

    def init(self):
        """ Tasks that are could usually be done in __init__, but might need to be redone by self.reinit"""
        self._old_config = copy.deepcopy(self.config)

    def reinit(self):
        """ Performs any tasks that need to be done in order to reinitialize a class that might have been loaded from
            a dill/pickle file
        """

        self.update_config()
        self.init()

    def update_config(self) -> dict:
        """ Combines the default and user provided configurations into one dictionary.

        User parameters override defaults"""

        if self.user_config is None:
            self._config = copy.deepcopy(self.default_config)
        else:
            self._config = {**self.default_config, **self.user_config}

        # Even though it will be copied all over the place, I think it is best to include as many refs to the version
        #   a config was made under. So we add it here to the defaults.
        self._config['TidalPy_version'] = __version__

        self.config_constructed = True
        return self.config

    def save_config(self, save_default: bool = False) -> str:
        """ Saves final configuration to a JSON file """

        if not auto_write:
            warnings.warn('Tried to write config to JSON file but auto_write set to False.')

        if save_default:
            config_filepath = os.path.join(inner_save_dir, f'{self.pyname}_default.cfg')
            config_to_save = self.default_config
            with open(config_filepath, 'w') as config_file:
                json5.dump(config_to_save, config_file)

        config_filepath = os.path.join(inner_save_dir, f'{self.pyname}.cfg')
        config_to_save = self.config
        with open(config_filepath, 'w') as config_file:
            json5.dump(config_to_save, config_file)

        return config_filepath

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


class ModelHolder(ConfigHolder):

    def __init__(self, model_name: str = None, user_config: dict = None, function_searcher = None, automate: bool = False):

        super().__init__(replacement_config=user_config, automate=automate)
        if model_name is None:
            model_name = self.config['model']

        self.model = model_name
        self.pyname = f'{self.__class__}_{self.model}'
        self.searcher = None
        self.func = None
        self.inputs = None

        if function_searcher is not None and automate:
            if debug_mode:
                assert isinstance(function_searcher, ModelSearcher)
            self.searcher = function_searcher
            self.searcher.defaults_require_key = False
            self.searcher.default_config = self.config
            self.func, self.inputs = self.searcher.find_model(self.model)

        # Switch between calculate and calculate debug. Generally, _calculate_debug is a much slower function that
        #    includes additional checks
        if debug_mode:
            if '_calculate_debug' in self.__dict__:
                self.calculate = getattr(self, '_calculate_debug')
                if self.calculate.__doc__ is None and self._calculate.__doc__ is not None:
                    self.calculate.__doc__ = """DEBUG VERSION OF: """ + self._calculate.__doc__
            else:
                # Subclass did not implement a debug mode calculator
                self.calculate = self._calculate
        else:
            self.calculate = self._calculate

    def reinit(self):
        raise ReinitNotAllowedError

    def _calculate(self, *args, **kwargs):
        raise ImplementedBySubclassError



class LayerModel(ModelHolder):

    config_key = None

    def __init__(self, layer: LayerType, function_searcher = None, model_name = None, automate: bool = False):

        self.layer = layer
        self.layer_type = layer.type
        if self.config_key is not None:
            config = layer.config[self.config_key]
        else:
            config = None

        # Select the sub-dictionary of the default config (based on layer type)
        if self.default_config is not None:
            if self.layer_type in self.default_config:
                self.default_config = copy.deepcopy(self.default_config[self.layer_type])


        super().__init__(model_name=model_name, user_config=config,
                         function_searcher=function_searcher, automate=automate)

        self.pyname = f'{self.__class__}_{self.layer_type}_{self.model}'