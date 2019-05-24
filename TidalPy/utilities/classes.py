from TidalPy.exceptions import ImproperAttributeHandling, ParameterMissingError, ImplementedBySubclassError
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

    def __init__(self, replacement_config: dict = None, default_config: dict = None, automate: bool = False):

        super().__init__()

        if debug_mode:
            assert type(replacement_config) == dict
            assert type(default_config) == dict

        # State Variables
        self._config_default = default_config
        self._config_replacement = replacement_config
        self._config = None

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
        pass

    def reinit(self):
        """ Performs any tasks that need to be done in order to reinitialize a class that might have been loaded from
            a dill/pickle file
        """

        self.update_config()
        self.init()

    def update_config(self) -> dict:
        """ Combines the default and user provided configurations into one dictionary.

        User parameters override defaults"""

        if self._config_replacement is None:
            self._config = copy.deepcopy(self.config_default)
        else:
            self._config = {**self.config_default, **self.config_user}

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
            config_to_save = self.config_default
            with open(config_filepath, 'w') as config_file:
                json5.dump(config_to_save, config_file)

        config_filepath = os.path.join(inner_save_dir, f'{self.pyname}.cfg')
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
    def config_replacement(self) -> dict:
        return self._config_replacement

    @config_replacement.setter
    def config_replacement(self, value: dict):

        if debug_mode:
            assert type(value) == dict

        self._config_replacement = copy.deepcopy(value)
        self.update_config()

    @property
    def config(self) -> dict:
        return self._config

    @config.setter
    def config(self, value):
        raise ImproperAttributeHandling('To change configurations set the "config_user" attribute '
                                        'or run "update_config"')


class ModelHolder(ConfigHolder):

    def __init__(self, model_name: str = None, user_config: dict = None, default_config: dict = None, function_searcher = None,
                 automate: bool = False):

        super().__init__(user_config, default_config, automate)
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
            self.func, self.inputs = self.searcher.find_model(self.model, parameters=self.config)

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

    def _calculate(self, *args, **kwargs):
        raise ImplementedBySubclassError



class LayerModel(ModelHolder):


    def __init__(self, layer: LayerType, default_config: dict = None,
                 function_searcher = None, model_name = None, automate: bool = False, config_key: str = None, ):

        self.layer = layer
        self.layer_type = layer.type
        if config_key is not None:
            config = layer.config[config_key]
        else:
            config = None
        if default_config is not None:
            default_config = default_config[self.layer_type]

        super().__init__(model_name=model_name, user_config=config, default_config=default_config,
                         function_searcher=function_searcher, automate=automate)
