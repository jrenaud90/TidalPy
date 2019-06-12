import operator

from TidalPy.exceptions import (ImproperAttributeHandling, ParameterMissingError, ImplementedBySubclassError,
                                ReinitNotAllowedError, TidalPyException)
from .. import debug_mode, auto_write, __version__, log
from ..io import inner_save_dir
import copy
import warnings
import os
import json5
from TidalPy.types import list_like
from TidalPy.utilities.dict_tools import nested_get
from ..exceptions import ParameterMissingError, UnknownModelError, MissingArgumentError

from inspect import getmembers, isfunction, isclass
from numba.targets.registry import CPUDispatcher

from typing import Union, List


class TidalPyClass:
    """ All functional classes used in TidalPy inherit from this class
    """

    def __init__(self):

        self.pyname = f'{self.__class__.__name__}'


class ConfigHolder(TidalPyClass):
    """ A helper class to take in and hold onto configurations that are both user provided and a default (hardcoded)
    """

    default_config = None

    def __init__(self, replacement_config: dict = None, automate: bool = False):

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
                raise ParameterMissingError(f'Parameter {param_name} missing from user provided and default configurations.')
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
        self.pyname = f'{self.__class__.__name__}_{self.model}'
        self.searcher = None
        self.func = None
        self.inputs = tuple()
        self.live_inputs = tuple()
        self.live_input_func = None

        if function_searcher is not None and automate:
            if debug_mode:
                assert isinstance(function_searcher, ModelSearcher)
            self.searcher = function_searcher
            self.searcher.defaults_require_key = False
            self.searcher.default_config = self.config
            self.func, self.inputs, self.live_input_func = self.searcher.find_model(self.model)

        # Switch between calculate and calculate debug. Generally, _calculate_debug is a much slower function that
        #    includes additional checks
        self._calc = self._calculate

        # TODO: it would be nice to carry the documentation like the below, but it is not allowed as written.
        # self.calculate.__doc__ = self._calculate.__doc__
        if debug_mode:
            if '_calculate_debug' in self.__dict__:
                self._calc = getattr(self, '_calculate_debug')

                # TODO: see above
                # self.calculate.__doc__ = f'-DEBUG VERSION- {self.calculate.__doc__}'
            else:
                # Subclass did not implement a debug mode calculator, use regular one
                pass

    def calculate(self, *args, **kwargs):

        # Some models have inputs that need to be updated at each call
        if self.live_input_func is not None:
            self.live_inputs = self.live_input_func(self)

        for input_ in self.live_inputs:
            if input_ is None:
                raise ParameterMissingError

        return self._calc(*args, **kwargs)

    def reinit(self):
        raise ReinitNotAllowedError

    def _calculate(self, *args, **kwargs):
        raise ImplementedBySubclassError



class LayerModel(ModelHolder):

    config_key = None

    def __init__(self, layer, function_searcher = None, model_name = None, automate: bool = False):

        self.layer = layer
        self.layer_type = layer.type
        if self.config_key is not None:
            try:
                config = layer.config[self.config_key]
            except KeyError:
                log(f"User provided no model information for layer {layer.name}'s {self.__class__.__name__}, "
                    f'using defaults instead.', level='debug')
                config = None
        else:
            config = None

        # Select the sub-dictionary of the default config (based on layer type)
        if self.default_config is not None:
            if self.layer_type in self.default_config:
                self.default_config = copy.deepcopy(self.default_config[self.layer_type])

        if config is None and self.default_config is None:
            raise ParameterMissingError(f'Config was not provided for {self.__class__} and no defaults are set.')

        super().__init__(model_name=model_name, user_config=config,
                         function_searcher=function_searcher, automate=automate)

        self.pyname = f'{self.__class__}_{self.layer_type}_{self.model}'

general_func_reject_list = [
    'njit',
    ]


class ModelSearcher(ConfigHolder):

    additional_reject_list = None

    def __init__(self, module, default_parameters: dict = None, defaults_require_key: bool = True):

        self.default_config = default_parameters
        self.defaults_require_key = defaults_require_key
        self.default_key = None
        if self.defaults_require_key:
            automate = False
        else:
            automate = True

        super().__init__(replacement_config=None, automate=automate)

        # Generate a dictionary of functions
        self.known_models, self.args_needed, self.live_args = self.find_known_models(module)

    def is_function(self, potential_func) -> bool:
        """ Checks if a function is a python or numba function """

        func_check = False
        if isfunction(potential_func):
            func_check = True
        if isinstance(potential_func, CPUDispatcher):
            func_check = True

        if func_check:
            if potential_func.__name__ in general_func_reject_list:
                func_check = False
            elif self.additional_reject_list is not None:
                if type(self.additional_reject_list) in list_like:
                    if potential_func.__name__ in self.additional_reject_list:
                        func_check = False
                else:
                    raise TypeError

        return func_check

    def find_known_models(self, module):
        # Generate a dictionary of functions
        func_list = getmembers(module, self.is_function)
        known_models = {name: func for name, func in func_list}

        # Parse the functions' doc string for function-call information
        args_needed = dict()
        live_args = dict()

        # See documentation file "Model Building.md" for more information on how custom TidalPy functions work.
        for model, func in known_models.items():

            args_needed[model] = None
            live_args[model] = None
            oargs_found = False
            largs_found = False
            for line in func.__doc__.split('\n'):
                if oargs_found and largs_found:
                    break

                # The doc string format for all functions may contain a line that is
                #     "live args: self.<name>.<name>...etc.arg1, self.<name>.<name>...etc.arg2, ..."
                if 'live args:' in line:
                    if oargs_found:
                        raise TidalPyException('Live args must come before other args in doc string (and in function signature)')
                    args = line.split('live args:')[-1].split(',')
                    cleaned_args = [arg.strip() for arg in args]
                    new_live_args = None
                    if len(cleaned_args) == 1:
                        if cleaned_args[0].lower() in ['none']:
                            cleaned_args = None
                    if cleaned_args is not None:
                        live_arg_list = list()
                        for live_arg in cleaned_args:
                            live_arg_sigs = live_arg.split('.')
                            if live_arg_sigs[0] != 'self':
                                raise TidalPyException('Live args must start with "self".')
                            live_ara_signature = '.'.join(live_arg_sigs[1:])
                            live_arg_list.append(operator.attrgetter(live_ara_signature))
                            new_live_args = lambda _self: tuple([_input(_self) for _input in live_arg_list])
                    live_args[model] = new_live_args
                    largs_found = True

                # The doc string format for all functions may contain a line that is "other_args: arg1, arg2, ..."
                if 'other args:' in line:
                    args = line.split('other args:')[-1].split(',')
                    cleaned_args = [arg.strip() for arg in args]
                    if len(cleaned_args) == 1:
                        if cleaned_args[0].lower() in ['none']:
                            cleaned_args = None

                    if type(cleaned_args) == list:
                        cleaned_args = tuple(cleaned_args)
                    args_needed[model] = cleaned_args
                    oargs_found = True

        return known_models, args_needed, live_args

    def find_model(self, model_name: str = None, parameters: dict = None, default_key: Union[str,List[str]] = None):

        # Update self.config based on function input
        if default_key is None:
            default_key = self.default_key
        if default_key is None:
            if self.defaults_require_key:
               raise MissingArgumentError
            else:
                defaults = self.default_config
        else:
            defaults = nested_get(self.default_config, default_key, raiseon_nolocate=True)

        if parameters is not None:
            user = parameters
        else:
            user = dict()
        config = {**defaults, **user}

        # Find Model
        try:
            if model_name is None:
                model_name = config['model']
        except KeyError:
            raise MissingArgumentError('No user provided Model and no fallback found in defaults')

        if model_name not in self.known_models:
            if model_name.lower() in self.known_models:
                model_name = model_name.lower()
            else:
                raise UnknownModelError(f'Unknown model provided to {self.pyname}.')
        model_func = self.known_models[model_name]
        needed_args = self.args_needed[model_name]
        live_args = self.live_args[model_name]

        # Build tuple of function inputs
        old_config = self.config
        self._config = config
        inputs = self.build_inputs(needed_args)
        self._config = old_config

        return model_func, inputs, live_args

    def build_inputs(self, needed_args) -> tuple:
        """ Builds an input tuple based on a function's needed arguments along with default
            and user-provided parameters"""

        inputs = list()
        if needed_args is not None:
            for arg_name in needed_args:
                # Arguments will default to the default parameter list (if present)
                # They will then be overridden by user input (if provided)
                # If neither of these sources have the required parameter, then an exception will be raised.
                arg = self.get_param(arg_name)
                if arg is None:
                    raise ParameterMissingError(f'required argument:"{arg_name}" not found in model defaults or in '
                                                f'user provided parameters.')
                inputs.append(arg)

        return tuple(inputs)

    def __call__(self, model_name: str, parameters: dict = None, default_key: Union[str, List[str]] = None):

        # Wrapper for self.find_model
        return self.find_model(model_name, parameters, default_key)