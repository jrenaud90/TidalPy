import copy
import operator
from inspect import getmembers
from typing import List, Union, Tuple
import operator

from numba import njit
from .classes import ConfigHolder
from .dict_tools import nested_get
from .. import debug_mode
from ..exceptions import (ImplementedBySubclassError, MissingArgumentError, ParameterMissingError, TidalPyException,
                          UnknownModelError, AttributeNotSetError)
from ..initialize import log
from .functional_utils import is_function

general_func_reject_list = [
    njit,
]


class ModelHolder(ConfigHolder):

    """ Classes which are used in OOP calculation scheme
    """

    config_key = None
    known_models = None
    known_model_const_args = None
    known_model_live_args = None

    def __init__(self, model_name: str = None, replacement_config: dict = None, outer_key: str = None,
                 call_reinit: bool = True):

        super().__init__(replacement_config=replacement_config, call_reinit=call_reinit)

        if model_name is None and 'model' in self.config:
            model_name = self.config['model']

        # Store model name information
        self.model = model_name
        self.pyname = f'{self.__class__.__name__}_{self.model}'

        # Attempt to find function information
        self.func = None
        self._constant_arg_names = None
        self._live_arg_names = None
        try:
            self.func = self.known_models[self.model]
            self._constant_arg_names = self.known_model_const_args[self.model]
            self._live_arg_names = self.known_model_live_args[self.model]
        except KeyError:
            raise UnknownModelError(f'Unknown model: {self.model} for {self.__class__.__name__}')

        # Pull out constant arguments
        if self.config_key is not None:
            config = self.config[self.config_key]
        else:
            config = self.config
        self.inputs = self.build_args(self._constant_arg_names, parameter_dict=config, is_live_args=False)

        # Build live argument functions and calls
        self.live_inputs = None
        if len(self._live_arg_names) == 0:
            self.get_live_args = None
        else:
            live_funcs = self.build_args(self._live_arg_names, is_live_args=True)
            self.get_live_args = lambda: tuple([live_func(self) for live_func in live_funcs])

        # Switch between calculate and calculate debug. Generally, _calculate_debug is a much slower function that
        #    includes additional checks
        if debug_mode:
            self._calc = self._calculate_debug
        else:
            self._calc = self._calculate

    def calculate(self, *args, **kwargs):

        # Some models have inputs that need to be updated at each call
        if self.get_live_args is not None:
            try:
                self.live_inputs = self.get_live_args()
            except AttributeError:
                raise AttributeNotSetError('One or more live arguments are not set.')

        return self._calc(*args, **kwargs)

    def _calculate(self, *args, **kwargs):
        raise ImplementedBySubclassError

    def _calculate_debug(self, *args, **kwargs):
        raise ImplementedBySubclassError

    @staticmethod
    def build_args(arg_names: Tuple[str, ...], parameter_dict: dict = None, is_live_args: bool = False):
        """ Build an input tuple based on the required, constant, arguments and a parameter dictionary

        Parameters
        ----------
        arg_names : Tuple[str, ...]
            List of required constant argument names needed for this model's function.
        parameter_dict : dict
            Dictionary of parameters.
        is_live_args : bool
            Flag if this call is looking for live args, rather than constant args.

        Returns
        -------
        args : Tuple[Any, ...]
            List of default argument parameters
        """

        args = list()

        if not is_live_args:
            # Constant Arguments
            if parameter_dict is None:
                raise MissingArgumentError('Parameter configuration dictionary is required to build constant args.')

            for arg_name in arg_names:
                if arg_name not in parameter_dict:
                    raise ParameterMissingError(f'Parameter: {arg_name} is missing from configuration dictionary.')
                args.append(parameter_dict[arg_name])
        else:
            # Live Arguments
            for live_arg_signature in arg_names:

                # Create a attrgetter for the live argument.
                if 'self.' in live_arg_signature:
                    live_arg_signature = live_arg_signature.split('self.')[1]

                getter_func = operator.attrgetter(live_arg_signature)
                args.append(getter_func)

        return tuple(args)




class LayerModel(ModelHolder):
    config_key = None

    def __init__(self, layer, function_searcher=None, model_name=None, call_reinit: bool = True,
                 store_config_in_layer: bool = True):

        self.layer = layer
        self.layer_type = layer.type
        self.world = self.layer.world
        if self.world is not None:
            world_name = self.world.name
        else:
            world_name = 'Unknown_World'
        if self.config_key is not None:
            try:
                config = layer.config[self.config_key]
            except KeyError:
                log(f"User provided no model information for [layer: {layer.name} in world: {world_name}]'s "
                    f"{self.__class__.__name__}, using defaults instead.", level='debug')
                config = None
        else:
            config = None

        # Select the sub-dictionary of the default config (based on layer type)
        if self.default_config is not None:
            if self.layer_type in self.default_config:
                self.default_config = copy.deepcopy(self.default_config[self.layer_type])

        if config is None and self.default_config is None:
            raise ParameterMissingError(f"Config was not provided for [layer: {layer.name} in world: {world_name}]'s "
                                        f"{self.__class__} and no defaults are set.")

        super().__init__(model_name=model_name, user_config=config,
                         function_searcher=function_searcher, call_reinit=call_reinit)

        self.pyname = f'{self.__class__}_{self.layer_type}_{self.model}'


        if store_config_in_layer and self.config_key is not None:
            self.layer._config[self.config_key] = copy.deepcopy(self.config)


class ModelSearcher(ConfigHolder):
    additional_reject_list = None

    def __init__(self, module, default_parameters: dict = None, defaults_require_key: bool = True,
                 call_reinit: bool = True):

        self.default_config = default_parameters
        self.defaults_require_key = defaults_require_key
        self.default_key = None

        super().__init__(replacement_config=None, call_reinit=call_reinit)

        # Generate a dictionary of functions
        self.known_models, self.known_models_docs, self.args_needed, self.live_args = self.find_known_models(module)

    def is_function(self, potential_func) -> bool:
        """ Checks if a function is a python or numba function """

        func_check = is_function(potential_func)

        if func_check:
            if potential_func.__name__ in general_func_reject_list:
                func_check = False
            elif self.additional_reject_list is not None:
                if potential_func.__name__ in self.additional_reject_list:
                    func_check = False

        return func_check

    def find_known_models(self, module):
        # Generate a dictionary of functions
        func_list = getmembers(module, self.is_function)
        known_models = {name: func for name, func in func_list}
        known_models_docs = {name: func.__doc__ for name, func in func_list}

        # Parse the functions' doc string for function-call information
        args_needed = dict()
        live_args = dict()

        # See documentation file "Model Building.md" for more information on how custom TidalPy functions work.
        for model, func in known_models.items():

            args_needed[model] = None
            live_args[model] = None
            oargs_found = False
            largs_found = False
            if func.__doc__ is None:
                docstr = ''
            else:
                docstr = func.__doc__

            for line in docstr.split('\n'):
                if oargs_found and largs_found:
                    break

                # The doc string format for all functions may contain a line that is
                #     "live args: self.<name>.<name>...etc.arg1, self.<name>.<name>...etc.arg2, ..."
                if '!TPY_args live:' in line:
                    if oargs_found:
                        raise TidalPyException(
                                'Live args must come before other args in doc string (and in function signature)')
                    args = line.split('!TPY_args live:')[-1].split(',')
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
                            live_arg_signature = '.'.join(live_arg_sigs[1:])
                            live_arg_list.append(operator.attrgetter(live_arg_signature))
                            new_live_args = lambda _self: tuple([_input(_self) for _input in live_arg_list])
                    live_args[model] = new_live_args
                    largs_found = True

                # The doc string format for all functions may contain a line that is "other_args: arg1, arg2, ..."
                if '!TPY_args const:' in line:
                    args = line.split('!TPY_args const:')[-1].split(',')
                    cleaned_args = [arg.strip() for arg in args]
                    if len(cleaned_args) == 1:
                        if cleaned_args[0].lower() in ['none']:
                            cleaned_args = None

                    if type(cleaned_args) == list:
                        cleaned_args = tuple(cleaned_args)
                    args_needed[model] = cleaned_args
                    oargs_found = True

        return known_models, known_models_docs, args_needed, live_args

    def find_model(self, model_name: str = None, parameters: dict = None, default_key: Union[str, List[str]] = None):

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
