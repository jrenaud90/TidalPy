from TidalPy.types import list_like
from TidalPy.utilities.classes import ConfigHolder
from TidalPy.utilities.dict_tools import nested_get
from ..exceptions import ParameterMissingError, UnknownModelError, MissingArgumentError

from inspect import getmembers, isfunction, isclass
from numba.targets.registry import CPUDispatcher

from typing import Union, List

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
        self.known_models, self.args_needed = self.find_known_models(module)

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
        for model, func in self.known_models:
            args_needed[model] = None
            for line in func.__doc__.split('\n'):
                # The doc string format for all functions should contain a line that is "other_args: arg1, arg2, ..."
                if 'other args:' in line:
                    args = line.split('other args:')[-1].split(',')
                    cleaned_args = [arg.strip() for arg in args]
                    if len(args) == 1:
                        if args[0].lower() in ['None']:
                            cleaned_args = None
                    args_needed[model] = cleaned_args
                    break

        return known_models, args_needed

    def find_model(self, model_name: str = None, parameters: dict = None, default_key: Union[str,List[str]] = None):

        # Update self.config based on function input
        if default_key is None:
            default_key = self.default_key
        if default_key is None:
            if self.defaults_require_key:
                raise MissingArgumentError
        else:
            self.default_config = nested_get(self.default_config, default_key, raiseon_nolocate=True)
        if parameters is not None:
            self._user_config = parameters
        self.update_config()

        # Find Model
        try:
            if model_name is None:
                model_name = self.config['model']
        except KeyError:
            raise MissingArgumentError('No user provided Model and no fallback found in defaults')

        if model_name not in self.known_models:
            raise UnknownModelError
        model_func = self.known_models[model_name]
        needed_args = self.args_needed[model_name]

        # Build tuple of function inputs
        inputs = self.build_inputs(needed_args)

        return model_func, inputs

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