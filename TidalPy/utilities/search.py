from ..exceptions import ParameterMissingError, UnknownModelError

from inspect import getmembers, isfunction, isclass
from numba.targets.registry import CPUDispatcher

general_func_reject_list = [
    'njit',
    ]


class ModelSearcher:

    reject_list = tuple()

    def __init__(self, module, default_parameters: dict = None):

        self.model_name = module

        # Generate a dictionary of functions
        func_list = getmembers(module, self.is_function)
        self.known_models = {name: func for name, func in func_list}

        # Parse the functions' doc string for calling information
        self.default_parameters = dict()
        if default_parameters is not None:
            self.default_parameters = default_parameters

        self.args_needed = dict()
        for model, func in self.known_models:
            self.args_needed[model] = None
            for line in func.__doc__.split('\n'):
                # The doc string format for all functions should contain a line that is "other_args: arg1, arg2, ..."
                if 'other_params:' in line:
                    args = line.split('other_params:')[-1].split(',')
                    cleaned_args = [arg.strip() for arg in args]
                    self.args_needed[model] = cleaned_args
                    break

        # Last user provided parameters
        self.user_parameters = None

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
            elif potential_func.__name__ in self.reject_list:
                func_check = False

        return func_check

    def find_model(self, model_name: str, parameters: dict = None, default_key: str = None):
        """ Searched known models for model_name and returns the function and required inputs """

        if model_name not in self.known_models:
            raise UnknownModelError
        model_func = self.known_models[model_name]
        needed_args = self.args_needed[model_name]
        inputs = None

        self.user_parameters = parameters
        if default_key is None:
            defaults = self.default_parameters
        else:
            defaults = self.default_parameters[default_key]

        # Build tuple of function inputs
        if needed_args is not None:
            inputs = list()
            for arg_name in needed_args:
                # Arguments will default to the default parameter list (if present)
                # They will then be overridden by user input (if provided)
                # If neither of these sources have the required parameter, then an exception will be raised.
                arg = defaults.get(arg_name, None)
                if parameters is not None:
                    arg = parameters.get(arg_name, arg)
                if arg is None:
                    raise ParameterMissingError(f'required argument:"{arg_name}" not found in model defaults or in '
                                                f'user provided parameters.')
                inputs.append(arg)
            inputs = tuple(inputs)

        return model_func, inputs

    def __call__(self, model_name: str, parameters: dict = None, default_key: str = None):

        # Wrapper for self.find_model
        return self.find_model(model_name, parameters, default_key)



