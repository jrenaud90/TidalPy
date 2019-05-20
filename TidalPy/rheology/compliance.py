from TidalPy.exceptions import UnknownModelError, ParameterMissingError, IncompatibleModelError
from ..utilities.search import ModelSearcher
from . import compliance_models
from . import andrade_frequency_models
from inspect import getmembers

class ComplianceModelSearcher(ModelSearcher):

    def __init__(self, module, frequency_module, default_parameters: dict = None):

        self.frequency_module = frequency_module
        # Generate a dictionary of functions
        func_list = getmembers(module, self.is_function)
        self.known_frequency_models = {name: func for name, func in func_list}

        self.frequency_args_needed = dict()
        for model, func in self.known_models:
            self.frequency_args_needed[model] = None
            for line in func.__doc__.split('\n'):
                # The doc string format for all functions should contain a line that is "other_args: arg1, arg2, ..."
                if 'other_params:' in line:
                    args = line.split('other_params:')[-1].split(',')
                    cleaned_args = [arg.strip() for arg in args]
                    if len(args) == 1:
                        if args[0].lower() in ['None']:
                            cleaned_args = None
                    self.frequency_args_needed[model] = cleaned_args
                    break

        super().__init__(module, default_parameters=default_parameters)

    def find_model(self, model_name: str, parameters: dict = None, default_key: str = None):
        """ Searches known models for model_name and returns the function and required inputs """

        if model_name not in self.known_models:
            raise UnknownModelError
        model_func = self.known_models[model_name]
        needed_args = self.args_needed[model_name]

        self.user_parameters = parameters
        if default_key is None:
            defaults = self.default_parameters
        else:
            defaults = self.default_parameters[default_key]

        # Add in frequency check to the parameters
        frequency_model = defaults.get('andrade_frequency_model', None)
        frequency_model = parameters.get('andrade_frequency_model', frequency_model)
        use_frequency = False
        if frequency_model is not None:
            if frequency_model.lower() not in ['off']:
                # If these two checks fail then the frequency model is turned off
                use_frequency = True

        if use_frequency and model_name[-5:] != '_freq':
            raise IncompatibleModelError
        if not use_frequency and model_name[-5:] == '_freq':
            raise IncompatibleModelError
        if model_name[-5:] == '_freq':
            if 'andrade' not in model_name or 'sundberg' not in model_name:
                raise IncompatibleModelError('Only the Andrade and Sundberg-Cooper rheologies are allowed to have'
                                             'additional frequency dependency.')
        if use_frequency:
            frequency_func = self.known_frequency_models[model_name]
            needed_frequency_args = self.frequency_args_needed[model_name]
            frequency_inputs = self.build_inputs(needed_frequency_args, defaults, user_provided=parameters)
            # Add the frequency model information to the parameters dict. It will be used in the next setup step.
            parameters['andrade_freq_params'] = frequency_inputs
            parameters['andrade_freq_func'] = frequency_func

        # Build tuple of function inputs
        inputs = self.build_inputs(needed_args, defaults, user_provided=parameters)

        return model_func, inputs

    def __call__(self, model_name: str, parameters: dict = None, default_key: str = None):

        # Wrapper for self.find_model
        return self.find_model(model_name, parameters, default_key)


find_compliance_func = ComplianceModelSearcher(compliance_models, andrade_frequency_models, compliance_param_defaults)