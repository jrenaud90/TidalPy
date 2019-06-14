from TidalPy.exceptions import UnknownModelError, ParameterMissingError, IncompatibleModelError, MissingArgumentError
from TidalPy.utilities.dict_tools import nested_get
from ..utilities.model import ModelSearcher
from . import compliance_models
from . import andrade_frequency_models
from inspect import getmembers
from typing import Union, List
from .defaults import rheology_param_defaults

class ComplianceModelSearcher(ModelSearcher):

    def __init__(self, compliance_module, frequency_module,
                 default_parameters: dict = None, defaults_require_key: bool = True):

        super().__init__(compliance_module, default_parameters=default_parameters,
                         defaults_require_key=defaults_require_key)

        # Also find the frequency models and their arguments
        self.known_frequency_models, self.frequency_args_needed, _ = self.find_known_models(frequency_module)

    def find_model(self, model_name: str = None, parameters: dict = None, default_key: Union[str, List[str]] = None):
        """ Searches known models for model_name and returns the function and required inputs """

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

        # Find Compliance Model
        try:
            if model_name is None:
                model_name = config['model']
        except KeyError:
            raise MissingArgumentError('No user provided Model and no fallback found in defaults')

        # Add in frequency check to the parameters
        frequency_model = config['andrade_frequency_model']
        use_frequency = False
        if model_name[-5:] == '_freq':
            use_frequency = True

        if use_frequency:
            if 'andrade' not in model_name or 'sundberg' not in model_name:
                raise IncompatibleModelError('Only the Andrade and Sundberg-Cooper rheologies can have '
                                             'additional frequency dependency.')
            frequency_func = self.known_frequency_models[frequency_model]
            needed_frequency_args = self.frequency_args_needed[frequency_model]

            old_config = self.config
            self._config = config
            frequency_inputs = self.build_inputs(needed_frequency_args)
            self._config = old_config

            # Add the frequency model information to the parameters dict. It will be used in the next setup step.
            self.config['andrade_freq_params'] = frequency_inputs
            self.config['andrade_freq_func'] = frequency_func

        return super().find_model(model_name, parameters, default_key)

    def __call__(self, model_name: str, parameters: dict = None, default_key: Union[str, List[str]] = None):

        # Wrapper for self.find_model
        return self.find_model(model_name, parameters, default_key)


find_compliance_func = ComplianceModelSearcher(compliance_models, andrade_frequency_models, rheology_param_defaults)