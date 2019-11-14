from typing import List, Union

from . import andrade_frequency_models, compliance_models
from .defaults import rheology_param_defaults
from ..exceptions import IncompatibleModelConfigError, MissingArgumentError
from ..utilities.dict_tools import nested_get
from ..utilities.model import ModelSearcher


class ComplianceModelSearcher(ModelSearcher):
    additional_reject_list = (compliance_models.find_factorial.__name__,)

    def __init__(self, compliance_module, frequency_module,
                 default_parameters: dict = None, defaults_require_key: bool = True):

        super().__init__(compliance_module, default_parameters=default_parameters,
                         defaults_require_key=defaults_require_key)

        # Also find the frequency models and their arguments
        self.known_frequency_models, self.frequency_args_needed, _, _ = self.find_known_models(frequency_module)

    def find_model(self, model_name: str = None, parameters: dict = None, default_key: Union[str, List[str]] = None):
        """ Searches known models for model_name and returns the function and required inputs """

        if default_key is None:
            default_key = self.default_key
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
        self.config['andrade_freq_params'] = None
        self.config['andrade_freq_func'] = None
        use_frequency = False
        if model_name.lower() in ['andrade_freq', 'sundberg_freq', 'sundberg_cooper_freq']:

            use_frequency = True
            frequency_model = config['andrade_frequency_model']

        if use_frequency:
            if 'andrade' not in model_name or 'sundberg' not in model_name:
                raise IncompatibleModelConfigError('Only the Andrade and Sundberg-Cooper rheologies can have '
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

        model_func, inputs, live_args = super().find_model(model_name, parameters, default_key)

        if use_frequency:
            # Need to remove the old zeta from the inputs list
        inputs = list(inputs)

        return super().find_model(model_name, parameters, default_key)

    def __call__(self, model_name: str, parameters: dict = None, default_key: Union[str, List[str]] = None):

        # Wrapper for self.find_model
        return self.find_model(model_name, parameters, default_key)


find_complex_compliance = ComplianceModelSearcher(compliance_models, andrade_frequency_models, rheology_param_defaults)

# Build Tides Plotting Style
rheology_styles = dict()
for known_model, model_docs in find_complex_compliance.known_models_docs.items():
    # Initialize style with defaults
    style = {'color': 'black', 'ls': '-'}
    if model_docs is None:
        # If the model has no docstring then use the default
        pass
    else:
        for line in model_docs.split('\n'):
            if 'color' in line.lower():
                col = line.split('color:')[-1].strip().lower()
                style['color'] = col
            elif 'line style' in line.lower():
                ls = line.split('line style:')[-1].strip().lower()
                style['ls'] = ls
    rheology_styles[known_model] = style
