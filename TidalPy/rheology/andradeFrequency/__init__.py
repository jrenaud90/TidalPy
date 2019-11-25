from ...utilities.model_utils import find_all_models, build_model_default_inputs
from ..defaults import rheology_param_defaults
from . import andrade_frequency_models

parameter_info_loc = ('complex_compliance', 'andrade_frequency_dependence')

known_models, known_model_const_args, known_model_live_args = find_all_models(andrade_frequency_models)

solid_viscosity_model_defaults = build_model_default_inputs(known_model_const_args, rheology_param_defaults,
                                                            inner_keys=parameter_info_loc)