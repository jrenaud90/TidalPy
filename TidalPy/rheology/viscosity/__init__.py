from ...utilities.model_utils import find_all_models, build_model_defaults
from ..defaults import rheology_param_defaults
from . import viscosity_models

solid_parameter_info_loc = ('solid_viscosity',)
liquid_parameter_info_loc = ('liquid_viscosity',)

known_models, known_model_const_args, known_model_live_args = find_all_models(viscosity_models)

solid_viscosity_model_defaults = build_model_defaults(known_model_const_args, rheology_param_defaults,
                                                      inner_keys=solid_parameter_info_loc)

liquid_viscosity_model_defaults = build_model_defaults(known_model_const_args, rheology_param_defaults,
                                                       inner_keys=liquid_parameter_info_loc)