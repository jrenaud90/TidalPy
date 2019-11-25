from ...utilities.model_utils import find_all_models, build_model_defaults
from ..defaults import rheology_param_defaults
from . import melting_models as partial_melting_models
from .partialmelt import PartialMelt
from .partialmelt import calc_partial_melting

parameter_info_loc = ('partial_melting',)

known_models, known_model_const_args, known_model_live_args = find_all_models(partial_melting_models)

known_model_defaults = build_model_defaults(known_model_const_args, rheology_param_defaults,
                                            inner_keys=parameter_info_loc)
