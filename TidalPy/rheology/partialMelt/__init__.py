from ...utilities.model_utils import find_all_models, build_model_default_inputs
from .defaults import partial_melt_defaults
from . import melting_models as partial_melting_models

parameter_info_loc = ('rheology', 'partial_melting')

known_models, known_model_const_args, known_model_live_args = find_all_models(partial_melting_models)

def get_partial_melt_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args, partial_melt_defaults, inner_keys=layer_type)

from .partialmelt import PartialMelt
from .partialmelt import calc_partial_melting