from . import melting_models as partial_melting_models
from .defaults import partial_melt_defaults
from ...utilities.classes.model.modelUtils import build_model_default_inputs, find_all_models

known_models, known_model_const_args, known_model_live_args = find_all_models(partial_melting_models)


def get_partial_melt_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args, partial_melt_defaults, inner_keys=layer_type)


from .partialmelt import PartialMelt
from .partialmelt import calculate_melt_fraction, calculate_temperature_frommelt, \
    calculate_temperature_frommelt_array, calculate_melt_fraction_array
