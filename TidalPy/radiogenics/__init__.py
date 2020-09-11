from . import radiogenic_models
from .defaults import standard_isotope_input, radiogenics_defaults, known_isotope_data
from ..utilities.classes.model.modelUtils import find_all_models, build_model_default_inputs

known_models, known_model_const_args, known_model_live_args = find_all_models(radiogenic_models)

def get_radiogenic_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args, radiogenic_models, inner_keys=layer_type)

from .radiogenics import Radiogenics