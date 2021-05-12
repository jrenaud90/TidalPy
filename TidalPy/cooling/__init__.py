from . import cooling_models
from .defaults import cooling_defaults
from ..utilities.classes.model.modelUtils import find_all_models, build_model_default_inputs

parameter_info_loc = ('cooling',)

known_models, known_model_const_args, known_model_live_args = find_all_models(cooling_models)

def get_cooling_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args, cooling_defaults, inner_keys=layer_type)

from .cooling import CoolingModel, CoolingOutputType