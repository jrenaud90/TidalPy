from TidalPy.utilities.classes.model.model_utils import build_model_default_inputs, find_all_models

from . import radiogenic_models

known_models, known_model_const_args, known_model_live_args = find_all_models(radiogenic_models)


def get_radiogenic_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args, radiogenic_models, inner_keys=layer_type)

from .radiogenics import Radiogenics
