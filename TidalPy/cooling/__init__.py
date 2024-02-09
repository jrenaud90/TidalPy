from TidalPy.utilities.classes.model.model_utils import build_model_default_inputs, find_all_models

import TidalPy
from . import cooling_models


parameter_info_loc = ('cooling',)

known_models, known_model_const_args, known_model_live_args = find_all_models(cooling_models)


def get_cooling_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args,
                                      TidalPy.config['layers'],
                                      inner_keys=(layer_type, 'cooling'))


from .cooling import CoolingModel
from .cooling_models import CoolingOutputType
