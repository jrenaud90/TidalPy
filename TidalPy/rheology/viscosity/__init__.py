from TidalPy.utilities.classes.model.model_utils import build_model_default_inputs, find_all_models

import TidalPy
from . import viscosity_models

known_models, known_model_const_args, known_model_live_args = find_all_models(viscosity_models)


def get_solid_viscosity_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args,
                                      TidalPy.config['layers'],
                                      inner_keys=(layer_type, 'solid_viscosity'))


def get_liquid_viscosity_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args,
                                      TidalPy.config['layers'],
                                      inner_keys=(layer_type, 'liquid_viscosity'))


from .viscosity import LiquidViscosity, SolidViscosity
