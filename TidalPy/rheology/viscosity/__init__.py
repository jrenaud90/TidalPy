from . import viscosity_models
from .defaults import liquid_viscosity_defaults, solid_viscosity_defaults
from ...utilities.classes.model.modelUtils import build_model_default_inputs, find_all_models

known_models, known_model_const_args, known_model_live_args = find_all_models(viscosity_models)


def get_solid_viscosity_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args, solid_viscosity_defaults, inner_keys=layer_type)


def get_liquid_viscosity_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args, liquid_viscosity_defaults, inner_keys=layer_type)


from .viscosity import LiquidViscosity, SolidViscosity
