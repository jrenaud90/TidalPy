from . import compliance_models as complex_compliances
from .compliance_models import find_factorial
from .defaults import complex_compliance_defaults
from ...utilities.classes.model.modelUtils import find_all_models, build_model_default_inputs

known_models, known_model_const_args, known_model_live_args = \
    find_all_models(complex_compliances, ignore_functional_types=(find_factorial,))

def get_complex_comp_model_default_inputs(layer_type: str):
    return build_model_default_inputs(known_model_const_args, complex_compliance_defaults, inner_keys=layer_type)

from .complex_compliance import ComplexCompliance