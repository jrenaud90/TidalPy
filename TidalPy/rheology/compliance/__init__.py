from ...utilities.model_utils import find_all_models, build_model_defaults
from ..defaults import rheology_param_defaults
from . import compliance_models as complex_compliances
from .compliance_models import find_factorial

parameter_info_loc = ('complex_compliance',)

known_models, known_model_const_args, known_model_live_args = \
    find_all_models(complex_compliances, ignore_functional_types=(find_factorial,))

known_model_defaults = build_model_defaults(known_model_const_args, rheology_param_defaults,
                                            inner_keys=parameter_info_loc)
