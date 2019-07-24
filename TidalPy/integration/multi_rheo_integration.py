import copy
from typing import Tuple, List

import numpy as np

from ..initialize import log
from .integration import ivp_integration


def multi_rheo_integration(diff_eq, time_span: Tuple[float, float], initial_conditions: Tuple[np.ndarray],
                           dependent_variable_names: Tuple[str],
                           compliance_func_list: list, compliance_input_list: list, compliance_names: List[str],
                           **ivp_kwargs):

    universal_input_args = ivp_kwargs.get('other_input_args', tuple())
    all_results = dict()
    for comp_func, comp_input, comp_name in zip(compliance_func_list, compliance_input_list, compliance_names):
        log(f'Starting integration on Rheology: {comp_name}')
        integration_input_args = (comp_func, comp_input) + universal_input_args
        integration_input_kwargs = copy.deepcopy(ivp_kwargs)
        integration_input_kwargs['other_input_args'] = integration_input_args
        rheo_result = ivp_integration(diff_eq, time_span, initial_conditions, dependent_variable_names,
                                      **integration_input_kwargs)
        rheo_result[comp_name] = rheo_result

    return all_results
