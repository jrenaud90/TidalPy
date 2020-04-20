import os
import copy
from typing import Tuple, List

import numpy as np

from TidalPy.initialize import log
from TidalPy.io import inner_save_dir
from .integration import ivp_integration


def multi_rheo_integration(diff_eq, time_span: Tuple[float, float], initial_conditions: Tuple[np.ndarray],
                           dependent_variable_names: Tuple[str],
                           compliance_func_list: list, compliance_input_list: list, compliance_names: List[str],
                           alternative_save_dir: str = None,
                           **ivp_kwargs):

    universal_input_args = ivp_kwargs.get('other_input_args', tuple())
    all_results = dict()
    rheo_success_track = dict()
    for comp_func, comp_input, comp_name in zip(compliance_func_list, compliance_input_list, compliance_names):
        log(f'Starting integration on Rheology: {comp_name}')
        integration_input_args = (comp_func, comp_input) + universal_input_args
        integration_input_kwargs = copy.deepcopy(ivp_kwargs)
        integration_input_kwargs['other_input_args'] = integration_input_args
        if alternative_save_dir is not None:
            rheo_dir = os.path.join(alternative_save_dir, comp_name.lower())
        else:
            rheo_dir = os.path.join(inner_save_dir, comp_name.lower())
        integration_input_kwargs['alternative_save_dir'] = rheo_dir

        rheo_success, rheo_result = ivp_integration(diff_eq, time_span, initial_conditions, dependent_variable_names,
                                                    **integration_input_kwargs)
        all_results[comp_name] = rheo_result
        rheo_success_track[comp_name] = rheo_success

    return rheo_success_track, all_results
