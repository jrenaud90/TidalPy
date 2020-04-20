""" Tools to calculate a success index for a given integration simulation.

"""
from typing import Tuple, Dict

import numpy as np

from ..utilities.numpy_help import find_nearest, normalize_dict

def calc_success_index(data_dict: Dict[str, np.ndarray], integration_success: int,
                       variables_to_check: Tuple[str], expected_values: Tuple[float],
                       time_var_name: str, time_to_check_at: float):
    """ Calculates the success index of multiple variables at a specific time and returns a single value as an
    estimation of that simulation's success at matching expectations.

    Parameters
    ----------
    data_dict : Dict[str, np.ndarray]
        Dictionary of data produced by an integration simulation.
    integration_success : int
        Output of the integrator used to override the success index if the integration timeout'd or failed for another
        reason.
    variables_to_check : Tuple[str]
        List-like object of strings of variables to test.
        The strings correspond to the keys of the variables in the data_dict.
    expected_values : Tuple[float]
        List-like object of expected values for each test variable.
        Must be the same length and order as variables_to_check.
    time_var_name : str
        Name of the key for the time variable stored in data_dict.
    time_to_check_at : float
        Time at which to check the variables.
        Must be given in the same dimensional units as the array provided with time_var_name.

    Returns
    -------
    total_perc_diff : float
        Numerical measurement of this simulation's ability to match expectations.

    """

    # Check input
    if type(variables_to_check) not in [list, tuple]:
        variables_to_check = (variables_to_check, )
    if type(expected_values) not in [list, tuple]:
        expected_values = (expected_values, )
    if len(variables_to_check) != len(expected_values):
        raise ValueError

    if integration_success < 0:
        return integration_success

    # Find index at which to check (this will be the same for all variables).
    time_array = data_dict[time_var_name]
    check_index = find_nearest(time_array, time_to_check_at)

    max_perc_diff = len(variables_to_check)
    total_perc_diff = 0.
    for var_name, expected_value in zip(variables_to_check, expected_values):

        experimental_value = data_dict[var_name][check_index]
        # TODO: Is this the correct way to handle these statistics
        if expected_value == 0.:
            perc_diff = abs(experimental_value)
        else:
            perc_diff = abs((expected_value - experimental_value) / expected_value)

        total_perc_diff += perc_diff

    success = max_perc_diff - total_perc_diff
    if success < 0.:
        success = 0.
    return success


def calc_success_index_multirheo(data_dict_byrheo: Dict[str, Dict[str, np.ndarray]],
                                 integration_success_byrheo: Dict[str, int],
                                 variables_to_check: Tuple[str], expected_values: Tuple[float],
                                 time_var_name: str, time_to_check_at: float, normalize_results: bool = False,
                                 **normalize_kwargs):
    """ Calculates the success index of multiple variables, for multiple rheologies, at a specific time and
    returns a single value as an estimation of that rheology simulation's success at matching expectations.

    Parameters
    ----------
    data_dict_byrheo : Dict[str, Dict[str, np.ndarray]]
        Dictionary of data produced by an integration simulation.
        Stored as rheology_name: data_dict where data_dict is stored as variable_name: variable_data
    integration_success_byrheo : Dict[str, int]
        Dictionary of integrator success results for each rheology.
        Output of the integrator used to override the success index if the integration timeout'd or failed for another
        reason.
    variables_to_check : Tuple[str]
        List-like object of strings of variables to test.
        The strings correspond to the keys of the variables in the data_dict.
    expected_values : Tuple[float]
        List-like object of expected values for each test variable.
        Must be the same length and order as variables_to_check.
    time_var_name : str
        Name of the key for the time variable stored in data_dict.
    time_to_check_at : float
        Time at which to check the variables.
        Must be given in the same dimensional units as the array provided with time_var_name.

    Returns
    -------
    success_index_by_rheo : Dict[str, float]
        Numerical measurement of this simulation's ability to match expectations for each rheology.

    """

    success_index_by_rheo = dict()
    for rheo_name, data_dict in data_dict_byrheo.items():
        int_success = integration_success_byrheo[rheo_name]
        perc_diff = calc_success_index(data_dict, int_success, variables_to_check, expected_values,
                                       time_var_name, time_to_check_at)
        success_index_by_rheo[rheo_name] = perc_diff


    if normalize_results:
        return normalize_dict(success_index_by_rheo, pass_negatives=True, **normalize_kwargs)
    else:
        return success_index_by_rheo