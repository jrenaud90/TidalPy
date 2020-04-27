import os
import gc
from typing import Tuple
import numpy as np
from time import time as timer

from pathos.multiprocessing import ProcessingPool as Pool

from TidalPy.exceptions import TidalPyIntegrationException
from TidalPy.initialize import log
from TidalPy.utilities.io.io import inner_save_dir, unique_path
from TidalPy.graphics.grid_plot import success_grid_plot
from .integration import ivp_integration
from .success import calc_success_index
from TidalPy.utilities.types import TupleNone

MAX_PROCS_DUE_TO_RAM = 6

def _mp_wrap(mp_input: Tuple[int, str, tuple, dict, TupleNone]):
    """ Wrapper for multiprocessing

    Designed to execute an integration run within a multi-processing scheme.

    Parameters
    ----------
    mp_input : Tuple[int, str, tuple, dict, TupleNone]
        Tuple of various inputs sent to the integrator:
            * run_i: <int> Run number, printed to console.
            * run_dir: <str> Directory to save data and figures.
            * solve_ivp_args: <tuple> args required by ivp_integration function
            * solve_ivp_kwargs: <dict> kwargs required by ivp_integration function
            * success_check_params: <tuple> various inputs used for calc_success_index

    Returns
    -------
    run_i : int
        Run number
    success :
        An indication to if the integration was successful
    integration_result : Dict[str: np.ndarray]
        Result of integration
    """

    # Pull out information for integration
    run_i, run_dir, solve_ivp_args, solve_ivp_kwargs, success_check_params = mp_input
    solve_ivp_kwargs['alternative_save_dir'] = run_dir

    print(f'Working on mp run {run_i}')
    run_init_time = timer()
    integration_success, integration_result = ivp_integration(*solve_ivp_args, **solve_ivp_kwargs)
    print(f'Completed run {run_i}, taking {timer() - run_init_time:.1f} seconds.')

    if success_check_params is not None:
        success = calc_success_index(integration_result, integration_success, *success_check_params)
    else:
        success = integration_success

    return run_i, success, integration_result


def solve_ivp_mp_ic(diff_eq, time_span: Tuple[float, float], initial_conditions: Tuple[np.ndarray],
                    dependent_variable_names: Tuple[str],
                    mp_ic_name_1: str, mp_ic_domain_1: np.ndarray, mp_ic_index_1: int, mp_ic_log_1: bool = False,
                    mp_ic_name_2: str = None, mp_ic_domain_2: np.ndarray = None, mp_ic_index_2: int = None, mp_ic_log_2: bool = False,
                    timeout: float = 300., max_procs: int = None,
                    show_success_plot: bool = False, success_plot_title: str = 'Success Plot',
                    variables_to_check: Tuple[str] = None, expected_values: Tuple[float] = None,
                    time_var_name: str = None, time_to_check_at: float = None,
                    **ivp_integrator_kwargs):

    log('Setting up Multi-Processor Integration Run')
    # Setup Directory
    mp_dir = unique_path(os.path.join(inner_save_dir, 'MultiProcessor-Run'), is_dir=True)
    runs_dir = os.path.join(mp_dir, 'mp_runs')
    if not os.path.exists(runs_dir):
        os.makedirs(runs_dir)

    # Determine phase space of multi-processor integration
    if mp_ic_name_1 not in dependent_variable_names:
        raise TidalPyIntegrationException('Multiprocessor names must match dependent variable names.')
    if mp_ic_name_2 is not None:
        if mp_ic_name_2 not in dependent_variable_names:
            raise TidalPyIntegrationException('Multiprocessor names must match dependent variable names.')
    initial_condition_storage = list()
    mp_name_storage = list()
    run_dir_storage = list()
    run_ic_storage = list()
    run_icnum_storage = list()
    run_i = 0
    for icnum_1, ic_1 in enumerate(mp_ic_domain_1):
        if mp_ic_domain_2 is None:
            temp_ic = list(initial_conditions)
            temp_ic[mp_ic_index_1] = ic_1
            initial_condition_storage.append(tuple(temp_ic))
            run_dir = os.path.join(runs_dir, f'Run_{run_i}')
            run_dir_storage.append(run_dir)
            if not os.path.exists(run_dir):
                os.makedirs(run_dir)
            with open(os.path.join(run_dir, 'run.dat'), 'w') as run_datfile:
                run_datfile.write(f'# MP Run {run_i}\n\n{mp_ic_name_1}\n{ic_1}\n')
            mp_name_storage.append(f'MP Run {run_i}. {mp_ic_name_1} = {ic_1:0.4E}.')
            run_ic_storage.append(ic_1)
            run_icnum_storage.append(icnum_1)
            run_i += 1
        else:
            for icnum_2, ic_2 in enumerate(mp_ic_domain_2):
                temp_ic = list(initial_conditions)
                temp_ic[mp_ic_index_1] = ic_1
                temp_ic[mp_ic_index_2] = ic_2
                initial_condition_storage.append(tuple(temp_ic))
                run_dir = os.path.join(runs_dir, f'Run_{run_i}')
                run_dir_storage.append(run_dir)
                if not os.path.exists(run_dir):
                    os.makedirs(run_dir)
                with open(os.path.join(run_dir, 'run.dat'), 'w') as run_datfile:
                    run_datfile.write(f'# MP Run {run_i}\n\n{mp_ic_name_1}, {mp_ic_name_2}\n{ic_1}, {ic_2}\n')
                mp_name_storage.append(f'MP Run {run_i}. {mp_ic_name_1} = {ic_1:0.4E}; {mp_ic_name_2} = {ic_2:0.4E}.')
                run_ic_storage.append((ic_1, ic_2))
                run_icnum_storage.append((icnum_1, icnum_2))
                run_i += 1

    # Build mp_run input
    ivp_integrator_kwargs['save_data'] = True
    ivp_integrator_kwargs['use_timeout_kill'] = True
    ivp_integrator_kwargs['timeout_time'] = timeout
    ivp_integrator_kwargs['use_progress_bar'] = False
    ivp_integrator_kwargs['suppress_logging'] = True
    mp_run_inputs = list()

    if variables_to_check is None or expected_values is None or time_var_name is None or time_to_check_at is None:
        # User did not provide enough information to perform a success check
        success_check_params = None
    else:
        success_check_params = (variables_to_check, expected_values, time_var_name, time_to_check_at)

    for run_i, init_cond in enumerate(initial_condition_storage):
        mp_run_input = (
            run_i,
            run_dir_storage[run_i],
            (diff_eq, time_span, init_cond, dependent_variable_names),
            ivp_integrator_kwargs,
            success_check_params
        )
        mp_run_inputs.append(mp_run_input)

    if max_procs is None:
        max_procs = max(1, os.cpu_count() - 1)
    max_procs = min(MAX_PROCS_DUE_TO_RAM, max_procs)

    mp_length = len(mp_run_inputs)
    log(f'Phase Space Length = {mp_length}')
    chunksize, extra = divmod(mp_length, max_procs)
    if extra:
        chunksize += 1

    # Main MP Call
    log(f'Processors = \t{max_procs} (of {os.cpu_count()})')
    log(f'Chunksize =  \t{chunksize}')
    log(f'Timeout =    \t{timeout:0.1f}')
    log(f'Starting Pool...')
    success_by_run = dict()
    result_by_run = dict()
    with Pool(max_procs) as mp_pool:
        mp_output = mp_pool.uimap(_mp_wrap, mp_run_inputs, chunksize=chunksize)
        for run_i, integration_success, integration_result in mp_output:
            success_by_run[run_i] = integration_success
            result_by_run[run_i] = integration_result
    del mp_output, mp_pool
    gc.collect()

    # Store results and put in a more readable container
    if mp_ic_name_2 is None:
        ic2_len = 1
    else:
        ic2_len = len(mp_ic_domain_2)
    success_by_ic = np.zeros((len(mp_ic_domain_1), ic2_len))
    run_to_index = dict()
    for run_i, integration_success in success_by_run.items():
        run_info = run_icnum_storage[run_i]
        if type(run_info) == tuple:
            icnum_1, icnum_2 = run_info
            success_by_ic[icnum_1, icnum_2] = integration_success
            run_to_index[run_i] = (icnum_1, icnum_2)
        else:
            icnum_1 = run_info
            success_by_ic[icnum_1, 0] = integration_success
            run_to_index[run_i] = (icnum_1, 0)

    # Plot success grid
    if show_success_plot:
        if mp_ic_domain_2 is None:
            mp_ic_domain_2 = np.asarray([0])
        fig = success_grid_plot({success_plot_title: success_by_ic}, mp_ic_domain_1, mp_ic_domain_2,
                                xname=mp_ic_name_1, yname=mp_ic_name_2, xlog=mp_ic_log_1, ylog=mp_ic_log_2,
                                auto_show=True)
        fig_path = os.path.join(mp_dir, 'success_plot.pdf')
        log(f'Saving success figure to {fig_path}')
        fig.savefig(fig_path)

    return success_by_run, result_by_run, success_by_ic, mp_dir, run_to_index

