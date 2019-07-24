from typing import Tuple, Union, Dict

from time import time as timer
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import METHODS
from datetime import datetime

from .. import verbose_level
from ..utilities.conversions import convert_to_hms
from ..initialize import log
from ..exceptions import IntegrationTimeOut, UnknownModelError, TidalPyIntegrationException


def diffeq_wrap(loop_limits: Tuple[float, float],
                other_input_args: tuple = tuple(),
                use_timeout_kill: bool = False, kill_time: int = 1000, integration_run: bool = False,
                use_progress_bar: bool = True, progress_bar_poll_time: float = 1.5):

    init_time = timer()

    time_i = loop_limits[0]
    time_f = loop_limits[-1]
    time_delta = time_f - time_i
    time_max = 0.997 * time_delta

    prev_delta = 1.
    prev_prev_delta = 1.
    prev_perc = 0.

    # Check if the TidalPy settings override printing
    if verbose_level == 0:
        use_progress_bar = False

    def outer_wrap(diffeq):

        def inner_wrap(time: np.ndarray, dependent_variables: Tuple[np.ndarray]):

            current_int_time = timer() - init_time
            if use_timeout_kill and integration_run:
                if current_int_time >= kill_time:
                    raise IntegrationTimeOut

            if use_progress_bar and integration_run and \
                    current_int_time - inner_wrap.clock_old > progress_bar_poll_time:
                if time < time_i:
                    prog_time = time_i
                elif time > time_f:
                    prog_time = time_max
                else:
                    prog_time = time
                percent_done = (prog_time - time_i) / time_delta
                percent_done = np.nan_to_num(percent_done)
                perc_delta = (percent_done - inner_wrap.prev_perc) / (current_int_time - inner_wrap.clock_old)

                # Use an average of the last 3 (arbitrary) time_left's.
                time_left = 3. * (1. - percent_done) / \
                            (perc_delta + inner_wrap.prev_delta + inner_wrap.prev_prev_delta)
                if time_left > 99999999.0:
                    time_left = 99999999.0
                print('\rPercent Done: {:0>5.2f}%. Approx. Time Remaining: {:0>2} days, '
                      '{:0>2.0f}:{:0>2.0f}::{:0>4.1f}'.format(100. * percent_done, *convert_to_hms(time_left)),
                      flush=True, end='')
                inner_wrap.prev_perc = percent_done
                inner_wrap.clock_old = current_int_time
                inner_wrap.prev_prev_delta = max(inner_wrap.prev_delta, 1.e-5)
                inner_wrap.prev_delta = max(perc_delta, 1.e-5)

            # Main diffeq Call
            result = diffeq(time, dependent_variables, *other_input_args)

            # Return Results
            num_dependent = len(dependent_variables)
            dependent_var_results = result[:num_dependent]
            aux_results = result[num_dependent:]
            if integration_run:
                return dependent_var_results
            else:
                return aux_results

        inner_wrap.clock_old = init_time
        inner_wrap.prev_delta = prev_delta
        inner_wrap.prev_prev_delta = prev_prev_delta
        inner_wrap.prev_perc = prev_perc

        return inner_wrap

    return outer_wrap


def ivp_integration(diff_eq, time_span: Tuple[float, float], initial_conditions: Tuple[np.ndarray],
                    dependent_variable_names: Tuple[str], aux_data_names: Tuple[str] = None,
                    other_input_args: tuple = tuple(), integration_method: str = 'LSODA',
                    integration_rtol: float = 1.0e-4, return_aux_data: bool = True,
                    suppress_integration_errors: bool = False,
                    use_timeout_kill: bool = False, timeout_time: int = 1000,
                    use_progress_bar: bool = True, progress_poll_time: float = 1.5) \
        -> Union[type(None), Dict[str, np.ndarray]]:
    """ Integrates a provided differential equation and returns dependent variables and auxiliary data

    The differential equation can be njit'd. Indeed, it can increase performance, sometimes significantly.

    Parameters
    ----------
    diff_eq : function
        Differential equation should have the following signature:
            diffeq(time, dependent_variables, *other_input_args, diff_eq: bool)
        The output of the diffeq must be a tuple structured as (*dependent_variables, *aux_data)
    time_span : Tuple[float, float]
        Tuple containing start and stop time of the integration.
            Note that, depending on the integration method, the integrator may look before or after the time span
    initial_conditions : Tuple[np.ndarray]
        Tuple of numpy arrays with the initial dependent_variables at the start time
    dependent_variable_names : Tuple[str]
        Names of the dependent variables (same as initial conditions). Used for data storage
    aux_data_names : Tuple[str]
        Names of the auxiliary data. Used for data storage
    other_input_args : tuple = tuple()
        Optional inputs to the diffeq
    integration_method : str = LSODA
        Name of the scipy solve_ivp solver
    integration_rtol : float = 1.0e-4
        Integration adaptive time step tolerance
    return_aux_data : bool = True
        Determines if the diffeq is called a second time to grab non-diffeq data
    suppress_integration_errors : bool = False
        Determines if errors found during integration are raised
        It is advised that for multi-processor runs set this to True
    use_timeout_kill : bool = False
        Determines if integrator checks for max time and, if exceeded, an error raised
    timeout_time : int = 1000
        The max time, in seconds, used for timeout_kill
    use_progress_bar : bool = True
        Determines if a visual progress bar is used during the integration
    progress_poll_time : float = 1.5
        The time, in seconds, after which the progress bar updates

    Returns
    -------
    integration_result : Dict[str, np.ndarray]

    See Also
    --------
    scipy.integrate.solve_ivp

    """

    # Check for issues
    if type(integration_rtol) not in [int, float]:
        raise TypeError
    if integration_rtol <= 0.:
        raise ValueError
    if integration_rtol < 1.e-10:
        # TODO: Is this call to log safe for multiprocessing?
        log.warn('User provided integration value is very small. Integration times may be long')
    for known_int_method in METHODS:
        if integration_method.lower() == known_int_method.lower():
            # Ensure that capitalization is the same
            integration_method = known_int_method
            break
    else:
        raise UnknownModelError('Integration method not found in the scipy package.')
    if len(dependent_variable_names) != len(initial_conditions):
        raise TidalPyIntegrationException('Number of initial conditions must be the same as the number of '
                                          'dependent data names.')

    # Record information to log
    now = datetime.now()
    now_str = now.strftime('%x at %X')
    # TODO: Again, are these log calls MP safe?
    log('Initializing integration using...', level='info')
    log(f'\tDifferential Equation:  {diff_eq.__name__}', level='info')
    log(f'\tIntegration Method:     {integration_method}', level='info')
    log(f'\tIntegration Tolerance:  {integration_rtol}', level='info')
    log(f'\tIntegration Started at: {now_str}', level='info')

    # Setup the differential equation with wrappers
    int_diffeq = diffeq_wrap(time_span, other_input_args=other_input_args,
                             use_timeout_kill=use_timeout_kill, kill_time=timeout_time, integration_run=True,
                             use_progress_bar=use_progress_bar,
                             progress_bar_poll_time=progress_poll_time)(diff_eq)
    aux_diffeq = diffeq_wrap(time_span, other_input_args=other_input_args,
                             use_timeout_kill=use_timeout_kill, kill_time=timeout_time, integration_run=False,
                             use_progress_bar=False,
                             progress_bar_poll_time=progress_poll_time)(diff_eq)

    # Main integration
    log('Integration Starting...', level='info')
    output = None
    try:
        # Main integration
        solution = solve_ivp(int_diffeq, time_span, initial_conditions,
                             method=integration_method, vectorized=True, rtol=integration_rtol)
    except IntegrationTimeOut:
        log('Integration was stopped short due to a forced timeout', level='debug')
        success = False
    except Exception as e:
        success = False
        if not suppress_integration_errors:
            raise e
        else:
            log(f'Integration failed due to an exception: {e}', level='warning')
    else:
        success = solution.success
        output = dict()
        if not success:
            log('WARNING - Integration may not have been preformed correctly.', level='info')
            log(f'ODEINT MSG: {solution.message}', level='info')
        else:
            log('No obvious issues with integration found!', level='info')

        # Pull out dependent variables
        if len(solution.t) > 50000:
            log('Solution data size is very large. Reducing to avoid memory errors in auxiliary grab.', info='debug')
            solution_time_domain = np.linspace(solution.t[0], solution.t[-1], 50000)
            solution_y = np.asarray([np.interp(solution_time_domain, solution.t, solution.y[i, :])
                                     for i in range(solution.y.shape[0])])
        else:
            solution_time_domain = solution.t
            solution_y = solution.y
        del solution
        solution_y_list = [solution_y[i, :] for i in range(solution_y.shape[0])]

        output['time_domain'] = solution_time_domain
        for dependent_data, dependent_data_name in zip(solution_y_list, dependent_variable_names):
            output[dependent_data_name] = dependent_data

        # Calculate auxiliary data
        if return_aux_data and aux_data_names is not None:
            # Calculate auxiliary data
            aux_data_results = aux_diffeq(solution_time_domain, solution_y)
            for aux_data, aux_data_name in zip(aux_data_results, aux_data_names):
                output[aux_data_name] = aux_data

    return output
