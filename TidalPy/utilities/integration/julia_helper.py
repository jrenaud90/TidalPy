""" Helper functions to interface with Julia / Diffeqpy's integration suite """

from typing import Tuple

import numpy as np

from . import _de, _ode, julia_installed

# Read more about Julia's ode solvers here: https://diffeq.sciml.ai/dev/solvers/ode_solve/

known_integration_methods = (
    'rk4', 'rk45', 'tsit5', 'rko65', 'tsitpap8', 'feagin10', 'feagin12', 'feagin14', 'bs3', 'bs5', 'vern6', 'vern7',
    'vern8', 'vern9', 'kuttaprk2p5', 'rosenbrock23', 'rodas4', 'rodas5'
    )


def get_julia_solver(solver_name: str):
    """ Find a ODE solver in the Julia diffeq package.

    Read more about Julia's ode solvers here: https://diffeq.sciml.ai/dev/solvers/ode_solve/

    Parameters
    ----------
    solver_name : str
        Name of the Julia ode solver

    Returns
    -------
    ode_system :
    solver :

    """

    if not julia_installed:
        raise ImportError('Julia package not found. Can not load ODE solver.')

    non_stiff_solvers = {
        # The canonical Runge-Kutta Order 4 method. Uses a defect control for adaptive stepping using maximum error over the whole interval.
        'rk4'        : _ode.RK4,
        'rk45'       : _ode.RK4,
        # Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).
        'tsit5'      : _ode.Tsit5,
        # Tsitouras' Runge-Kutta-Oliver 6 stage 5th order method. This method is robust on problems which have a singularity at t=0.
        'rko65'      : _ode.RKO65,
        # Tsitouras-Papakostas 8/7 Runge-Kutta method.
        'tsitpap8'   : _ode.TsitPap8,
        # Feagin's 10th-order Runge-Kutta method.
        'feagin10'   : _ode.Feagin10,
        # Feagin's 12th-order Runge-Kutta method.
        'feagin12'   : _ode.Feagin12,
        # Feagin's 14th-order Runge-Kutta method.
        'feagin14'   : _ode.Feagin14,
        # Additionally, the following algorithms have a lazy interpolant:
        # BS5 - Bogacki-Shampine 3/2 Runge-Kutta method. (lazy 5th order interpolant).
        'bs3'        : _ode.BS3,
        # BS5 - Bogacki-Shampine 5/4 Runge-Kutta method. (lazy 5th order interpolant).
        'bs5'        : _ode.BS5,
        # Verner's "Most Efficient" 6/5 Runge-Kutta method. (lazy 6th order interpolant).
        'vern6'      : _ode.Vern6,
        # Verner's "Most Efficient" 7/6 Runge-Kutta method. (lazy 7th order interpolant).
        'vern7'      : _ode.Vern7,
        # Verner's "Most Efficient" 8/7 Runge-Kutta method. (lazy 8th order interpolant).
        'vern8'      : _ode.Vern8,
        # Verner's "Most Efficient" 9/8 Runge-Kutta method. (lazy 9th order interpolant).
        'vern9'      : _ode.Vern9,
        # A 5 parallel, 2 processor explicit Runge-Kutta method of 5th order.
        # These methods utilize multithreading on the f calls to parallelize the problem. This requires that simultaneous calls to f are thread-safe.
        'kuttaprk2p5': _ode.KuttaPRK2p5
        }

    stiff_solvers = {
        'rosenbrock23': _de.Rosenbrock23,
        # The ODEInterface algorithms are the classic Fortran algorithms. While the non-stiff algorithms are superseded by the more featured and higher performance Julia implementations from OrdinaryDiffEq.jl, the stiff solvers such as radau are some of the most efficient methods available (but are restricted for use on arrays of Float64).
        # Rosenbrock 4(3) method.
        'rodas4'      : _de.Rodas4,
        'rodas5'      : _de.Rodas5
        }

    try:
        solver = non_stiff_solvers[solver_name.lower()]
        ode_system = _ode
    except KeyError:
        try:
            solver = stiff_solvers[solver_name.lower()]
            ode_system = _de
        except KeyError:
            raise KeyError(f'Unknown Julia Integration Model: {solver_name}.')

    def julia_integrator(
        diffeq, time_span: Tuple[float, float], initial_condition: np.ndarray, args: Tuple = None,
        rtol: float = 1.0e-6, atol: float = 1.0e-8, max_step: float = np.inf,
        first_step: float = 0., method: int = 1, t_eval: np.ndarray = np.empty((0,), dtype=np.float64)
        ):

        # Some inputs are unused. Input structure is kept for consistency
        del method, max_step, first_step

        # Change the diffeq to match the desired format
        def diffeq_julia(u, p, r):
            # Julia integrator flips the order of the variables for the differential equation.
            output = diffeq(r, u, *p)
            return list(output)

        # Setup Julia ODE problem
        problem = ode_system.ODEProblem(diffeq_julia, initial_condition, time_span, args)

        # Solve the ode
        solution = ode_system.solve(problem, solver(), abstol=atol, reltol=rtol)

        # Find integration codes
        message = solution.retcode
        success = message.lower() in ['default', 'success']

        if success:
            # Pull out y values and tranpose so they have the same shape as scipys
            y_results = np.transpose(solution.u)
            time_domain = solution.t

            if t_eval is not None:
                if t_eval.size > 1:
                    # Julia does not have the same t_eval. There is the "saveat" keyword but can cause issues.
                    #    So perform an interpolation for the desired radii
                    y_results_reduced = np.zeros((y_results.shape[0], t_eval.size), dtype=initial_condition.dtype)
                    for i in range(y_results.shape[0]):
                        y_results_reduced[i, :] = np.interp(t_eval, time_domain, y_results[i, :])
                    time_domain = t_eval
                    y_results = y_results_reduced
        else:
            time_domain = None
            y_results = None

        return time_domain, y_results, success, message

    return julia_integrator
