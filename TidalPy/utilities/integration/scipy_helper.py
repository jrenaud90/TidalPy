""" Helper functions to interface with SciPy's integration suite """

from typing import Tuple

import numpy as np

from . import _solve_ivp


def solve_ivp(
    diffeq, time_span: Tuple[float, float], initial_condition: np.ndarray, args: Tuple = None,
    rtol: float = 1.0e-6, atol: float = 1.0e-8, max_step: float = np.inf,
    first_step: float = None, method: str = 'RK45', t_eval: np.ndarray = np.empty((0,), dtype=np.float64)
    ):

    # Solve the ode
    solution = _solve_ivp(
        diffeq, time_span, initial_condition, method=method, t_eval=t_eval, args=args,
        first_step=first_step, max_step=max_step, rtol=rtol, atol=atol
        )

    # Pull out information from SciPy's solution class
    time_domain = solution.t
    y_results = solution.y
    success = solution.success
    message = solution.message

    return time_domain, y_results, success, message
