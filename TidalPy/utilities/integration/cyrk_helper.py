""" Helper functions to interface with CyRK's integration suite """

from typing import Tuple

import numpy as np

from TidalPy.utilities.performance import njit

from . import _cyrk_ode, _nbrk_ode, _cy2nb, _nb2cy


def cyrk_solver(
    diffeq, time_span: Tuple[float, float], initial_condition: np.ndarray, args: Tuple = tuple(),
    rtol: float = 1.0e-6, atol: float = 1.0e-8, max_step: float = np.inf,
    first_step: float = 0., method: int = 1, t_eval: np.ndarray = np.empty((0,), dtype=np.float64),
    convert_func_to_cy: bool = False
    ):

    # Change the diffeq to match the desired format
    if convert_func_to_cy:
        diffeq = _nb2cy(diffeq, use_njit=True, cache_njit=True)

    return _cyrk_ode(
        diffeq, time_span, initial_condition, args=args,
        rtol=rtol, atol=atol, max_step=max_step, first_step=first_step,
        rk_method=method, t_eval=t_eval
        )


@njit(cacheable=True)
def nbrk_solver(
    diffeq, time_span: Tuple[float, float], initial_condition: np.ndarray, args: Tuple = None,
    rtol: float = 1.0e-6, atol: float = 1.0e-8, max_step: float = np.inf,
    first_step: float = None, method: int = 1, t_eval: np.ndarray = np.empty((0,), dtype=np.float64)
    ):

    return _nbrk_ode(
        diffeq, time_span, initial_condition, args,
        rtol=rtol, atol=atol, max_step=max_step, first_step=first_step,
        rk_method=method, t_eval=t_eval
        )
