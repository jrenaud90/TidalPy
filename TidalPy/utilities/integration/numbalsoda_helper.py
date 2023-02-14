from typing import Tuple

import numpy as np

from numba import cfunc, njit
import numba as nb

from . import lsoda_sig, lsoda, ns_dop853


def numbalsoda_solver(
        diffeq, time_span: Tuple[float, float], initial_condition: np.ndarray, args: Tuple = None,
        rtol: float = 1.0e-6, atol: float = 1.0e-8, max_step: float = np.inf,
        first_step: float = None, method: int = 1, t_eval: np.ndarray = None
        ):

    if t_eval is None:
        raise ValueError('t_eval required for NumbaLSOSA')

    shape = initial_condition.shape

    @cfunc(lsoda_sig)
    def new_diffeq(t, u, du, p):
        u_ = nb.carray(u, shape)
        out = diffeq(t, u_, *args)
        for i in range(u_.size):
            du[i] = out[i]

    new_diffeq_pointer = new_diffeq.address

    solution, success = lsoda(new_diffeq_pointer, initial_condition, t_eval, rtol=rtol, atol=atol)

    message = 'None'

    # Convert solution to transpose
    solution = solution.T

    return t_eval, solution, success, message
