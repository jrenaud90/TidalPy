from typing import Tuple

import numpy as np

from . import lsoda, lsoda_sig
from ..performance import njit

from numba import cfunc

from numbalsoda import dop853

def numbalsoda_solver(
    diffeq, time_span: Tuple[float, float], initial_condition: np.ndarray, args: Tuple = None,
    rtol: float = 1.0e-6, atol: float = 1.0e-8, max_step: float = np.inf,
    first_step: float = None, method: int = 1, t_eval: np.ndarray = None):

    if t_eval is None:
        raise ValueError('t_eval required for NumbaLSOSA')

    y_size = initial_condition.size

    # Set arg array to none since they will be passed to the function during construction below.
    arg_array = np.empty(0, dtype=np.float64)

    # Convert u from complex128 to float64
    initial_condition_in = np.empty(y_size*2, dtype=np.float64)
    for i in range(y_size):
        # Store real part
        initial_condition_in[i * 2] = np.real(initial_condition[i])
        # Store imag part
        initial_condition_in[(i * 2) + 1] = np.imag(initial_condition[i])

    # Need to modify diffeq before we build wrapper
    @njit
    def float_diffeq(t, y):

        # Convert float64 -> complex128
        y_ = np.empty(y_size, dtype=np.complex128)
        for i in range(y_size):
            real_y = y[i * 2]
            imag_y = y[(i * 2) + 1]

            y_[i] = real_y + 1.0j * imag_y

        # run diffeq
        dy = np.asarray(diffeq(t, y_, *args))

        # convert back
        dy_ = np.empty(y_size * 2, dtype=np.float64)
        for i in range(y_size):
            # Store real part
            dy_[i * 2] = np.real(dy[i])
            # Store imag part
            dy_[(i * 2) + 1] = np.imag(dy[i])

        return dy_


    @cfunc(lsoda_sig)
    def new_diffeq(t, u, du, p):
        du = float_diffeq(t, u)

    new_diffeq_pointer = new_diffeq.address

    solution, success = dop853(new_diffeq_pointer, initial_condition_in, t_eval, data=arg_array, rtol=rtol, atol=atol)

    # message = solution.message
    message = 'None'
    # t_ = solution.t

    # Convert u from float64 to complex128
    y_ = np.empty((y_size, t_eval.size), dtype=np.complex128)
    for i in range(y_size):
        # Also take transpose at the same time.
        real_y = solution[:, i * 2]
        imag_y = solution[:, (i * 2) + 1]

        y_[i, :] = real_y + 1.0j * imag_y

    import pdb; pdb.set_trace()
    return t_eval, y_, success, message
