""" A Numba.njit save integration tool utilizing a Explicit Runge-Kutta integration scheme.

    These are largely copied and modified from SciPy's solve_ivp integration function
    (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)

    References
    ----------
    Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski,
        Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson,
        K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey,
        İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman,
        Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa,
        Paul van Mulbregt, and SciPy 1.0 Contributors. (2020)
        SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.
    J. R. Dormand, P. J. Prince, “A family of embedded Runge-Kutta formulae”,
        Journal of Computational and Applied Mathematics, Vol. 6, No. 1, pp. 19-26, 1980.
    L. W. Shampine, “Some Practical Runge-Kutta Formulas”,
        Mathematics of Computation, Vol. 46, No. 173, pp. 135-150, 1986.
    P. Bogacki, L.F. Shampine, “A 3(2) Pair of Runge-Kutta Formulas”,
        Appl. Math. Lett. Vol. 2, No. 4. pp. 321-325, 1989.
"""

from typing import Tuple

import numpy as np

from .helpers import select_initial_step, norm
from .methods.rk import rk_stepper, n_stages_RK23, n_stages_RK45, error_estimator_order_RK23, error_estimator_order_RK45
from ..performance import njit

EPS = np.finfo(float).eps


@njit(cacheable=True)
def select_initial_step(func, t0, y0, f0, direction, order, rtol, atol, args: tuple = tuple()):
    """Empirically select a good initial step.
    The algorithm is described in [1]_.
    Parameters
    ----------
    func : callable
        Right-hand side of the system.
    t0 : float
        Initial value of the independent variable.
    y0 : ndarray, shape (n,)
        Initial value of the dependent variable.
    f0 : ndarray, shape (n,)
        Initial value of the derivative, i.e., ``fun(t0, y0)``.
    direction : float
        Integration direction.
    order : float
        Error estimator order. It means that the error controlled by the
        algorithm is proportional to ``step_size ** (order + 1)`.
    rtol : float
        Desired relative tolerance.
    atol : float
        Desired absolute tolerance.
    Returns
    -------
    h_abs : float
        Absolute value of the suggested initial step.
    References
    ----------
    .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
           Equations I: Nonstiff Problems", Sec. II.4.
    """
    if y0.size == 0:
        return np.inf

    scale = atol + np.abs(y0) * rtol
    d0 = norm(y0 / scale)
    d1 = norm(f0 / scale)
    if d0 < 1e-5 or d1 < 1e-5:
        h0 = 1e-6
    else:
        h0 = 0.01 * d0 / d1

    y1 = y0 + h0 * direction * f0
    f1 = func(t0 + h0 * direction, y1, *args)
    d2 = norm((f1 - f0) / scale) / h0

    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2))**(1 / (order + 1))

    return min(100 * h0, h1)


@njit(cacheable=True)
def rk_integrator(func: callable, t_span: Tuple[float, float], y0: np.ndarray,
                  args: tuple = tuple(),
                  rk_method: int = 0, t_eval_N: int = 3, t_eval_log: bool = False, use_teval: bool = False,
                  rtol=1.0e-3, atol=1.0e-6,
                  verbose: bool = True):
    # Determine integration domain and shape of y
    t0, tf = float(t_span[0]), float(t_span[1])
    if t_eval_log:
        # TODO: Numba does not support np.logspace currently. So this is a hack on logspace
        t_eval = 10**np.linspace(np.log10(t0), np.log10(tf), t_eval_N)
    else:
        t_eval = np.linspace(t0, tf, t_eval_N)
    y_size = y0.size
    dtype = y0.dtype

    # Pull out RK method info
    # Explicit Runge-Kutta Method of order 3(2)
    error_estimator_order = error_estimator_order_RK23
    n_stages = n_stages_RK23
    if rk_method == 0:
        # Explicit Runge-Kutta Method of order 3(2)
        #    These are already set above.
        pass
    elif rk_method == 1:
        # Explicit Runge-Kutta Method of order 5(4)
        n_stages = n_stages_RK45
        error_estimator_order = error_estimator_order_RK45
    elif rk_method == 2:
        # Explicit Runge-Kutta method of order 8.
        #     This is a Python implementation of "DOP853"
        raise NotImplementedError
    else:
        raise Exception

    # Validate tolerances
    if rtol < 100 * EPS:
        if verbose:
            print("`rtol` is too low, setting to 100 * EPS")
        rtol = 100 * EPS

    atol = np.asarray(atol)
    ts = [t0]
    ys = [y0]
    # The extra tuple inside the asarray makes ys_array 2D instead of 1D.
    # ys_array = y0
    # ys_array = np.asarray((y0,), dtype=y0.dtype)
    ys_array = y0.reshape(1, y_size)

    # Determine integration parameters
    max_step = np.abs(tf - t0) / 2.

    # Initialize integration variables
    t_old = 0.
    previous_step_size = 0.
    t = t0
    y = y0
    current_deriv_val = func(t, y, *args)
    if tf != t0:
        direction = np.sign(tf - t0)
    else:
        direction = 1.

    # Determine size of initial step
    step_size_abs = select_initial_step(func, t0, y0, current_deriv_val, direction, error_estimator_order, rtol, atol,
                                        args=args)

    # Setup the K variable (Storage array for putting RK stages).
    K = np.empty((n_stages + 1, y_size), dtype=dtype)

    # While integrator is running the status is 100.
    message = 'Integration completed successfully'
    status = 100
    while status == 100:
        # While step is being calculated its status is 100.
        step_status = 100
        if (y_size == 0) or (t == tf):
            # Integration is complete
            status = 0
            t_old = t
            t = tf
            step_status = 0
        else:
            # Integration is not complete.
            # Run the step function
            t_new, y_new, new_deriv_val, new_step_size_abs, y_old, t_old, old_step, step_status, message = \
                rk_stepper(func, t, y, current_deriv_val, tf, direction, step_size_abs, max_step, atol, rtol, K,
                           rk_method=rk_method, args=args)

            # Check out integration status from step response.
            t = t_new
            y = y_new
            current_deriv_val = new_deriv_val
            step_size_abs = new_step_size_abs
            # Check out integration status from step response.
            if step_status == 0:
                # The step solver has finished, check to see if the full integration is complete
                if direction * (t - tf) >= 0:
                    # Integration is complete.
                    status = 0
            elif step_status < 0:
                # There was an issue with the step solver
                status = step_status
                message = f'There was a problem with the integration.'
                break

        ts.append(t)
        ys.append(y)

        # Numba does not support np.stack(x) if x is a list. So we have to continuously hstack as we go.
        y_new_array = y.reshape(1, y_size)
        ys_array = np.concatenate((ys_array, y_new_array))

    # Clean up output
    success = status >= 0
    ts = np.array(ts)

    # To match the format that scipy follows, we will take the transpose of y.
    ys_array = ys_array.T

    # FIX NOT
    if use_teval:
        t_out = t_eval
        # TODO Numba's version of np.interp does not allow for fp to be > 1D, so we have to
        y_out = np.empty((y_size, t_eval_N), dtype=dtype)
        for y_i in range(y_size):
            y_out[y_i, :] = np.interp(t_eval, ts, ys_array[y_i, :])

    else:
        t_out = ts
        y_out = ys_array

    return t_out, y_out, status, message, success
