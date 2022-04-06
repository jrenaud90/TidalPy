from typing import Tuple

import numpy as np

from ..performance import njit

# Multiply steps computed from asymptotic behaviour of errors by this.
SAFETY = 0.9

MIN_FACTOR = 0.2  # Minimum allowed decrease in a step size.
MAX_FACTOR = 10.  # Maximum allowed increase in a step size.

EPS = np.finfo(np.float64).eps

RK23_order = 3
RK23_error_estimator_order = 2
RK23_n_stages = 3
RK23_C = np.array([0, 1 / 2, 3 / 4], order='C')
RK23_A = np.array(
    [
        [0, 0, 0],
        [1 / 2, 0, 0],
        [0, 3 / 4, 0]
        ], order='C'
    )
RK23_B = np.array([2 / 9, 1 / 3, 4 / 9], order='C')
RK23_E = np.array([5 / 72, -1 / 12, -1 / 9, 1 / 8], order='C')
RK23_P = np.array(
    [[1, -4 / 3, 5 / 9],
     [0, 1, -2 / 3],
     [0, 4 / 3, -8 / 9],
     [0, -1, 1]], order='C'
    )

RK45_order = 5
RK45_error_estimator_order = 4
RK45_n_stages = 6
RK45_C = np.array([0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1], order='C')
RK45_A = np.array(
    [
        [0, 0, 0, 0, 0],
        [1 / 5, 0, 0, 0, 0],
        [3 / 40, 9 / 40, 0, 0, 0],
        [44 / 45, -56 / 15, 32 / 9, 0, 0],
        [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0],
        [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656]
        ], order='C'
    )
RK45_B = np.array([35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84], order='C')
RK45_E = np.array(
    [-71 / 57600, 0, 71 / 16695, -71 / 1920, 17253 / 339200, -22 / 525,
     1 / 40], order='C'
    )
# Corresponds to the optimum value of c_6 from [2]_.
RK45_P = np.array(
    [
        [1, -8048581381 / 2820520608, 8663915743 / 2820520608,
         -12715105075 / 11282082432],
        [0, 0, 0, 0],
        [0, 131558114200 / 32700410799, -68118460800 / 10900136933,
         87487479700 / 32700410799],
        [0, -1754552775 / 470086768, 14199869525 / 1410260304,
         -10690763975 / 1880347072],
        [0, 127303824393 / 49829197408, -318862633887 / 49829197408,
         701980252875 / 199316789632],
        [0, -282668133 / 205662961, 2019193451 / 616988883, -1453857185 / 822651844],
        [0, 40617522 / 29380423, -110615467 / 29380423, 69997945 / 29380423]], order='C'
    )


@njit(cacheable=True)
def _norm(x):
    return np.linalg.norm(x) / np.sqrt(x.size)


@njit(cacheable=True)
def rk_integrate(
    diffeq: callable, t_span: Tuple[float, float], y0: np.ndarray, args: tuple = tuple(),
    rtol: float = 1.e-6, atol: float = 1.e-8,
    max_step: float = np.inf, first_step: float = None,
    rk_method: int = 1, t_eval: np.ndarray = np.empty(0, dtype=np.float64)
    ):
    """ A Numba-safe Rugge-Kutta Integrator based on Scipy's solve_ivp RK integrator.

    Parameters
    ----------
    diffeq : callable
        An njit-compiled function that defines the derivatives of the problem.
    t_span : Tuple[float, float]
        A tuple of the beginning and end of the integration domain's dependent variables.
    y0 : np.ndarray
        1D array of the initial values of the problem at t_span[0]
    args : tuple = tuple()
        Any additional arguments that are passed to dffeq.
    rtol : float = 1.e-6
        Integration relative tolerance used to determine optimal step size.
    atol : float = 1.e-8
        Integration absolute tolerance used to determine optimal step size.
    max_step : float = np.inf
        Maximum allowed step size.
    first_step : float = None
        Initial step size. If `None`, then the function will attempt to determine an appropriate initial step.
    rk_method : int = 1
        The type of RK method used for integration
            0 = RK23
            1 = RK45
            # TODO NotImplemented 2 = DOP
    t_eval : np.ndarray = None
        If provided, then the function will interpolate the integration results to provide them at the
            requested t-steps.

    Returns
    -------
    time_domain : np.ndarray
        The final time domain. This is equal to t_eval if it was provided.
    y_results : np.ndarray
        The solution of the differential equation provided for each time_result.
    success : bool
        Final integration success flag.
    message : str
        Any integration messages, useful if success=False.

    """

    # Clean up and interpret inputs
    t_start = t_span[0]
    t_end = t_span[1]
    direction = np.sign(t_end - t_start) if t_end != t_start else 1
    direction_inf = direction * np.inf
    y0 = np.asarray(y0)
    y_size = y0.size
    dtype = y0.dtype
    time_domain = [t_start]
    y_results = y0.reshape(1, y_size)

    # Integrator Status Codes
    #   0  = Running
    #   -1 = Failed
    #   1  = Finished with no obvious issues
    status = 0

    # Determine RK constants
    if rk_method == 0:
        # RK23 Method
        rk_order = RK23_order
        error_order = RK23_error_estimator_order
        rk_n_stages = RK23_n_stages
        C = RK23_C
        A = RK23_A
        B = RK23_B
        E = RK23_E
        P = RK23_P
    else:
        # RK45 Method
        rk_order = RK45_order
        error_order = RK45_error_estimator_order
        rk_n_stages = RK45_n_stages
        C = RK45_C
        A = RK45_A
        B = RK45_B
        E = RK45_E
        P = RK45_P

    # Recast some constants into the correct dtype so they can be used with y.
    A = np.asarray(A, dtype=dtype)
    B = np.asarray(B, dtype=dtype)
    E = np.asarray(E, dtype=dtype)

    error_expo = 1. / (error_order + 1.)

    # Check tolerances
    if rtol < 100. * EPS:
        rtol = 100. * EPS

    atol = np.asarray(atol)
    if atol.ndim > 0 and atol.shape != (y_size,):
        # atol must be either the same for all y or must be provided as an array, one for each y.
        raise Exception

    # Initialize variables for start of integration
    t_now = t_start
    y_now = np.asarray(y0)
    dydt_now = np.asarray(diffeq(t_now, y_now, *args), dtype=dtype)
    if first_step is None:
        # Select an initial step size based on the differential equation.
        # .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
        #        Equations I: Nonstiff Problems", Sec. II.4.
        if y_size == 0:
            step_size = np.inf
        else:
            scale = atol + np.abs(y_now) * rtol
            d0 = _norm(y_now / scale)
            d1 = _norm(dydt_now / scale)
            if d0 < 1.e-5 or d1 < 1.e-5:
                h0 = 1.e-6
            else:
                h0 = 0.01 * d0 / d1

            y1 = y_now + h0 * direction * dydt_now
            t1 = t_now + h0 * direction

            dydt1 = np.asarray(diffeq(t1, y1, *args), dtype=dtype)
            d2 = _norm((dydt1 - dydt_now) / scale) / h0

            if d1 <= 1.e-15 and d2 <= 1.e-15:
                h1 = max(1.e-6, h0 * 1.e-3)
            else:
                h1 = (0.01 / max(d1, d2))**error_expo

            step_size = min(100. * h0, h1)
    else:
        if first_step <= 0.:
            # Step size must be a positive number
            raise Exception
        elif first_step > np.abs(t_end - t_start):
            # Step size can not exceed bounds
            raise Exception
        step_size = first_step

    # Initialize RK-K variable
    K = np.empty((rk_n_stages + 1, y_size), dtype=dtype)

    # Main integration loop
    message = 'Running...'
    # # Time Loop
    while status == 0:

        if t_now == t_end or y_size == 0:
            t_prev = t_now
            t_now = t_end
            message = 'Finished'
            status = 1
            break

        # Run RK integration step
        # Determine step size based on previous loop
        min_step = 10. * np.abs(np.nextafter(t_now, direction_inf) - t_now)
        # Look for over/undershoots in previous step size
        if step_size > max_step:
            step_size = max_step
        elif step_size < min_step:
            step_size = min_step

        # Determine new step size
        step_accepted = False
        step_rejected = False
        step_error = False
        rejected_message = 'Proper step size not found.'
        # # Step Loop
        while not step_accepted:

            if step_size < min_step:
                step_error = True
                rejected_message = 'Required step size is less than spacing between numbers.'
                break

            # Move time forward for this particular step size
            step = step_size * direction
            t_new = t_now + step

            # Check that we are not at the end of integration with that move
            if direction * (t_new - t_end) > 0:

                t_new = t_end

            # Correct the step if we were at the end of integration
            step = t_new - t_now
            step_size = np.abs(step)

            # Calculate derivative using RK method
            K[0] = dydt_now
            s = 1
            for a, c in zip(A[1:], C[1:]):
                K_ = np.ascontiguousarray(K[:s].T)
                a_ = np.ascontiguousarray(a[:s])
                dy = np.dot(K_, a_) * step
                K[s] = np.asarray(diffeq(t_now + c * step, y_now + dy, *args), dtype=dtype)
                s += 1

            y_new = y_now + step * np.dot(K[:-1].T, B)
            dydt_new = np.asarray(diffeq(t_now + step, y_new, *args), dtype=dtype)
            K[-1] = dydt_new

            # Check how well this step performed
            scale = atol + np.maximum(np.abs(y_now), np.abs(y_new)) * rtol
            error_norm = _norm(np.dot(K.T, E) * step / scale)

            if error_norm < 1.:
                # The error is low! Let's update this step for the next time loop
                if error_norm == 0.:
                    step_factor = MAX_FACTOR
                else:
                    step_factor = min(
                        MAX_FACTOR,
                        SAFETY * error_norm**-error_expo
                        )

                if step_rejected:
                    # There were problems with this step size on the previous step loop. Make sure factor does
                    #    not exasperate them.
                    step_factor = min(step_factor, 1.)

                step_size = step_size * step_factor
                step_accepted = True
            else:
                step_size = step_size * max(
                    MIN_FACTOR,
                    SAFETY * error_norm**-error_expo
                    )
                step_rejected = True

        if not step_accepted or step_error:
            # Issue with step convergence
            status = -1
            message = 'Error in step size calculation:\n' + rejected_message
            break

        # End of step loop. Update the _now variables
        t_old = t_now
        y_old = y_now
        dydt_old = dydt_now

        t_now = t_new
        y_now = y_new
        dydt_now = dydt_new

        # Save data
        time_domain.append(t_now)

        # Numba does not support np.stack(x) if x is a list. So we have to continuously hstack as we go.
        y_new_array = y_now.reshape(1, y_size)
        y_results = np.concatenate((y_results, y_new_array))

    time_domain = np.asarray(time_domain, dtype=np.float64)
    # To match the format that scipy follows, we will take the transpose of y.
    y_results = y_results.T

    if t_eval.size > 0:
        # User only wants data at specific points.
        # The current version of this function has not implemented sicpy's dense output, so we must use an interpolation.
        t_eval = np.asarray(t_eval, dtype=np.float64)
        y_results_reduced = np.empty((y_size, t_eval.size), dtype=dtype)

        for i in range(y_size):
            # np.interp only works on 1D arrays so we must loop through each of the variables:
            y_results_reduced[i, :] = np.interp(t_eval, time_domain, y_results[i, :])

        y_results = y_results_reduced
        time_domain = t_eval

    success = status == 1

    # Make sure arrays are C-contiguous
    y_results = np.ascontiguousarray(y_results)
    time_domain = np.ascontiguousarray(time_domain)

    return time_domain, y_results, success, message
