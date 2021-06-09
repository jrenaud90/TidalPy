""" Functions and constants used during Explicit Runge-Kutta integration.

    These are largely copied and edited from SciPy's solve_ivp integration function
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

import numpy as np

from ..helpers import norm
from ...performance import njit

# Multiply steps computed from asymptotic behaviour of errors by this.
SAFETY = 0.9

MIN_FACTOR = 0.2  # Minimum allowed decrease in a step size.
MAX_FACTOR = 10.  # Maximum allowed increase in a step size.

# Parameters for 3rd order RK
order_RK23 = 3
error_estimator_order_RK23 = 2
n_stages_RK23 = 3
C_RK23 = np.array([0, 1 / 2, 3 / 4])
A_RK23 = np.array([
    [0, 0, 0],
    [1 / 2, 0, 0],
    [0, 3 / 4, 0]
])
B_RK23 = np.array([2 / 9, 1 / 3, 4 / 9])
E_RK23 = np.array([5 / 72, -1 / 12, -1 / 9, 1 / 8])
P_RK23 = np.array([[1, -4 / 3, 5 / 9],
                   [0, 1, -2 / 3],
                   [0, 4 / 3, -8 / 9],
                   [0, -1, 1]])

# Parameters for 5th order RK
order_RK45 = 5
error_estimator_order_RK45 = 4
n_stages_RK45 = 6
C_RK45 = np.array([0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1])
A_RK45 = np.array([
    [0, 0, 0, 0, 0],
    [1 / 5, 0, 0, 0, 0],
    [3 / 40, 9 / 40, 0, 0, 0],
    [44 / 45, -56 / 15, 32 / 9, 0, 0],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656]
])
B_RK45 = np.array([35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84])
E_RK45 = np.array([-71 / 57600, 0, 71 / 16695, -71 / 1920, 17253 / 339200, -22 / 525,
                   1 / 40])
# Corresponds to the optimum value of c_6 from
#    L. W. Shampine, "Some Practical Runge-Kutta Formulas", Mathematics
#    of Computation,, Vol. 46, No. 173, pp. 135-150, 1986
P_RK45 = np.array([
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
    [0, 40617522 / 29380423, -110615467 / 29380423, 69997945 / 29380423]])


# Error estimation for RK23 and RK45 (DOP-8 requires a different method)
@njit(cacheable=True)
def _RK23_45_estimate_error_norm(K, E, step_to_use, scale):
    abs_error = np.dot(K.T, E) * step_to_use
    return norm(abs_error / scale)


@njit(cacheable=True)
def rk_stepper(func, t, y, current_deriv_val, tf, direction, step_size_abs, max_step, atol, rtol, K, rk_method: int = 0,
               args: tuple = tuple()):
    # Explicit Runge-Kutta Method of order 3(2)
    order = order_RK23
    error_estimator_order = error_estimator_order_RK23
    error_func_norm = _RK23_45_estimate_error_norm
    n_stages = n_stages_RK23
    C = C_RK23
    A = A_RK23
    B = B_RK23
    E = E_RK23
    P = P_RK23
    if rk_method == 0:
        # Explicit Runge-Kutta Method of order 3(2)
        #    These are already set above.
        pass
    elif rk_method == 1:
        # Explicit Runge-Kutta Method of order 5(4)
        order = order_RK45
        error_estimator_order = error_estimator_order_RK45
        error_func_norm = _RK23_45_estimate_error_norm
        n_stages = n_stages_RK45
        C = C_RK45
        A = A_RK45
        B = B_RK45
        E = E_RK45
        P = P_RK45
    elif rk_method == 2:
        # Explicit Runge-Kutta method of order 8.
        #     This is a Python implementation of "DOP853"
        raise NotImplementedError
    else:
        raise Exception

    # Determine step size
    min_step = 10 * np.abs(np.nextafter(t, direction * np.inf) - t)
    # Check for over/under shoots, otherwise proceed with the provided step size.
    if step_size_abs > max_step:
        step_size_abs = max_step
    elif step_size_abs < min_step:
        step_size_abs = min_step

    # Integration step parameters
    error_exponent = -1 / (error_estimator_order + 1)
    dtype = y.dtype

    # The A, B, and E parameters must be recast in the same dtype as y.
    A_array = np.asarray(A, dtype=dtype)
    B_array = np.asarray(B, dtype=dtype)
    E_array = np.asarray(E, dtype=dtype)

    # Integration step variables
    t_new = 0.
    y_new = y
    new_deriv_val = current_deriv_val
    current_step_size_abs = step_size_abs
    step = current_step_size_abs * direction
    status = 100
    message = 'Step was completed successfully.'
    step_accepted = False
    step_rejected = False

    while not step_accepted:

        if current_step_size_abs < min_step:
            # Step size too small
            status = -1
            message = 'Step Fail: Step size too small.'
            break

        # Advance time with the current step size
        step = current_step_size_abs * direction
        t_new = t + step

        # Check if the time has surpassed the end of the desired domain
        if direction * (t_new - tf) > 0:
            # If it has, correct the new time to be the last time and update the step size accordingly
            t_new = tf
            step = t_new - t
            current_step_size_abs = np.abs(step)

        # Start RK Specific step calculations
        K[0] = current_deriv_val
        s = 1
        for a, c in zip(A_array[1:], C[1:]):
            dy = np.dot(K[:s].T, a[:s]) * step
            K[s] = func(t + c * step, y + dy, *args)
            s += 1

        y_new = y + step * np.dot(K[:-1].T, B_array)
        new_deriv_val = func(t + step, y_new, *args)

        K[-1] = new_deriv_val
        # End RK Specific step calculations

        scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
        error_norm = error_func_norm(K, E_array, step, scale)

        if error_norm < 1.:
            if error_norm == 0.:
                factor = MAX_FACTOR
            else:
                factor = min(MAX_FACTOR, SAFETY * error_norm**error_exponent)

            if step_rejected:
                factor = min(1., factor)

            current_step_size_abs *= factor

            step_accepted = True
            status = 0
        else:
            current_step_size_abs *= max(MIN_FACTOR, SAFETY * error_norm**error_exponent)
            step_rejected = True

    previous_step = step
    new_step_size_abs = current_step_size_abs
    y_old = y
    t_old = t

    return t_new, y_new, new_deriv_val, new_step_size_abs, y_old, t_old, previous_step, status, message
