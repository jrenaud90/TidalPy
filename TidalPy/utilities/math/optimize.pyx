# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from libc.math cimport fabs
from libcpp.string cimport string as cpp_string
from TidalPy.constants cimport d_NAN_DBL, d_EPS_DBL

cdef struct BracketResult:
    double xa
    double xb
    double fa
    double fb
    size_t iters_p1
    size_t iters_p2

ctypedef double (*BracketFunc)(double x, char* param_ptr) noexcept nogil

cdef BracketResult bracket(
        BracketFunc fn,
        double x0,
        double dx, char* param_ptr, double ratio, size_t maxiter, int& error_state, cpp_string& message) noexcept nogil:
    """
    Given a function and a starting guess, find two
    inputs for the function that bracket a root.

    Parameters
    ----------
    fn : callable
        The function to bracket.
    x0 : double
        The starting guess
    dx : double
        Small step for starting the search
    param_ptr : char*
        Pointer to a parameter struct of any additional parameters required by fn.
    ratio : double
        The step size increases by this ratio every step in the search.
    maxiter : size_t
        The maximum number of steps before giving up.

    Return
    ------
    result : BracketResult
        Struct that contains:
        xa : double
            Left bracket of a root of fn.
        xb : double
            Right bracket of a root of fn.
        fa : double
            Value of fn at xa.
        fb : double
            Value of fn at xb.
    """
    dx = fabs(dx)

    cdef BracketResult result
    result.xa = d_NAN_DBL
    result.xb = d_NAN_DBL
    result.fa = d_NAN_DBL
    result.fb = d_NAN_DBL
    result.iters_p1 = 0
    result.iters_p2 = 0

    cdef const char* message_ptr = "No issues detected."
    if ratio <= 1.0:
        message_ptr = "ERROR: `bracket` - ratio must be greater than 1."
        message = message_ptr
        error_state = -1
        return result
        
    # Get the starting positions
    cdef double f0 = fn(x0, param_ptr)
    cdef double x_left = x0 - dx
    cdef double x_right = x0 + dx
    cdef double f_left = fn(x_left, param_ptr)
    cdef double f_right = fn(x_right, param_ptr)

    # Overshot zero, try making dx smaller
    if (f0 - f_left) * (f_right - f0) < 0.0:
        while (
            (f0 - f_left) * (f_right - f0) < 0.0
            and dx > d_EPS_DBL
            and result.iters_p1 < maxiter
        ):
            dx /= ratio
            x_left = x0 - dx
            x_right = x0 + dx
            f_left = fn(x_left, param_ptr)
            f_right = fn(x_right, param_ptr)
            result.iters_p1 += 1
        if (result.iters_p1 == maxiter):
            # Couldn't find something with same slope in both directions
            message_ptr = "ERROR: `bracket` - Cannot find zero; Maxiter's hit (part 1)."
            message = message_ptr
            error_state = -2
            return result

    cdef double slope = f_right - f0
    cdef double x1 = d_NAN_DBL
    cdef double f1 = d_NAN_DBL
    if slope > 0.0 and f0 > 0.0:  # Walk left
        dx = -dx
        x1 = x_left
        f1 = f_left
    elif slope > 0.0 and f0 < 0.0:  # Walk right
        x1 = x_right
        f1 = f_right
    elif slope < 0.0 and f0 > 0:  # Walk right
        x1 = x_right
        f1 = f_right
    else:  # Walk left
        dx = -dx
        x1 = x_left
        f1 = f_left

    # Do the walking
    cdef double xnew = d_NAN_DBL
    cdef double fnew = d_NAN_DBL
    while f0 * f1 > 0.0 and result.iters_p2 < maxiter:
        dx *= ratio
        xnew = x1 + dx
        fnew = fn(xnew, param_ptr)
        x0 = x1
        f0 = f1
        x1 = xnew
        f1 = fnew
        result.iters_p2 += 1

        if (result.iters_p2 == maxiter):
            # Couldn't find something with same slope in both directions
            message_ptr = "ERROR: `bracket` - Cannot find zero; Maxiter's hit (part 2)."
            message = message_ptr
            error_state = -2
            return result

    if f0 * f1 > 0.0:
        message_ptr = "ERROR: `bracket` - Cannot find zero; f0 * f1 > 0."
        message = message_ptr
        error_state = -3
        return result
    else:
        result.xa = x0
        result.xb = x1
        result.fa = f0
        result.fb = f1
        return result


cdef double test_func(double x, char* param_ptr) noexcept nogil:
    return x**5 - x**3


def test_bracket(double x0, double dx, double ratio=1.618, size_t maxiter=100):

    cdef const char* message_ptr = "No issues detected."
    cdef cpp_string message = cpp_string(message_ptr)
    cdef int error_code = 0

    cdef char* param_ptr = NULL
    cdef BracketResult result = bracket(test_func, x0, dx, param_ptr, ratio, maxiter, error_code, message)

    return result.xa, result.xb, result.fa, result.fb, result.iters_p1, result.iters_p2, error_code, str(message)