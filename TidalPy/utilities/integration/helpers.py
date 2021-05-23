import numpy as np

from ..performance import njit


@njit(cacheable=True)
def norm(x: np.ndarray) -> np.ndarray:
    """Compute RMS norm."""
    return np.linalg.norm(x) / x.size**0.5


@njit(cacheable=True)
def select_initial_step(func: callable, t0: float, y0: np.ndarray, f0: np.ndarray, direction: float, order: float,
                        rtol: float, atol: float) -> float:
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
    f1 = func(t0 + h0 * direction, y1)
    d2 = norm((f1 - f0) / scale) / h0

    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2))**(1 / (order + 1))

    return min(100 * h0, h1)
